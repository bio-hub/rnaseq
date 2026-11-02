#!/usr/bin/env bash
#
# 0.3 — Tabela de traducao de contigs NCBI (RefSeq GRCh38.p14) -> interomics (chr1, chr2, ...).
#
#   sbatch misc/make-contig-map.sbatch
#
# Regra de pareamento, em duas camadas:
#   1. M5 IDENTICO (md5 da sequencia) e 1-para-1: a traducao e por conteudo, nao por nome.
#   2. Contigs sem M5 igual, mas com MESMO COMPRIMENTO (unico entre os sem par) e cuja unica
#      diferenca, base a base, e o interomics ter 'N' onde o NCBI tem base. Descoberto no piloto:
#      o FASTA do interomics e do tipo "analysis set" e mascara com N as falsas duplicacoes do
#      GRCh38 (chr5, chr14, chr19, chr21, chr22) e as PARs do chrY. Os comprimentos sao iguais,
#      entao as coordenadas sao as mesmas; o que muda e a callability (N na referencia).
#      Qualquer outra diferenca (base != base) recusa o par e reprova o job.
# O que nao tem par, por nenhuma das duas, e o que a rna_conform_bam vai descartar.
#
# Saidas em /mnt/shared/projetos/tese_natalia/referencia_ncbi/:
#   rna_contig_map.tsv              2 colunas TAB: nome_NCBI  nome_interomics (sem cabecalho:
#                                   o awk da rna_conform_bam le $1 e $2), na ordem do dict do interomics
#   rna_contig_map.relatorio.txt    contagens, contigs sem par, pares aceitos pela camada 2
#   interomics_mascarado.bed        regioes N do interomics onde o NCBI tem base (nomes do interomics)
#
# Idempotente: se o mapa e o BED existem, sao mais novos que os dois .dict e passam na validacao,
# nao refaz. Temporarios mantem o sufixo final (CLAUDE.md 4.1.1). Le ~500 Mb de cada FASTA:
# roda como job (4.7).

#SBATCH --job-name=contig-map
#SBATCH --partition=gpu,cpu,bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%j.out

set -uo pipefail
export LC_ALL=C                      # CLAUDE.md 4.8: o locale dos nos nao e uniforme
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH
for t in python3 samtools; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERRO: '$t' ausente em $(hostname -s); nada foi verificado" >&2; exit 127; }
done

OUTDIR=/mnt/shared/projetos/tese_natalia/referencia_ncbi
NCBI_DICT=$OUTDIR/GRCh38_p14.dict
NCBI_FA=$OUTDIR/GRCh38_p14.fa
IO_DIR=/mnt/reference/interomics_snakemake/reference_genome
IO_DICT=$IO_DIR/GRCh38_character_genome.dict
IO_FA=$IO_DIR/GRCh38_character_genome.fasta
for f in "$NCBI_DICT" "$IO_DICT" "$NCBI_FA" "$IO_FA" "$NCBI_FA.fai" "$IO_FA.fai"; do
    [ -s "$f" ] || { echo "ERRO: ausente: $f" >&2; exit 2; }
done

python3 - "$NCBI_DICT" "$IO_DICT" "$NCBI_FA" "$IO_FA" "$OUTDIR" <<'PY'
import os, sys, collections, subprocess
import numpy as np
ncbi_dict, io_dict, ncbi_fa, io_fa, outdir = sys.argv[1:6]
MAP = os.path.join(outdir, "rna_contig_map.tsv")
REP = os.path.join(outdir, "rna_contig_map.relatorio.txt")
MASK = os.path.join(outdir, "interomics_mascarado.bed")

def read_dict(p):
    out = []
    for l in open(p):
        if not l.startswith("@SQ"): continue
        f = dict(x.split(":", 1) for x in l.rstrip("\n").split("\t")[1:])
        out.append((f["SN"], int(f["LN"]), f.get("M5", "")))
    return out
N = read_dict(ncbi_dict); I = read_dict(io_dict)
assert all(m for _, _, m in N) and all(m for _, _, m in I), "dict sem M5: nao da para juntar por conteudo"
ncbi_len = {n: l for n, l, _ in N}; io_len = {n: l for n, l, _ in I}

# 25 primarios, escritos a mao e independentes da juncao
vers = [11,12,12,12,10,12,14,11,12,11,10,12,11,9,10,10,11,10,10,11,9,11]
EXPECT = {f"NC_{i:06d}.{v}": f"chr{i}" for i, v in zip(range(1, 23), vers)}
EXPECT.update({"NC_000023.11": "chrX", "NC_000024.10": "chrY", "NC_012920.1": "chrM"})

def validate(rows):
    probs = []
    ks = [a for a, _ in rows]; vs = [b for _, b in rows]
    if len(set(ks)) != len(ks) or len(set(vs)) != len(vs): probs.append("mapa nao e 1-para-1")
    for a, b in rows:
        if a not in ncbi_len or b not in io_len: probs.append(f"nome desconhecido: {a} {b}")
        elif ncbi_len[a] != io_len[b]: probs.append(f"comprimento diverge: {a} {b}")
    got = {a: b for a, b in rows if a.startswith("NC_")}
    for a, b in EXPECT.items():
        if got.get(a) != b: probs.append(f"primario {a}: esperado {b}, mapa diz {got.get(a)}")
    if len(rows) < 25: probs.append(f"menos de 25 linhas ({len(rows)})")
    return probs

# --- idempotencia
newer = max(os.path.getmtime(ncbi_dict), os.path.getmtime(io_dict))
if all(os.path.exists(p) and os.path.getmtime(p) > newer for p in (MAP, MASK, REP)):
    rows = [tuple(l.rstrip("\n").split("\t")) for l in open(MAP)]
    if all(len(r) == 2 for r in rows) and not validate(rows):
        print(f"ja existe e e valido: {MAP} ({len(rows)} linhas)"); sys.exit(0)

# --- camada 1: M5 identico, 1-para-1
by_io = collections.defaultdict(list); by_nc = collections.defaultdict(list)
for n, l, m in I: by_io[m].append(n)
for n, l, m in N: by_nc[m].append(n)
amb_io = sum(1 for v in by_io.values() if len(v) > 1); amb_nc = sum(1 for v in by_nc.values() if len(v) > 1)
rows = [(n, by_io[m][0]) for n, l, m in N if m in by_io and len(by_io[m]) == 1 and len(by_nc[m]) == 1]
n1 = len(rows)

# --- camada 2: mesmo comprimento unico + diferenca so de N-mascaramento no interomics
def seq(fa, name):
    r = subprocess.run(["samtools", "faidx", fa, name], capture_output=True, check=True)
    return np.frombuffer(b"".join(r.stdout.split(b"\n")[1:]).upper(), dtype=np.uint8)
matched_io = {b for _, b in rows}; matched_nc = {a for a, _ in rows}
io_left = [(n, l) for n, l, _ in I if n not in matched_io]
nc_left = [(n, l) for n, l, _ in N if n not in matched_nc]
nc_by_len = collections.defaultdict(list)
for n, l in nc_left: nc_by_len[l].append(n)
io_by_len = collections.Counter(l for _, l in io_left)
accepted = []; refused = []; mask_rows = []
NB = ord("N")
for n, l in io_left:
    cands = nc_by_len.get(l, [])
    if len(cands) != 1 or io_by_len[l] != 1:
        continue
    a = cands[0]
    x = seq(ncbi_fa, a); y = seq(io_fa, n)
    d = x != y
    so_io_n = d & (y == NB) & (x != NB)
    outra = d & ~so_io_n
    if outra.any():
        refused.append((a, n, int(outra.sum()))); continue
    rows.append((a, n)); accepted.append((a, n, int(so_io_n.sum())))
    idx = np.flatnonzero(so_io_n)
    if len(idx):
        brk = np.flatnonzero(np.diff(idx) > 1)
        st = np.concatenate(([idx[0]], idx[brk + 1])); en = np.concatenate((idx[brk], [idx[-1]])) + 1
        mask_rows += [(n, int(s), int(e)) for s, e in zip(st, en)]

pos = {n: i for i, (n, _, _) in enumerate(I)}
rows.sort(key=lambda r: pos[r[1]]); mask_rows.sort(key=lambda r: (pos[r[0]], r[1]))
probs = validate(rows)
mio = {b for _, b in rows}; mnc = {a for a, _ in rows}
sem_io = [(n, l) for n, l, _ in I if n not in mio]
sem_nc = [(n, l) for n, l, _ in N if n not in mnc]
cls = collections.Counter(n[:2] for n, _ in sem_nc)

L = []
L.append(f"dict NCBI: {len(N)} contigs   dict interomics: {len(I)} contigs")
L.append(f"CONTIGS NO MAPA: {len(rows)}  (camada 1, M5 identico: {n1};  camada 2, so N-mascaramento: {len(accepted)})")
L.append(f"  primarios NC_ (esperado 25): {sum(1 for a,_ in rows if a.startswith('NC_'))}")
L.append(f"  nao-NC_ (unlocalized/unplaced NCBI == *_random/chrUn do interomics): {sum(1 for a,_ in rows if not a.startswith('NC_'))}")
L.append(f"M5 repetido dentro de um dict: interomics={amb_io} ncbi={amb_nc}")
L.append("")
L.append("PARES ACEITOS PELA CAMADA 2 (M5 difere; unica diferenca = interomics tem N onde o NCBI tem base):")
for a, n, k in accepted: L.append(f"  {a}\t{n}\tbases_mascaradas={k}")
L.append(f"  total mascarado no interomics nesses contigs: {sum(k for _,_,k in accepted)} bp em {len(mask_rows)} regioes (-> interomics_mascarado.bed)")
if refused:
    L.append("PARES RECUSADOS (diferenca que NAO e N-mascaramento):")
    for a, n, k in refused: L.append(f"  {a}\t{n}\tdiferencas_reais={k}")
L.append("")
L.append(f"CONTIGS DO INTEROMICS SEM PAR: {len(sem_io)}")
for n, l in sem_io: L.append(f"  {n}\t{l}")
L.append("")
L.append(f"CONTIGS DO NCBI SEM PAR (descartados pela rna_conform_bam): {len(sem_nc)}  por prefixo {dict(cls)}  bp={sum(l for _, l in sem_nc)}")
L.append("")
L.append("VALIDACAO: " + ("OK" if not probs and not refused else "FALHOU: " + "; ".join(probs) + (f"; recusados={len(refused)}" if refused else "")))
text = "\n".join(L) + "\n"
print(text)
if probs or refused:
    print("ERRO: validacao falhou; nada foi gravado", file=sys.stderr); sys.exit(1)

def write_atomic(path, content):
    tmp = os.path.join(os.path.dirname(path), "tmp." + os.path.basename(path))   # mesmo sufixo do final
    open(tmp, "w").write(content); os.replace(tmp, path)
write_atomic(MAP, "".join(f"{a}\t{b}\n" for a, b in rows))
write_atomic(MASK, "".join(f"{c}\t{s}\t{e}\n" for c, s, e in mask_rows))
write_atomic(REP, text)
print(f"gravado: {MAP} ({len(rows)} linhas), {MASK} ({len(mask_rows)} regioes)")
PY
rc=$?
[ $rc -eq 0 ] || { echo "ERRO: falhou (rc=$rc)" >&2; exit $rc; }
M="$OUTDIR/rna_contig_map.tsv"
[ -s "$M" ] || { echo "ERRO: mapa ausente apos o job" >&2; exit 1; }
echo "linhas: $(wc -l < "$M")   colunas distintas: $(awk -F'\t' '{print NF}' "$M" | sort -u | tr '\n' ' ')   primarios NC_: $(grep -c '^NC_' "$M")"
