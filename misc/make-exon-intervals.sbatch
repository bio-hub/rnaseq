#!/usr/bin/env bash
#
# 0.4 — Intervalos de exons para o variant calling de RNA-seq (BED + .interval_list).
#
#   sbatch --dependency=afterok:<JOB_DO_0.3> misc/make-exon-intervals.sbatch
#
# Do GTF que o STAR usou: exons dos contigs PRIMARIOS (NC_), renomeados pelo mapa do 0.3,
# com slop de 50 pb (limitado ao comprimento do contig), ordenados na ordem do dict do interomics
# e fundidos. O .interval_list e gerado com BedToIntervalList contra o dict do INTEROMICS
# (namespace de destino do pipeline).
#
# Saidas em /mnt/shared/projetos/tese_natalia/referencia_ncbi/intervalos/ :
#   exons_pad50.bed.gz (+ .tbi)   alvo do variant calling (BED 0-based, bgzip+tabix)
#   exons_pad50.interval_list     o mesmo, formato Picard/GATK (usar em refgenome.interval)
#   exons_nopad.bed               exons fundidos SEM slop, nomes do interomics
#   exons_resumo.txt              tamanhos, intervalos por contig, checagens, genes Tier 1
#
# Verificacoes (o job falha se alguma reprovar; nao se confia so no codigo de saida, CLAUDE.md 4.2.1):
#   - exons sem slop, no namespace NCBI, == 146.738.699 pb em 302.874 intervalos (analise de juncoes)
#   - BED relido do disco: ordenado, sem sobreposicao, dentro dos contigs, contigs no dict
#   - .interval_list: 195 @SQ, n de intervalos e bases iguais ao BED
#   - os 10 genes Tier 1: 100% das bases de exon dentro do alvo (BLOQUEANTE se nao)
# Idempotente: se tudo existe, e mais novo que as entradas e revalida, nao refaz.
# Escreve numa pasta de staging com os nomes FINAIS (extensoes intactas, CLAUDE.md 4.1.1) e
# so move para o destino depois de validar. Le o GTF (2 GB) uma vez: um leitor so (CLAUDE.md 5.3).

#SBATCH --job-name=exon-intervals
#SBATCH --partition=gpu,cpu,bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%j.out

set -uo pipefail
export LC_ALL=C                      # CLAUDE.md 4.8
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH

GTF=/mnt/reference/rnaseq/GRCh38_p14/gtf/GRCh38_p14.gtf
REFDIR=/mnt/shared/projetos/tese_natalia/referencia_ncbi
MAP=$REFDIR/rna_contig_map.tsv
MASK=$REFDIR/interomics_mascarado.bed
IO_DICT=/mnt/reference/interomics_snakemake/reference_genome/GRCh38_character_genome.dict
DEST=$REFDIR/intervalos
STAGE=$DEST/.stage_${SLURM_JOB_ID:-manual}

GATK="${GATK:-}"
if [ -z "$GATK" ]; then
    for g in /mnt/shared/conda/snakemake_envs/interomics/*/bin/gatk; do [ -x "$g" ] && { GATK="$g"; break; }; done
fi
for t in python3 bgzip tabix; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERRO: '$t' ausente em $(hostname -s); nada foi verificado" >&2; exit 127; }
done
[ -n "$GATK" ] && [ -x "$GATK" ] || { echo "ERRO: gatk nao encontrado (defina GATK=...)" >&2; exit 127; }
export PATH="$(dirname "$GATK"):$PATH"          # o wrapper gatk precisa do java do proprio ambiente
command -v java >/dev/null 2>&1 || { echo "ERRO: java ausente apos incluir $(dirname "$GATK")" >&2; exit 127; }
for f in "$GTF" "$MAP" "$MASK" "$IO_DICT"; do [ -s "$f" ] || { echo "ERRO: ausente: $f (o 0.3 rodou?)" >&2; exit 2; }; done
mkdir -p "$DEST" "$STAGE" || exit 1
trap 'rm -rf "$STAGE"' EXIT

# --- validador (usado tambem na idempotencia): relê os arquivos do disco ----------------
cat > "$STAGE/validar.py" <<'PY'
import sys, os, gzip, subprocess
d, bed_path, il_path, nopad_path = sys.argv[1:5]       # d = dict do interomics
dic = []
for l in open(d):
    if l.startswith("@SQ"):
        f = dict(x.split(":", 1) for x in l.rstrip("\n").split("\t")[1:]); dic.append((f["SN"], int(f["LN"])))
ln = dict(dic); order = {n: i for i, (n, _) in enumerate(dic)}
probs = []
def check_bed(p, gz):
    n = bp = 0; prev = None
    for l in (gzip.open(p, "rt") if gz else open(p)):
        c, s, e = l.split("\t")[:3]; s = int(s); e = int(e); n += 1; bp += e - s
        if c not in ln: probs.append(f"{p}: contig fora do dict: {c}"); break
        if not (0 <= s < e <= ln[c]): probs.append(f"{p}: intervalo invalido {c}:{s}-{e}"); break
        if prev:
            pc, pe = prev
            if order[c] < order[pc] or (c == pc and s < pe):
                probs.append(f"{p}: desordenado ou sobreposto em {c}:{s}"); break
            if c == pc and s == pe: probs.append(f"{p}: intervalos encostados nao fundidos em {c}:{s}"); break
        prev = (c, e)
    return n, bp
nb, bb = check_bed(bed_path, True); nn, bn = check_bed(nopad_path, False)
# interval_list: cabecalho @SQ, n de intervalos e bases (1-based inclusivo)
nsq = ni = ib = 0
for l in open(il_path):
    if l.startswith("@SQ"): nsq += 1
    elif l.startswith("@"): continue
    else:
        f = l.split("\t"); ni += 1; ib += int(f[2]) - int(f[1]) + 1
if nsq != len(dic): probs.append(f"interval_list: {nsq} @SQ, esperado {len(dic)}")
if ni != nb: probs.append(f"interval_list: {ni} intervalos, BED tem {nb}")
if ib != bb: probs.append(f"interval_list: {ib} bases, BED tem {bb}")
# tabix abre e responde
r = subprocess.run(["tabix", bed_path, "chr1:1-2000000"], capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip(): probs.append("tabix nao respondeu em chr1:1-2000000")
print(f"VALIDACAO: bed_pad={nb} intervalos/{bb} pb ; bed_nopad={nn}/{bn} ; interval_list={ni}/{ib}/{nsq}@SQ ; " + ("OK" if not probs else "FALHOU: " + "; ".join(probs)))
sys.exit(1 if probs else 0)
PY

OUT_BED=$DEST/exons_pad50.bed.gz; OUT_IL=$DEST/exons_pad50.interval_list; OUT_NP=$DEST/exons_nopad.bed; OUT_RS=$DEST/exons_resumo.txt
if [ -s "$OUT_BED" ] && [ -s "$OUT_BED.tbi" ] && [ -s "$OUT_IL" ] && [ -s "$OUT_NP" ] && [ -s "$OUT_RS" ] \
   && [ "$OUT_BED" -nt "$MAP" ] && [ "$OUT_BED" -nt "$GTF" ] && [ "$OUT_BED" -nt "$IO_DICT" ] \
   && python3 "$STAGE/validar.py" "$IO_DICT" "$OUT_BED" "$OUT_IL" "$OUT_NP"; then
    echo "ja existe e e valido: $DEST"; exit 0
fi

# --- 1. construir o alvo ----------------------------------------------------------------
python3 - "$GTF" "$MAP" "$MASK" "$IO_DICT" "$STAGE" <<'PY'
import sys, re, collections, bisect
gtf, mapf, maskf, dictf, stage = sys.argv[1:6]
FLANK = 50
TIER1 = ["NCSTN","PSENEN","PSEN1","PSEN2","APH1A","APH1B","PSTPIP1","NOD2","MEFV","SREBF1"]
REF_BP, REF_N = 146738699, 302874          # da analise de juncoes (exons NC_ fundidos, namespace NCBI)

cmap = dict(l.rstrip("\n").split("\t") for l in open(mapf))
dic = []
for l in open(dictf):
    if l.startswith("@SQ"):
        f = dict(x.split(":", 1) for x in l.rstrip("\n").split("\t")[1:]); dic.append((f["SN"], int(f["LN"])))
ln = dict(dic); order = {n: i for i, (n, _) in enumerate(dic)}

def merge(iv):                              # iv: lista de (s,e) ordenavel; funde sobrepostos e encostados
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            if e > out[-1][1]: out[-1][1] = e
        else: out.append([s, e])
    return out

exons = collections.defaultdict(list); gene_exons = {g: collections.defaultdict(list) for g in TIER1}
gene_span = collections.defaultdict(list)
rx = re.compile(r'gene_id "([^"]+)"'); nex = 0
t1 = set(TIER1)
for l in open(gtf):
    if l.startswith("#"): continue
    f = l.split("\t", 9)
    c = f[0]
    if not c.startswith("NC_"): continue
    ft = f[2]
    if ft != "exon" and ft != "gene": continue
    s = int(f[3]) - 1; e = int(f[4])
    g = None
    m = rx.search(f[8]); g = m.group(1) if m else None
    if ft == "exon":
        exons[c].append((s, e)); nex += 1
        if g in t1: gene_exons[g][c].append((s, e))
    elif g in t1:
        gene_span[g].append((c, s, e))

# --- cruzamento com a analise de juncoes: exons sem slop, namespace NCBI
tot = sum(e - s for c in exons for s, e in merge(exons[c])); n_iv = sum(len(merge(exons[c])) for c in exons)
cross_ok = (tot == REF_BP and n_iv == REF_N)
print(f"exons NC_ no GTF: {nex}; fundidos sem slop (namespace NCBI): {n_iv} intervalos, {tot} pb  "
      f"(esperado {REF_N} / {REF_BP}) -> {'CONFERE' if cross_ok else 'DIVERGE'}")

miss = [c for c in exons if c not in cmap]
if miss: print("ERRO: contigs NC_ do GTF sem entrada no mapa:", miss, file=sys.stderr); sys.exit(1)

# --- renomear, slop (limitado ao contig), ordenar na ordem do interomics, fundir
nopad = []; pad = []; clamped = 0
for c in sorted(exons, key=lambda k: order[cmap[k]]):
    ic = cmap[c]
    for s, e in merge(exons[c]): nopad.append((ic, s, e))
    raw = []
    for s, e in exons[c]:
        s2 = max(0, s - FLANK); e2 = min(ln[ic], e + FLANK)
        if s2 != s - FLANK or e2 != e + FLANK: clamped += 1
        raw.append((s2, e2))
    for s, e in merge(raw): pad.append((ic, s, e))
with open(f"{stage}/exons_nopad.bed", "w") as o:
    for c, s, e in nopad: o.write(f"{c}\t{s}\t{e}\n")
with open(f"{stage}/exons_pad50.bed", "w") as o:
    for c, s, e in pad: o.write(f"{c}\t{s}\t{e}\n")
pad_bp = sum(e - s for _, s, e in pad)

# --- exons dentro de regioes mascaradas (N) do interomics
mk = collections.defaultdict(list)
for l in open(maskf):
    c, s, e = l.split()[:3]; mk[c].append((int(s), int(e)))
def overlap_bp(iv_by_c, mask_by_c):
    tot = 0
    for c, ivs in iv_by_c.items():
        ms = sorted(mask_by_c.get(c, []))
        if not ms: continue
        starts = [s for s, _ in ms]
        for s, e in ivs:
            i = max(0, bisect.bisect_right(starts, e) - 1)
            while i >= 0 and ms[i][1] > s:
                tot += max(0, min(e, ms[i][1]) - max(s, ms[i][0])); i -= 1
    return tot
np_by_c = collections.defaultdict(list)
for c, s, e in nopad: np_by_c[c].append((s, e))
masked_exon_bp = overlap_bp(np_by_c, mk)

# --- por contig
per = collections.OrderedDict()
for c, s, e in pad: per.setdefault(c, [0, 0, 0, 0]); per[c][0] += 1; per[c][1] += e - s
for c, s, e in nopad: per.setdefault(c, [0, 0, 0, 0]); per[c][2] += 1; per[c][3] += e - s

# --- genes Tier 1: 100% das bases de exon dentro do alvo
pad_by_c = collections.defaultdict(list)
for c, s, e in pad: pad_by_c[c].append((s, e))
def covered(c, ivs):                       # bases de ivs (fundidos) dentro do alvo
    t = pad_by_c.get(c, []); starts = [s for s, _ in t]; got = 0
    for s, e in ivs:
        i = max(0, bisect.bisect_right(starts, s) - 1)
        while i < len(t) and t[i][0] < e:
            got += max(0, min(e, t[i][1]) - max(s, t[i][0])); i += 1
    return got
tier_rows = []; tier_fail = []
for g in TIER1:
    spans = gene_span.get(g, [])
    if len(spans) != 1: tier_fail.append(f"{g}: {len(spans)} features 'gene' em NC_ (esperado 1)"); tier_rows.append((g, "?", 0, 0, 0)); continue
    gc, gs, ge = spans[0]; ic = cmap[gc]
    eb = 0; cv = 0
    for c, ivs in gene_exons[g].items():
        m = merge(ivs); eb += sum(e - s for s, e in m)
        cv += covered(ic, [(s, e) for s, e in m]) if c == gc else 0
    n_ex = sum(len(merge(v)) for v in gene_exons[g].values())
    tier_rows.append((g, f"{ic}:{gs}-{ge}", n_ex, eb, cv))
    if eb == 0 or cv != eb: tier_fail.append(f"{g}: exons={eb} pb, dentro do alvo={cv} pb")

L = []
L.append(f"GTF: {gtf}")
L.append(f"exons NC_ lidos: {nex}")
L.append(f"CRUZAMENTO (sem slop, namespace NCBI): {n_iv} intervalos / {tot} pb  vs esperado {REF_N} / {REF_BP}  -> {'CONFERE' if cross_ok else 'DIVERGE'}")
L.append("")
L.append(f"ALVO FINAL (slop {FLANK} pb, fundido, namespace interomics): {len(pad)} intervalos, {pad_bp} pb  (slop limitado ao contig em {clamped} exons)")
L.append(f"exons sem slop (interomics): {len(nopad)} intervalos, {sum(e-s for _,s,e in nopad)} pb")
L.append(f"bases de exon dentro de regioes N-mascaradas do interomics: {masked_exon_bp} pb")
L.append("")
L.append("INTERVALOS POR CONTIG (alvo final | sem slop):")
L.append("contig\tintervalos_pad\tpb_pad\tintervalos_nopad\tpb_nopad")
for c, (a, b, cc, d) in per.items(): L.append(f"{c}\t{a}\t{b}\t{cc}\t{d}")
L.append("")
L.append("GENES TIER 1 (bases de exon dentro do alvo):")
L.append("gene\tlocus(interomics)\texons\tpb_exon\tpb_dentro\tstatus")
for g, loc, n_ex, eb, cv in tier_rows: L.append(f"{g}\t{loc}\t{n_ex}\t{eb}\t{cv}\t{'OK' if eb and eb == cv else 'FORA'}")
L.append("")
L.append("TIER1: " + ("TODOS DENTRO DO ALVO" if not tier_fail else "BLOQUEANTE: " + "; ".join(tier_fail)))
open(f"{stage}/exons_resumo.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L[:6])); print("\n".join(L[-14:]))
if not cross_ok or tier_fail: sys.exit(3)
PY
rc=$?
[ $rc -eq 0 ] || { echo "ERRO: construcao/checagens reprovaram (rc=$rc); nada foi publicado" >&2; cp -f "$STAGE/exons_resumo.txt" "$DEST/exons_resumo.FALHOU.txt" 2>/dev/null; exit $rc; }

# --- 2. bgzip + tabix e BedToIntervalList (contra o dict do interomics) ----------------------
bgzip -c "$STAGE/exons_pad50.bed" > "$STAGE/exons_pad50.bed.gz" || { echo "ERRO: bgzip" >&2; exit 1; }
tabix -f -p bed "$STAGE/exons_pad50.bed.gz" || { echo "ERRO: tabix" >&2; exit 1; }
"$GATK" --java-options "-Xmx4g" BedToIntervalList -I "$STAGE/exons_pad50.bed" -O "$STAGE/exons_pad50.interval_list" -SD "$IO_DICT" \
    || { echo "ERRO: BedToIntervalList" >&2; exit 1; }

# --- 3. validar o que vai ser publicado, depois publicar -------------------------------------
python3 "$STAGE/validar.py" "$IO_DICT" "$STAGE/exons_pad50.bed.gz" "$STAGE/exons_pad50.interval_list" "$STAGE/exons_nopad.bed" \
    || { echo "ERRO: validacao reprovou; nada foi publicado" >&2; exit 1; }
for f in exons_pad50.bed.gz exons_pad50.bed.gz.tbi exons_pad50.interval_list exons_nopad.bed exons_resumo.txt; do
    mv -f "$STAGE/$f" "$DEST/$f" || { echo "ERRO: mv $f" >&2; exit 1; }
done
rm -f "$DEST/exons_resumo.FALHOU.txt"
python3 "$STAGE/validar.py" "$IO_DICT" "$OUT_BED" "$OUT_IL" "$OUT_NP" || { echo "ERRO: revalidacao do publicado reprovou" >&2; exit 1; }
ls -la "$DEST"
