#!/usr/bin/env bash
#
# Investigacao do 0.3: por que 6 cromossomos primarios (chr5, chr14, chr19, chr21, chr22, chrY)
# NAO tem M5 igual entre o FASTA NCBI (usado pelo STAR) e o do interomics?
#
# Compara base a base, nesses contigs, as duas sequencias e classifica cada diferenca:
#   so_interomics_N : interomics tem N e o NCBI tem base  (mascaramento no interomics)
#   so_ncbi_N       : NCBI tem N e o interomics tem base
#   outra           : as duas tem base e divergem (mudanca REAL de sequencia: bloqueante)
# e grava as regioes mascaradas (BED, coordenadas NCBI == interomics, pois os comprimentos
# sao iguais) para cruzar com os genes Tier 1.
#
# Le ~500 Mb de cada FASTA: por isso roda como job (CLAUDE.md 4.7).
#   sbatch misc/compare-masked-contigs.sbatch
#
# Saidas (scratch duravel, nao toca /mnt/reference nem os results):
#   /mnt/scratch-shared/ronald/bai_piloto/contig_diff/{resumo.txt,masked_regions.bed}

#SBATCH --job-name=contig-diff
#SBATCH --partition=gpu,cpu,bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/scratch-shared/ronald/bai_piloto/slurm/%x_%j.out

set -uo pipefail
export LC_ALL=C
export PATH=/mnt/shared/conda/envs/biohub_rnaseq/bin:$PATH
for t in python3 samtools; do command -v "$t" >/dev/null 2>&1 || { echo "ERRO: '$t' ausente em $(hostname -s)" >&2; exit 127; }; done

NCBI_FA=/mnt/shared/projetos/tese_natalia/referencia_ncbi/GRCh38_p14.fa
IO_FA=/mnt/reference/interomics_snakemake/reference_genome/GRCh38_character_genome.fasta
OUT=/mnt/scratch-shared/ronald/bai_piloto/contig_diff
mkdir -p "$OUT"
[ -s "${IO_FA}.fai" ] || { echo "ERRO: sem .fai do interomics: ${IO_FA}.fai" >&2; exit 2; }

python3 - "$NCBI_FA" "$IO_FA" "$OUT" <<'PY'
import subprocess, sys, numpy as np
ncbi_fa, io_fa, out = sys.argv[1:4]
PAIRS = [("NC_000005.10","chr5"),("NC_000014.9","chr14"),("NC_000019.10","chr19"),
         ("NC_000021.9","chr21"),("NC_000022.11","chr22"),("NC_000024.10","chrY")]
def seq(fa, name):
    r = subprocess.run(["samtools","faidx",fa,name], capture_output=True, check=True)
    body = b"".join(r.stdout.split(b"\n")[1:])
    return np.frombuffer(body.upper(), dtype=np.uint8)
N = ord("N")
lines=[]; bed=open(f"{out}/masked_regions.bed","w")
for a, b in PAIRS:
    x = seq(ncbi_fa, a); y = seq(io_fa, b)
    if len(x) != len(y):
        lines.append(f"{a} vs {b}: COMPRIMENTOS DIFEREM {len(x)} {len(y)}"); continue
    d = x != y
    so_io_n = d & (y == N) & (x != N); so_nc_n = d & (x == N) & (y != N); outra = d & ~so_io_n & ~so_nc_n
    lines.append(f"{a} vs {b}: len={len(x)} diferencas={int(d.sum())}  so_interomics_N={int(so_io_n.sum())}  so_ncbi_N={int(so_nc_n.sum())}  OUTRA={int(outra.sum())}")
    # intervalos mascarados (so_interomics_N), fundidos
    idx = np.flatnonzero(so_io_n)
    if len(idx):
        brk = np.flatnonzero(np.diff(idx) > 1)
        starts = np.concatenate(([idx[0]], idx[brk+1])); ends = np.concatenate((idx[brk], [idx[-1]])) + 1
        for s, e in zip(starts, ends): bed.write(f"{a}\t{s}\t{e}\n")
        lines.append(f"   regioes mascaradas (so_interomics_N): {len(starts)}; maior={int((ends-starts).max())} bp; total={int((ends-starts).sum())} bp")
    if outra.sum():
        j = np.flatnonzero(outra)[:5]
        lines.append(f"   primeiras 5 'outra' (0-based): {j.tolist()} ncbi={[chr(x[k]) for k in j]} interomics={[chr(y[k]) for k in j]}")
bed.close()
open(f"{out}/resumo.txt","w").write("\n".join(lines)+"\n")
print("\n".join(lines))
PY
rc=$?; [ $rc -eq 0 ] || { echo "ERRO: falhou rc=$rc" >&2; exit $rc; }
wc -l "$OUT/masked_regions.bed"
