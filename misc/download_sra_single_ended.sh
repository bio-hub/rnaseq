#!/usr/bin/env bash
#
# Baixa corridas single-end do SRA listadas na primeira coluna de um CSV.
#
#   bash misc/download_sra_single_ended.sh <lista.csv> [dir_saida] [dir_temp]
#
# Ver nota sobre /tmp em misc/download_sra.sh.

set -euo pipefail

csv="${1:?uso: bash misc/download_sra_single_ended.sh <lista.csv> [dir_saida] [dir_temp]}"
outdir="${2:-samples}"
tmpdir="${3:-/mnt/scratch-shared/$USER/sra_tmp}"

ncores=$(nproc)
threads=$(( ncores > 16 ? 16 : ncores ))
memory_mega=$(free -m | awk 'NR==2 {print $2}')

mkdir -p "$outdir" "$tmpdir"

samples=$(awk -F"," 'NR>1 {print $1}' "$csv")

for id in ${samples}; do
    echo "==> ${id}"
    prefetch --output-directory "$outdir" "$id"
    fasterq-dump --progress \
        --threads "$threads" \
        --mem "${memory_mega}M" \
        --temp "$tmpdir" \
        --outdir "$outdir" \
        "$outdir/$id"
    pigz -p "$threads" "$outdir/${id}.fastq"
    rm -rf "$outdir/$id"
done

echo "Pronto. FASTQ em ${outdir}/"
