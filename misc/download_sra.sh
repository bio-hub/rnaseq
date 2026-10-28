#!/usr/bin/env bash
#
# Baixa corridas paired-end do SRA listadas na primeira coluna de um CSV.
#
#   bash misc/download_sra.sh <lista.csv> [dir_saida] [dir_temp]
#
# O dir_temp precisa ficar em disco: /tmp neste cluster e tmpfs em RAM
# (ver plano de rede, secao 13.7) e o fasterq-dump estoura a memoria la.

set -euo pipefail

csv="${1:?uso: bash misc/download_sra.sh <lista.csv> [dir_saida] [dir_temp]}"
outdir="${2:-samples}"
tmpdir="${3:-/mnt/scratch-shared/$USER/sra_tmp}"

# acima de ~16 threads o fasterq-dump nao ganha mais nada
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
    pigz -p "$threads" "$outdir/${id}_1.fastq" "$outdir/${id}_2.fastq"
    rm -rf "$outdir/$id"
done

echo "Pronto. FASTQ em ${outdir}/"
