#!/usr/bin/env bash
#
# Reescreve o read group de BAMs ja alinhados, para corridas antigas em que
# o STAR rodou sem --outSAMattrRGline. O Snakefile atual ja grava o RG no
# alinhamento, entao isto so serve para dados legados.
#
#   bash misc/fix_rg_flag.sh [dir_samples] [dir_mapped]

set -euo pipefail

samples_dir="${1:-./samples}"
mapped_dir="${2:-./output/mapped}"

ncores=$(nproc)
threads=$(( ncores > 2 ? ncores - 2 : 1 ))

for i in $(ls "$samples_dir" | grep _1.fastq.gz | sed 's/_1.fastq.gz//' | uniq); do
    echo "==> ${i}"
    for suffix in Aligned.sortedByCoord.out Aligned.toTranscriptome.out; do
        in_bam="${mapped_dir}/${i}_${suffix}.bam"
        out_bam="${mapped_dir}/${i}_${suffix}_ok.bam"
        [[ -f "$in_bam" ]] || { echo "   pulando, nao existe: $in_bam"; continue; }
        samtools addreplacerg \
            -r "ID:${i}\tSM:${i}\tLB:TRUSEQ\tPL:ILLUMINA" \
            -w \
            -o "$out_bam" \
            -@ "$threads" \
            "$in_bam"
    done
done

echo "Pronto."
