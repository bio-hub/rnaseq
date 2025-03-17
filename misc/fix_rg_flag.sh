#computational resources
ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)
memory_byte=$(free -b | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
    samtools addreplacerg \
    -r "ID:"$i"\\tSM:"$i"\\tLB:TRUSEQ\\tPL:ILLUMINA" \
    -w  \
    -o $i"_Aligned.sortedByCoord.out_ok.bam" \
    -@ $((ncores - 2))
    $i"_Aligned.sortedByCoord.out.bam"

        samtools addreplacerg \
    -r "ID:"$i"\\tSM:"$i"\\tLB:TRUSEQ\\tPL:ILLUMINA" \
    -w  \
    -o $i"_Aligned.toTranscriptome.out_ok.bam" \
    -@ $((ncores - 2))
    $i"_Aligned.toTranscriptome.out.bam"

done