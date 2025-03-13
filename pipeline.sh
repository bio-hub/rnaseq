#computational resources
ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)
memory_byte=$(free -b | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

#creating directories
mkdir reference
mkdir -p ./reference/GRCh38_p14_RSEM_index/    
mkdir -p ./output/clean
mkdir -p ./output/mapped
mkdir -p ./output/counts
mkdir -p ./output/logs

#dowload reference genome and gene annotation
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz --output-document ./reference/GRCh38_p14.fna.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz --output-document ./reference/GRCh38_p14.gtf.gz

gzip -d ./reference/GCF_000001405.40_GRCh38.p14_genomic.fna.gz 
gzip -d ./reference/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

#create index
# STAR --genomeDir ./reference/GRCh38_p14_STAR_index \
# --runMode genomeGenerate \
# --runThreadN $((ncores - 4)) \
# --limitGenomeGenerateRAM $((memory_byte - 8000000000)) \
# --genomeFastaFiles ./reference/GRCh38_p14.fa \
# --sjdbGTFfile ./reference/GRCh38_p14.gtf \
# --sjdbOverhang 99

rsem-prepare-reference -p $((ncores - 2)) \
--gtf ./reference/GRCh38_p14.gtf \
--star \
./reference/GRCh38_p14.fa \
./reference/GRCh38_p14_RSEM_index/GRCh38_p14

#trim adapters
for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
  fastp --thread $((ncores - 2)) \
  -i "./samples/"$i"_1.fastq.gz" \
  -I "./samples/"$i"_2.fastq.gz" \
  -o "./output/clean/"$i"_clean_1.fastq.gz" \
  -O "./output/clean/"$i"_clean_2.fastq.gz" \
  --html "./output/logs/"$i".html" \
  --json "./output/logs/"$i".json"
done

#mapping with the reference
for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
  ulimit -n 16384
  STAR --genomeDir ./reference/GRCh38_p14_RSEM_index \
  --readFilesIn \
  "./output/clean/"$i"_clean_1.fastq.gz" \
  "./output/clean/"$i"_clean_2.fastq.gz" \
  --readFilesCommand zcat \
  --sjdbScore 2 \
  --outSAMattributes NH HI AS nM XS \
  --outFilterIntronMotifs RemoveNoncanonical \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "./output/mapped/"$i"_" \
  --runThreadN $((ncores - 2)) \
  --quantMode TranscriptomeSAM GeneCounts \
  --limitGenomeGenerateRAM $((memory_byte - 8000000000)) \
  --limitBAMsortRAM $((memory_byte - 8000000000)) \
  --sjdbGTFfile ./reference/GRCh38_p14.gtf \
  --outSAMattrRGline "ID:"$id"\\tSM:"$i"\\tLB:TRUSEQ\\tPL:ILLUMINA"
done

for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
${RSEM}/rsem-calculate-expression --seed 1858 -p $nucleo --paired-end \
--alignments \
--no-bam-output \
 --append-names \
$i"_Aligned.toTranscriptome.out.bam" \
${db}/GRCh38_RSEM_index/GRCh38 \
$i
done

#remove tmp STAR folders
rm -r ./output/mapped/*__STARgenome
rm -r ./output/mapped/*__STARtemp

#move files to other folders
mv *.out* ../logs
mv *.stat ../counts

#calculate DEGs
Rscript ./DEG.R

