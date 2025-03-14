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
mkdir -p ./output/rdata

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

echo "Building index files for the reference genome..."
if [ ! -f ./reference/GRCh38_p14_RSEM_index/genomeParameters.txt ]; then
    rsem-prepare-reference -p $((ncores - 2)) \
    --gtf ./reference/GRCh38_p14.gtf \
    --star \
    ./reference/GRCh38_p14.fa \
    ./reference/GRCh38_p14_RSEM_index/GRCh38_p14
    echo "Done!"
else
    echo "reference genome index already exists. Skipping to next step."
fi

#trim adapters
echo "Trimming adapters and low-quality reads..."
for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
  echo $i
  fastp --thread $((ncores - 2)) \
  -i "./samples/"$i"_1.fastq.gz" \
  -I "./samples/"$i"_2.fastq.gz" \
  -o "./output/clean/"$i"_clean_1.fastq.gz" \
  -O "./output/clean/"$i"_clean_2.fastq.gz" \
  --html "./output/logs/"$i".html" \
  --json "./output/logs/"$i".json"
done
echo "Done!"

#mapping with the reference
echo "Mapping with the reference genome..."
for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
  echo $i
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
echo "Done!"

#counting transcripts
echo "Quantifiying transcripts..."
for i in $(ls ./samples | grep _1.fastq.gz | sed 's/_1.fastq.gz*//' | uniq); do
  rsem-calculate-expression \
  --seed 1858 \
  -p $((ncores - 2)) \
  --paired-end \
  --alignments \
  --no-bam-output \
  --append-names \
  "./output/mapped/"$i"_Aligned.toTranscriptome.out.bam" \
  ./reference/GRCh38_RSEM_index/GRCh38 \
  $i
done
echo "Done!"

#remove tmp STAR folders
rm -r ./output/mapped/*__STARgenome

#move files to other folders
mv ./output/mapped/*.out* ./output/logs
mv ./output/mapped/*.stat ./output/counts

#calculate DEGs
echo "calculating DEGs and Reactome analysis..."
#Rscript ./DEG.R
echo "Done!"
