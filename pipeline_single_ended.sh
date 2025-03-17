#!/bin/bash

#to run this script execute the command line:
#bash your_script.sh /path/to/fastq.gz/samples/folder /path/to/experiment_design.csv

# Get variables from arguments
samples_path=$1
exp_design_path=$2

echo "Samples folder: $samples_path"
echo "Experiment design file: $exp_design_path"

#computational resources
ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)
memory_byte=$(free -b | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

#creating directories
mkdir -p ./reference/GRCh38_p14_RSEM_index/    
mkdir -p ./output/clean
mkdir -p ./output/mapped
mkdir -p ./output/counts
mkdir -p ./output/logs
mkdir -p ./output/rdata

#dowload reference genome and gene annotation
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz --output-document ./reference/GRCh38_p14.fa.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz --output-document ./reference/GRCh38_p14.gtf.gz
wget -c "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz" --output-document ./reference/ncbi_gene_info.gz

gzip -d ./reference/GRCh38_p14.fa.gz 
gzip -d ./reference/GRCh38_p14.gtf.gz

#create index
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

#create gene_info file
echo "creating gene_info file..."
Rscript ./ETL_gene_info.R
echo "Done!"

#trimming adapters
echo "Trimming adapters and low-quality reads..."
for i in $(ls $samples_path | grep .fastq.gz | sed 's/.fastq.gz*//' | uniq); do
  echo $i
  fastp --thread $((ncores - 2)) \
  -i "$samples_path/"$i".fastq.gz" \
  -o "./output/clean/"$i"_clean.fastq.gz" \
  --html "./output/logs/"$i".html" \
  --json "./output/logs/"$i".json"
done
echo "Done!"

#mapping with the reference
echo "Mapping with the reference genome..."
for i in $(ls $samples_path  | grep .fastq.gz | sed 's/.fastq.gz*//' | uniq); do
  echo $i
  ulimit -n 16384
  STAR --genomeDir ./reference/GRCh38_p14_RSEM_index \
  --readFilesIn \
  "./output/clean/"$i"_clean.fastq.gz" \
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
  --outSAMattrRGline "ID:"$i"\\tSM:"$i"\\tLB:TRUSEQ\\tPL:ILLUMINA"

  #remove tmp STAR folders
  rm -r ./output/mapped/*__STARgenome

  #move files to other folders
  mv ./output/mapped/*_ReadsPerGene.out.tab ./output/counts
  mv ./output/mapped/*_SJ.out.tab ./output/counts
  mv ./output/mapped/*.out ./output/logs

done
echo "Done!"

#counting transcripts
echo "Quantifiying transcripts..."
for i in $(ls $samples_path | grep .fastq.gz | sed 's/.fastq.gz*//' | uniq); do
  rsem-calculate-expression \
  --seed 1954 \
  -p $((ncores - 2)) \
  --alignments \
  --no-bam-output \
  --append-names \
  "./output/mapped/"$i"_Aligned.toTranscriptome.out.bam" \
  ./reference/GRCh38_p14_RSEM_index/GRCh38_p14 \
  "./output/counts/"$i

  cp "./output/counts/"$i".stat/"$i".cnt" "./output/logs/"$i".cnt"
done
echo "Done!"

#generating fastq QC report
echo "Generating fastq QC report"
multiqc ./output/logs
mv multiqc_data multiqc_report.html ./output/logs
echo "Done!"

#calculate DEGs
echo "calculating DEGs and Reactome analysis..."
Rscript ./DEG.R $exp_design_path 0.05 1
echo "Done!"

#moving output to samples file
echo "moving output to samples file"
mv ./output/ $samples_path
echo "Done!"
