ncores=$(cat /proc/cpuinfo | grep processor | wc -l)
memory_mega=$(free -m | awk -F" " 'NR>1 {print $2}' | sed -n 1p)
memory_byte=$(free -b | awk -F" " 'NR>1 {print $2}' | sed -n 1p)

#dowload reference genome and gene annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz 
gzip -d GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

mv GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38_p14.fa
mv GCF_000001405.40_GRCh38.p14_genomic.gtf GRCh38_p14.gtf


STAR --genomeDir ./reference/GRCh38_p14_STAR_index \
--runMode genomeGenerate \
--runThreadN $((ncores - 2)) \
--limitGenomeGenerateRAM $((memory_byte - 500000000)) \
--genomeFastaFiles ./reference/GRCh38_p14.fa \
--sjdbGTFfile ./reference/GRCh38_p14.gtf \
--sjdbOverhang 99

mkdir ./reference/GRCh38_p14_RSEM_index/    

rsem-prepare-reference -p $((ncores - 2)) \
--gtf ./reference/GRCh38_p14.gtf \
--star \
./reference/GRCh38_p14.fa \
./reference/GRCh38_p14_RSEM_index/GRCh38_p14

###############################
###### Sampling trasncriptome #
###############################
cd ${home_dir}/users/${username}/${projectname}/transcriptome
cd rawdata

#creating symbolic links of the rawdata

nfiles=$(sed -n '$=' ${input_transcriptome})
for i in $(seq 2 ${nfiles})
do
  ln -s $(awk -F"\t" '{print $1}' ${input_transcriptome} | sed -n "$i"p) ./
done

#####dicovering adapter#####
export adapter=$(awk -F"\t" '{print $10}' ${input_transcriptome} | sed -n 2p)

################################
###### fastQC PRE CLEAN ########
################################
#Runing the fastQC in all files (before trimming)

ls | grep .fastq.gz | parallel --max-args=1 \
${FastQC}/fastqc --outdir  ../QC/fastqc_reports_before/www {1} 

### FASTQ_SUMMARY.SH ###

fastqc_report=fastqc_reports_before

cd ../QC/${fastqc_report}/www

chmod +x ${FastQC}/fastaq_summary_before.sh
sh ${FastQC}/fastaq_summary_before.sh

#display the outputs of fastqc in a browser with a different tab for each fastq file
#firefox '--new-window' ${home_dir}/users/${username}/${projectname}/transcriptome/QC/fastqc_reports_before/fastaqc_summary_before.html &

#####ITERATION START
cd ${home_dir}/users/${username}/${projectname}/transcriptome
cd rawdata

if [ $(echo $pair_ended) == true ]
  then
  echo detected pair-ended reads!
  
TEST=$(ls | grep _1.fastq.gz | sed 's/_.*//')
for i in $TEST; 
do

### TRIM ADAPTERS AND QC

java -Xmx64G -jar ${home_dir}/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $nucleo \
-basein $i"_1.fastq.gz" \
-baseout $i".fastq.gz" \
ILLUMINACLIP:${home_dir}/apps/Trimmomatic-0.39/adapters/${adapter}.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

###MAPPING TO GENOME REFERENCE

ulimit -n 16384
${home_dir}/apps/STAR --genomeDir ${db}/GRCh38_RSEM_index \
--readFilesIn \
$i"_1P.fastq.gz" \
$i"_2P.fastq.gz" \
--readFilesCommand zcat \
--sjdbScore 2 \
--outSAMattributes NH HI AS nM XS \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $i"_" \
--runThreadN $nucleo \
--quantMode TranscriptomeSAM GeneCounts \
--limitGenomeGenerateRAM 128000000000 \
--limitBAMsortRAM 128000000000 \
--sjdbGTFfile ${db}/GRCh38_p14.gtf \
--outSAMattrRGline "ID:"$i"\tSM:"$i


${RSEM}/rsem-calculate-expression --seed 1858 -p $nucleo --paired-end \
--alignments \
--no-bam-output \
 --append-names \
$i"_Aligned.toTranscriptome.out.bam" \
${db}/GRCh38_RSEM_index/GRCh38 \
$i

done

else #for single_ended

#####ITERATION START
TEST=$(ls | grep .fastq.gz | sed 's/.fastq.gz*//')
for i in $TEST; 
do

### TRIM ADAPTERS AND QC

java -Xmx64G -jar ${home_dir}/apps/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $nucleo \
$i".fastq.gz" \
$i"_1P.fastq.gz" \
ILLUMINACLIP:${home_dir}/apps/Trimmomatic-0.39/adapters/${adapter}.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

ulimit -n 16384
${home_dir}/apps/STAR --genomeDir ${db}/GRCh38_RSEM_index \
--readFilesIn \
$i"_1P.fastq.gz" \
--readFilesCommand zcat \
--sjdbScore 2 \
--outSAMattributes NH HI AS nM XS \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $i"_" \
--runThreadN $nucleo \
--quantMode TranscriptomeSAM GeneCounts \
--limitGenomeGenerateRAM 128000000000 \
--limitBAMsortRAM 128000000000 \
--sjdbGTFfile ${db}/GRCh38_p14.gtf \
--outSAMattrRGline "ID:"$i" SM:"$i

#"'@RG\\tID:{sample}\\tSM:{sample}\\tLB:SureSelectV7\\tPL:ILLUMINA'"
#"ID:"$i"\\tSM:"$i"\\tLB:SureSelectV7\\tPL:ILLUMINA"


${RSEM}/rsem-calculate-expression --seed 1858 -p $nucleo --paired-end \
--alignments \
--no-bam-output \
 --append-names \
$i"_Aligned.toTranscriptome.out.bam" \
${db}/GRCh38_RSEM_index/GRCh38 \
$i

samtools index $i"_Aligned.sortedByCoord.out.bam"
mv _Aligned.sortedByCoord.out.bam.bai _Aligned.sortedByCoord.out.bai

done
fi

#remove tmp STAR folders
rm -r *__STARgenome

#move files to other folders
mv *.out.bam ../mapped
mv *.out.bai ../mapped
mv *.out* ../transcriptome_report
mv *.results ../mapped_rsem
mv *.stat ../mapped_rsem

################################
###### fastQC POS CLEAN ########
################################
cd ${home_dir}/users/${username}/${projectname}/transcriptome
cd rawdata

#Runing the fastQC in all files (before trimming)
#mkdir fastqc_reports_before

ls | grep .fastq.gz | parallel --max-args=1 \
${FastQC}/fastqc --outdir ../QC/fastqc_reports_after/www/ {1}

### FASTQ_SUMMARY.SH ###

fastqc_report=fastqc_reports_after

cd ../QC/${fastqc_report}/www

chmod +x ${FastQC}/fastaq_summary_after.sh
sh ${FastQC}/fastaq_summary_after.sh

#display the outputs of fastqc in a browser with a different tab for each fastq file
#firefox '--new-window' ${home_dir}/users/${username}/${projectname}/transcriptome/QC/fastqc_reports_after/fastaqc_summary_after.html &

#####################	
#### MERGE LANES #### 
#####################

cd ${home_dir}/users/${username}/${projectname}/transcriptome
cd ./mapped

for i in $(awk -F"\t" '{print $3}'  ${input_transcriptome} | tail -n +2 | uniq);
do
if [ "$(awk -F"\t" '{print $2,$3}'  ${input_transcriptome} | tail -n +2 | uniq | grep $i | awk -F" " '{print $1}' | wc -l)" -ge 2 ]
then
awk -F"\t" '{print $2,$3}'  ${input_transcriptome} | tail -n +2 | uniq | grep $i | awk -F" " '{print $1"_Aligned.sortedByCoord.out.bam"}' > merge_lanes.txt
samtools merge $i"_Aligned.sortedByCoord.out.bam" -b merge_lanes.txt --threads $nucleo
samtools index $i"_Aligned.sortedByCoord.out.bai"
rm -rf `cat merge_lanes.txt`
else 
samtools index $i"_Aligned.sortedByCoord.out.bai"
fi
done

#####################	
#### DEGS ########### 
#####################

cd ${home_dir}/users/${username}/${projectname}/transcriptome
Rscript ${home_dir}/transcriptome/output/DEG.R

