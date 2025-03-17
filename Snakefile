# Snakefile for RNA-seq analysis pipeline

import os
import glob
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("7.0")

# Configuration
configfile: "config.yaml"

# Define default resource allocation
ncores = os.cpu_count()
memory_mb = int(os.popen("free -m | awk -F' ' 'NR>1 {print $2}' | sed -n 1p").read().strip())
memory_bytes = memory_mb * 1024 * 1024

# Get sample names from the samples directory
SAMPLES = [os.path.basename(f).replace("_1.fastq.gz", "") 
          for f in glob.glob("samples/*_1.fastq.gz")]

# Target rule to define the output files we want
rule all:
    input:
        # Trimmed reads
        expand("output/clean/{sample}_clean_1.fastq.gz", sample=SAMPLES),
        expand("output/clean/{sample}_clean_2.fastq.gz", sample=SAMPLES),
        # Mapped reads
        expand("output/mapped/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        # Count files
        expand("output/counts/{sample}_ReadsPerGene.out.tab", sample=SAMPLES),
        expand("output/counts/{sample}.genes.results", sample=SAMPLES),
        # MultiQC report
        "output/logs/multiqc_report.html",
        # DEG analysis results
        # "output/rdata/DEG.RData"

# Create directories
rule create_directories:
    output:
        sratookit_cache_dir = directory("sratookit_cache"),
        reference_dir = directory("reference"),
        rsem_index_dir = directory("reference/GRCh38_p14_RSEM_index"),
        clean_dir = directory("output/clean"),
        mapped_dir = directory("output/mapped"),
        counts_dir = directory("output/counts"),
        logs_dir = directory("output/logs"),
        rdata_dir = directory("output/rdata")
    shell:
        """
        mkdir -p {output.sratookit_cache_dir}
        mkdir -p {output.reference_dir}
        mkdir -p {output.rsem_index_dir}
        mkdir -p {output.clean_dir}
        mkdir -p {output.mapped_dir}
        mkdir -p {output.counts_dir}
        mkdir -p {output.logs_dir}
        mkdir -p {output.rdata_dir}
        """

# Download reference genome and gene annotation
rule download_reference:
    output:
        genome = "reference/GRCh38_p14.fa",
        annotation = "reference/GRCh38_p14.gtf",
        ncbi_gene_info = "reference/ncbi_gene_info",
        human_gene_info = "reference/GRCh38_p14_gene_info.tsv"
    log:
        "output/logs/download_reference.log"
    shell:
        """
        wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
        --output-document ./reference/GRCh38_p14.fa.gz

        wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
        --output-document ./reference/GRCh38_p14.gtf.gz
        
        wget -c "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz" \
        --output-document ./reference/ncbi_gene_info.gz

        gzip -d ./reference/GRCh38_p14.fa.gz
        gzip -d ./reference/GRCh38_p14.gtf.gz
        gzip -d ./reference/ncbi_gene_info.gz

        Rscript ./ETL_gene_info.R > {log} 2>&1
        """

# Create RSEM index
rule build_rsem_index:
    input:
        genome = "reference/GRCh38_p14.fna",
        annotation = "reference/GRCh38_p14.gtf"
    output:
        index = "reference/GRCh38_p14_RSEM_index/genomeParameters.txt"
    log:
        "output/logs/build_rsem_index.log"
    threads: lambda wildcards, attempt: min(ncores - 2, 20)
    shell:
        """
        rsem-prepare-reference -p {threads} \
            --gtf {input.annotation} \
            --star \
            {input.genome} \
            ./reference/GRCh38_p14_RSEM_index/GRCh38_p14 2> {log}
        """

# Trim adapters and low-quality reads
rule trim_reads:
    input:
        r1 = "samples/{sample}_1.fastq.gz",
        r2 = "samples/{sample}_2.fastq.gz"
    output:
        r1 = "output/clean/{sample}_clean_1.fastq.gz",
        r2 = "output/clean/{sample}_clean_2.fastq.gz",
        html = "output/logs/{sample}.html",
        json = "output/logs/{sample}.json"
    log:
        "output/logs/trim_{sample}.log"
    threads: lambda wildcards, attempt: min(ncores - 2, 8)
    shell:
        """
        fastp --thread {threads} \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --html {output.html} \
            --json {output.json} 2> {log}
        """

# Map reads to the reference genome
rule map_reads:
    input:
        r1 = "output/clean/{sample}_clean_1.fastq.gz",
        r2 = "output/clean/{sample}_clean_2.fastq.gz",
        index = "reference/GRCh38_p14_RSEM_index/genomeParameters.txt",
        gtf = "reference/GRCh38_p14.gtf"
    output:
        bam = "output/mapped/{sample}_Aligned.sortedByCoord.out.bam",
        transcriptome_bam = "output/mapped/{sample}_Aligned.toTranscriptome.out.bam",
        counts = "output/counts/{sample}_ReadsPerGene.out.tab",
        sj = "output/counts/{sample}_SJ.out.tab",
        logs = "output/logs/{sample}_Log.out"
    log:
        "output/logs/map_{sample}.log"
    threads: lambda wildcards, attempt: min(ncores - 2, 16)
    params:
        out_prefix = "output/mapped/{sample}_",
        memory_bytes = lambda wildcards, attempt: memory_bytes - 8000000000
    shell:
        """
        ulimit -n 16384
        
        STAR --genomeDir ./reference/GRCh38_p14_RSEM_index \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --sjdbScore 2 \
            --outSAMattributes NH HI AS nM XS \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.out_prefix} \
            --runThreadN {threads} \
            --quantMode TranscriptomeSAM GeneCounts \
            --limitGenomeGenerateRAM {params.memory_bytes} \
            --limitBAMsortRAM {params.memory_bytes} \
            --sjdbGTFfile {input.gtf} \
            --outSAMattrRGline ID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:TRUSEQ\\tPL:ILLUMINA \
            2> {log}
            
        # Move files to their correct locations (done by directing output to proper locations)
        """

# Quantify transcripts with RSEM
rule quantify_transcripts:
    input:
        bam = "output/mapped/{sample}_Aligned.toTranscriptome.out.bam",
        index = "reference/GRCh38_p14_RSEM_index/genomeParameters.txt"
    output:
        results = "output/counts/{sample}.genes.results",
        cnt = "output/logs/{sample}.cnt"
    log:
        "output/logs/quantify_{sample}.log"
    threads: lambda wildcards, attempt: min(ncores - 2, 12)
    shell:
        """
        rsem-calculate-expression \
            --seed 1954 \
            -p {threads} \
            --paired-end \
            --alignments \
            --no-bam-output \
            --append-names \
            {input.bam} \
            ./reference/GRCh38_p14_RSEM_index/GRCh38_p14 \
            ./output/counts/{wildcards.sample} 2> {log}
            
        mkdir -p ./output/counts/{wildcards.sample}.stat/
        
        # Copy the count file to logs
        cp ./output/counts/{wildcards.sample}.stat/{wildcards.sample}.cnt {output.cnt}
        """

# Run MultiQC to generate quality control report
rule run_multiqc:
    input:
        expand("output/logs/{sample}.html", sample=SAMPLES),
        expand("output/logs/{sample}.cnt", sample=SAMPLES)
    output:
        report = "output/logs/multiqc_report.html"
    log:
        "output/logs/multiqc.log"
    shell:
        """
        multiqc ./output/logs -o ./output/logs -f 2> {log}
        """

#Run DEG analysis with R
rule run_deg_analysis:
    input:
        counts = expand("output/counts/{sample}_ReadsPerGene.out.tab", sample=SAMPLES),
        rsem_counts = expand("output/counts/{sample}.genes.results", sample=SAMPLES)
    output:
        rdata = "output/rdata/DEG.RData"
    params: 
       exp_design = config["deg_analysis"]["experiment_design_path"]
       alpha = config["deg_analysis"]["alpha"]
       lfc_threshold = config["deg_analysis"]["lfc_threshold"]
    log:
        "output/logs/deg_analysis.log"
    threads: lambda wildcards, attempt: min(ncores - 2, 12)
    shell:
        """
        Rscript ./DEG.R {params.exp_design} {threads} {params.alpha} {params.lfc_threshold} > {log} 2>&1
        """