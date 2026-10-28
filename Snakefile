# Snakefile for RNA-seq analysis pipeline
#
# A referencia NAO e construida aqui. Construa uma vez com:
#     bash scripts/build_reference.sh
# e aponte config["reference"]["dir"] para onde ela vive (default /mnt/reference/rnaseq).

import glob
import os

from snakemake.utils import min_version
from snakemake.exceptions import WorkflowError

min_version("9.0")

configfile: "config.yaml"
include: "common.smk"

# A referencia e um insumo externo e somente-leitura. Falhar no parse, com
# mensagem clara, e melhor do que comecar a rodar e quebrar na metade.
# A sentinela carrega o sjdbOverhang, entao isso tambem detecta divergencia
# entre o indice construido e o config["star"]["sjdbOverhang"] atual.
if not config["reference"].get("skip_check", False) and not os.path.exists(REFP["done"]):
    raise WorkflowError(
        f"Referencia ausente ou incompleta: nao encontrei {REFP['done']}\n"
        f"  dir           = {REF_ROOT}\n"
        f"  genome_id     = {GENOME}\n"
        f"  overhang      = {OVH_ORIGEM}\n"
        "Construa uma vez com:  bash scripts/build_reference.sh\n"
        "Se as reads deste projeto tem outro comprimento, ajuste star.read_length\n"
        "em config.local.yaml e construa o indice correspondente.\n"
        "Para planejar sem a referencia montada: --config reference='{skip_check: true}'"
    )

SAMPLES = sorted(
    os.path.basename(f).replace("_1.fastq.gz", "")
    for f in glob.glob("samples/*_1.fastq.gz")
)


rule all:
    input:
        expand("output/clean/{sample}_clean_{r}.fastq.gz", sample=SAMPLES, r=[1, 2]),
        expand("output/mapped/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("output/counts/{sample}_ReadsPerGene.out.tab", sample=SAMPLES),
        expand("output/counts/{sample}.genes.results", sample=SAMPLES),
        "output/logs/multiqc_report.html",
        ["output/rdata/DEG.RData"] if config["deg_analysis"]["enabled"] else [],


# Trim adapters and low-quality reads
rule trim_reads:
    input:
        r1 = "samples/{sample}_1.fastq.gz",
        r2 = "samples/{sample}_2.fastq.gz",
    output:
        r1 = "output/clean/{sample}_clean_1.fastq.gz",
        r2 = "output/clean/{sample}_clean_2.fastq.gz",
        html = "output/logs/{sample}.html",
        json = "output/logs/{sample}.json",
    log:
        "output/logs/trim_{sample}.log",
    threads: config["threads"]["trim"]
    resources:
        mem_mb = config["mem_mb"]["trim"],
        runtime = config["runtime_min"]["trim"],
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
        index = REFP["index_flag"],
        gtf = REFP["gtf"],
    output:
        bam = "output/mapped/{sample}_Aligned.sortedByCoord.out.bam",
        transcriptome_bam = "output/mapped/{sample}_Aligned.toTranscriptome.out.bam",
        counts = "output/counts/{sample}_ReadsPerGene.out.tab",
        sj = "output/counts/{sample}_SJ.out.tab",
        final_log = "output/logs/{sample}_Log.final.out",
    log:
        "output/logs/map_{sample}.log",
    threads: config["threads"]["map"]
    resources:
        mem_mb = config["mem_mb"]["map"],
        runtime = config["runtime_min"]["map"],
    params:
        out_prefix = "output/mapped/{sample}_",
        genome_dir = REFP["index_dir"],
        overhang = OVH,
        sort_ram = int(config["star"]["limitBAMsortRAM_gb"]) * 10**9,
    shell:
        """
        ulimit -n 16384

        STAR --genomeDir {params.genome_dir} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --sjdbScore 2 \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --outSAMattributes NH HI AS nM XS \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.out_prefix} \
            --runThreadN {threads} \
            --quantMode TranscriptomeSAM GeneCounts \
            --limitBAMsortRAM {params.sort_ram} \
            --outSAMattrRGline ID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:TRUSEQ\\tPL:ILLUMINA \
            2> {log}

        # O STAR escreve tudo sob {params.out_prefix}, ou seja em output/mapped/.
        # As rules declaram counts/ e logs/, entao movemos como pipeline.sh fazia.
        mv {params.out_prefix}ReadsPerGene.out.tab {output.counts}
        mv {params.out_prefix}SJ.out.tab {output.sj}
        mv {params.out_prefix}Log.final.out {output.final_log}
        mv {params.out_prefix}Log.out output/logs/{wildcards.sample}_Log.out
        mv {params.out_prefix}Log.progress.out output/logs/{wildcards.sample}_Log.progress.out

        rm -rf {params.out_prefix}_STARgenome {params.out_prefix}_STARtmp
        """


# Quantify transcripts with RSEM
rule quantify_transcripts:
    input:
        bam = "output/mapped/{sample}_Aligned.toTranscriptome.out.bam",
        index = REFP["index_flag"],
    output:
        results = "output/counts/{sample}.genes.results",
        cnt = "output/logs/{sample}.cnt",
    log:
        "output/logs/quantify_{sample}.log",
    threads: config["threads"]["quantify"]
    resources:
        mem_mb = config["mem_mb"]["quantify"],
        runtime = config["runtime_min"]["quantify"],
    params:
        index_prefix = REFP["index_prefix"],
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
            {params.index_prefix} \
            output/counts/{wildcards.sample} 2> {log}

        cp output/counts/{wildcards.sample}.stat/{wildcards.sample}.cnt {output.cnt}
        """


# Run MultiQC to generate quality control report
rule run_multiqc:
    input:
        fastp = expand("output/logs/{sample}.json", sample=SAMPLES),
        star = expand("output/logs/{sample}_Log.final.out", sample=SAMPLES),
        rsem = expand("output/logs/{sample}.cnt", sample=SAMPLES),
    output:
        report = "output/logs/multiqc_report.html",
    log:
        "output/logs/multiqc.log",
    resources:
        mem_mb = 4000,
        runtime = 60,
    shell:
        "multiqc output/logs -o output/logs -f 2> {log}"


# Run DEG analysis with R
rule run_deg_analysis:
    input:
        counts = expand("output/counts/{sample}_ReadsPerGene.out.tab", sample=SAMPLES),
        rsem_counts = expand("output/counts/{sample}.genes.results", sample=SAMPLES),
        design = config["deg_analysis"]["experiment_design_path"],
        gene_info = REFP["gene_info"],
    output:
        rdata = "output/rdata/DEG.RData",
    params:
        alpha = config["deg_analysis"]["alpha"],
        lfc_threshold = config["deg_analysis"]["lfc_threshold"],
    log:
        "output/logs/deg_analysis.log",
    threads: config["threads"]["deg"]
    resources:
        mem_mb = config["mem_mb"]["deg"],
        runtime = config["runtime_min"]["deg"],
    shell:
        """
        Rscript {workflow.basedir}/DEG.R \
            {input.design} \
            {params.alpha} \
            {params.lfc_threshold} \
            {input.gene_info} > {log} 2>&1
        """
