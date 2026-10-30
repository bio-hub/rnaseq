# Snakefile for RNA-seq analysis pipeline
#
# A referencia NAO e construida aqui. Construa uma vez com:
#     bash scripts/build_reference.sh
#
# Intermediarios vao para config["output"]["work_dir"] (scratch) e apenas os
# entregaveis sao exportados para config["output"]["final_dir"].

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
# entre o indice construido e o config["star"]["read_length"] atual.
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

SHARED_GENOME = bool(config["star"].get("shared_genome", False))


# O alvo final e o export: tudo que importa precisa ter chegado ao final_dir.
# Os intermediarios no work_dir sao consequencia, nao objetivo.
rule all:
    input:
        f"{FINAL}/.export_done",


# Trim adapters and low-quality reads
rule trim_reads:
    input:
        r1 = "samples/{sample}_1.fastq.gz",
        r2 = "samples/{sample}_2.fastq.gz",
    output:
        r1 = f"{W['clean']}/{{sample}}_clean_1.fastq.gz",
        r2 = f"{W['clean']}/{{sample}}_clean_2.fastq.gz",
        html = f"{W['logs']}/{{sample}}.html",
        json = f"{W['logs']}/{{sample}}.json",
    log:
        f"{W['logs']}/trim_{{sample}}.log",
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
        r1 = f"{W['clean']}/{{sample}}_clean_1.fastq.gz",
        r2 = f"{W['clean']}/{{sample}}_clean_2.fastq.gz",
        index = REFP["index_flag"],
    output:
        bam = f"{W['mapped']}/{{sample}}_Aligned.sortedByCoord.out.bam",
        transcriptome_bam = f"{W['mapped']}/{{sample}}_Aligned.toTranscriptome.out.bam",
        counts = f"{W['counts']}/{{sample}}_ReadsPerGene.out.tab",
        sj = f"{W['counts']}/{{sample}}_SJ.out.tab",
        final_log = f"{W['logs']}/{{sample}}_Log.final.out",
    log:
        f"{W['logs']}/map_{{sample}}.log",
    threads: config["threads"]["map"]
    resources:
        mem_mb = config["mem_mb"]["map"],
        runtime = config["runtime_min"]["map"],
    params:
        out_prefix = f"{W['mapped']}/{{sample}}_",
        genome_dir = REFP["index_dir"],
        logs_dir = W["logs"],
        sort_ram = int(config["star"]["limitBAMsortRAM_gb"]) * 10**9,
        genome_load = "LoadAndKeep" if SHARED_GENOME else "NoSharedMemory",
    shell:
        # --sjdbGTFfile e --sjdbOverhang NAO aparecem aqui de proposito: o indice
        # ja foi construido com as juncoes (--star-sjdboverhang no build), entao
        # passa-las de novo apenas reinsere as mesmas juncoes a cada amostra. E
        # a insercao em tempo de mapeamento e incompatıvel com genoma em memoria
        # compartilhada, que precisa ser read-only.
        #
        # --outSAMattrRGline usa campos separados por ESPACO: o STAR os une com
        # tabs. Escrever "\\t" aqui produzia o literal 't' (o bash converte \t
        # fora de aspas), gerando um @RG com um unico campo ID: e sem SM:/PL:.
        """
        ulimit -n 16384

        STAR --genomeDir {params.genome_dir} \
            --genomeLoad {params.genome_load} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --sjdbScore 2 \
            --outSAMattributes NH HI AS nM XS \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.out_prefix} \
            --runThreadN {threads} \
            --quantMode TranscriptomeSAM GeneCounts \
            --limitBAMsortRAM {params.sort_ram} \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:TRUSEQ PL:ILLUMINA \
            2> {log}

        # O STAR escreve tudo sob {params.out_prefix}, ou seja em mapped/.
        # As rules declaram counts/ e logs/, entao movemos.
        mv {params.out_prefix}ReadsPerGene.out.tab {output.counts}
        mv {params.out_prefix}SJ.out.tab {output.sj}
        mv {params.out_prefix}Log.final.out {output.final_log}
        mv {params.out_prefix}Log.out {params.logs_dir}/{wildcards.sample}_Log.out
        mv {params.out_prefix}Log.progress.out {params.logs_dir}/{wildcards.sample}_Log.progress.out

        rm -rf {params.out_prefix}_STARgenome {params.out_prefix}_STARtmp
        """


# Garante @RG completo no BAM coordenado.
#
# O STAR acima ja grava o read group corretamente, entao para BAMs novos esta
# rule e so uma rede de seguranca e faz hardlink (instantaneo, sem duplicar
# disco). Ela existe porque BAMs de versoes antigas do workflow saem com o @RG
# malformado — um unico campo ID: contendo tudo, sem SM: nem PL: — e o
# MarkDuplicates e o BQSR do interomics_snakemake exigem esses dois campos.
#
# Reescrever o header custa uma passada de I/O; realinhar custaria horas.
rule fix_read_group:
    input:
        bam = f"{W['mapped']}/{{sample}}_Aligned.sortedByCoord.out.bam",
    output:
        bam = f"{W['mapped']}/{{sample}}_rg.bam",
    log:
        f"{W['logs']}/fix_rg_{{sample}}.log",
    threads: config["threads"].get("fix_rg", 4)
    resources:
        mem_mb = config["mem_mb"].get("fix_rg", 8000),
        runtime = config["runtime_min"].get("fix_rg", 120),
    shell:
        """
        set -euo pipefail
        exec 2> {log}

        # So reescreve se o header estiver de fato quebrado. Um @RG valido tem
        # SM: como campo proprio, delimitado por tab.
        if samtools view -H {input.bam} | grep -q '^@RG.*\\bSM:'; then
            echo "@RG ja valido; hardlink em vez de reescrever"
            ln -f {input.bam} {output.bam}
        else
            echo "@RG malformado; reescrevendo com addreplacerg"
            samtools view -H {input.bam} | grep '^@RG' || true
            samtools addreplacerg \
                -r "ID:{wildcards.sample}" \
                -r "SM:{wildcards.sample}" \
                -r "LB:TRUSEQ" \
                -r "PL:ILLUMINA" \
                -w -o {output.bam} -@ {threads} {input.bam}
        fi
        samtools quickcheck {output.bam}
        """


# Quantify transcripts with RSEM
rule quantify_transcripts:
    input:
        bam = f"{W['mapped']}/{{sample}}_Aligned.toTranscriptome.out.bam",
        index = REFP["index_flag"],
    output:
        results = f"{W['counts']}/{{sample}}.genes.results",
        cnt = f"{W['logs']}/{{sample}}.cnt",
    log:
        f"{W['logs']}/quantify_{{sample}}.log",
    threads: config["threads"]["quantify"]
    resources:
        mem_mb = config["mem_mb"]["quantify"],
        runtime = config["runtime_min"]["quantify"],
    params:
        index_prefix = REFP["index_prefix"],
        counts_dir = W["counts"],
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
            {params.counts_dir}/{wildcards.sample} 2> {log}

        cp {params.counts_dir}/{wildcards.sample}.stat/{wildcards.sample}.cnt {output.cnt}
        """


# Run MultiQC to generate quality control report
rule run_multiqc:
    input:
        fastp = expand(f"{W['logs']}/{{sample}}.json", sample=SAMPLES),
        star = expand(f"{W['logs']}/{{sample}}_Log.final.out", sample=SAMPLES),
        rsem = expand(f"{W['logs']}/{{sample}}.cnt", sample=SAMPLES),
    output:
        report = f"{W['logs']}/multiqc_report.html",
    log:
        f"{W['logs']}/multiqc.log",
    params:
        logs_dir = W["logs"],
    resources:
        mem_mb = 4000,
        runtime = 60,
    shell:
        "multiqc {params.logs_dir} -o {params.logs_dir} -f 2> {log}"


# Run DEG analysis with R
rule run_deg_analysis:
    input:
        counts = expand(f"{W['counts']}/{{sample}}_ReadsPerGene.out.tab", sample=SAMPLES),
        rsem_counts = expand(f"{W['counts']}/{{sample}}.genes.results", sample=SAMPLES),
        design = config["deg_analysis"]["experiment_design_path"],
        gene_info = REFP["gene_info"],
    output:
        rdata = f"{W['rdata']}/DEG.RData",
    params:
        alpha = config["deg_analysis"]["alpha"],
        lfc_threshold = config["deg_analysis"]["lfc_threshold"],
    log:
        f"{W['logs']}/deg_analysis.log",
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


def _entregaveis(wildcards):
    """O que precisa existir no work_dir antes de exportar."""
    alvos = [
        *expand(f"{W['counts']}/{{sample}}_ReadsPerGene.out.tab", sample=SAMPLES),
        *expand(f"{W['counts']}/{{sample}}.genes.results", sample=SAMPLES),
        f"{W['logs']}/multiqc_report.html",
    ]
    if KEEP_BAMS in ("coord", "all"):
        # o _rg.bam, nao o cru: e ele que tem @RG utilizavel downstream
        alvos += expand(f"{W['mapped']}/{{sample}}_rg.bam", sample=SAMPLES)
    if KEEP_BAMS == "all":
        alvos += expand(
            f"{W['mapped']}/{{sample}}_Aligned.toTranscriptome.out.bam", sample=SAMPLES
        )
    if config["deg_analysis"]["enabled"]:
        alvos.append(f"{W['rdata']}/DEG.RData")
    return alvos


# Copia os entregaveis do scratch para o final_dir e so entao escreve a
# sentinela. A sentinela e o alvo de `rule all`: enquanto ela nao existir, o
# workflow nao esta pronto, por mais que o work_dir pareca completo.
rule export_results:
    input:
        _entregaveis,
    output:
        done = touch(f"{FINAL}/.export_done"),
    log:
        f"{W['logs']}/export.log",
    params:
        work = WORK,
        final = FINAL,
        keep_bams = KEEP_BAMS,
    shell:
        """
        bash {workflow.basedir}/scripts/export_results.sh \
            {params.work} {params.final} {params.keep_bams} > {log} 2>&1
        """


#------------------------------------------------------------------------------
# Genoma compartilhado e limpeza do scratch
#------------------------------------------------------------------------------

onstart:
    if SHARED_GENOME:
        shell(
            "STAR --genomeLoad LoadAndExit --genomeDir %s "
            "--outFileNamePrefix %s/_shm_ > /dev/null 2>&1 || true"
            % (REFP["index_dir"], W["logs"])
        )
        print(f"[star] genoma carregado em memoria compartilhada: {REFP['index_dir']}")


def _descarrega_genoma():
    if SHARED_GENOME:
        shell(
            "STAR --genomeLoad Remove --genomeDir %s "
            "--outFileNamePrefix %s/_shm_ > /dev/null 2>&1 || true"
            % (REFP["index_dir"], W["logs"])
        )
        print("[star] memoria compartilhada liberada")


onsuccess:
    _descarrega_genoma()
    if config["output"].get("cleanup_scratch", False):
        # O script confere os entregaveis no final_dir ANTES de apagar.
        shell(
            "bash {workflow.basedir}/scripts/cleanup_scratch.sh %s %s %s"
            % (WORK, FINAL, KEEP_BAMS)
        )


onerror:
    # Em falha o work_dir e PRESERVADO: e nele que esta o material para
    # diagnosticar e retomar.
    _descarrega_genoma()
    print(f"[erro] work_dir preservado para diagnostico: {WORK}")
