# Construcao da referencia compartilhada. Rode UMA VEZ, via scripts/build_reference.sh
# (que segura um flock externo). Estas rules NAO existem no Snakefile principal:
# assim nenhuma execucao de analise pode reconstruir ou apagar a referencia.
#
#   bash scripts/build_reference.sh
#
# Escrita atomica em todos os passos: downloads em .part, indice em .staging_*,
# e a sentinela .done_sjdbN so e criada no fim. Nenhum consumidor ve estado parcial.

import os
from snakemake.utils import min_version

min_version("9.0")

configfile: "config.yaml"
include: "common.smk"

R = config["reference"]


rule reference_all:
    input:
        REFP["done"]


rule download_reference:
    output:
        genome = protected(REFP["fasta"]),
        annotation = protected(REFP["gtf"]),
        ncbi_gene_info = protected(REFP["ncbi_gene_info"]),
    log:
        REFP["logdir"] + "/download_reference.log",
    params:
        genome_url = R["genome_url"],
        annotation_url = R["annotation_url"],
        gene_info_url = R["gene_info_url"],
        ref = REF,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.ref}/fasta {params.ref}/gtf {params.ref}/gene_info {params.ref}/logs
        exec 2> {log}

        wget -c -O {output.genome}.gz.part "{params.genome_url}"
        mv {output.genome}.gz.part {output.genome}.gz

        wget -c -O {output.annotation}.gz.part "{params.annotation_url}"
        mv {output.annotation}.gz.part {output.annotation}.gz

        wget -c -O {output.ncbi_gene_info}.part "{params.gene_info_url}"
        mv {output.ncbi_gene_info}.part {output.ncbi_gene_info}

        gzip -dc {output.genome}.gz > {output.genome}.tmp
        mv {output.genome}.tmp {output.genome}
        rm -f {output.genome}.gz

        gzip -dc {output.annotation}.gz > {output.annotation}.tmp
        mv {output.annotation}.tmp {output.annotation}
        rm -f {output.annotation}.gz
        """


# Separado do download: o ETL le o .gz direto, entao o gene_info nunca e
# descomprimido (era o que quebrava a versao anterior desta etapa).
rule build_gene_info:
    input:
        gz = REFP["ncbi_gene_info"],
    output:
        tsv = protected(REFP["gene_info"]),
    log:
        REFP["logdir"] + "/build_gene_info.log",
    params:
        tax_id = R["tax_id"],
    shell:
        "Rscript {workflow.basedir}/ETL_gene_info.R {input.gz} {output.tsv} {params.tax_id} > {log} 2>&1"


rule build_rsem_index:
    input:
        genome = REFP["fasta"],
        annotation = REFP["gtf"],
    output:
        index = protected(REFP["index_flag"]),
    log:
        REFP["logdir"] + "/build_rsem_index.log",
    threads: config["threads"]["build_index"]
    resources:
        mem_mb = config["mem_mb"]["build_index"],
        runtime = config["runtime_min"]["build_index"],
    params:
        overhang = OVH,
        index_dir = REFP["index_dir"],
        genome_id = GENOME,
        ref = REF,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.ref}/logs
        # staging no mesmo filesystem (xfs) => o mv final e um rename atomico
        stg={params.ref}/.staging_$(hostname)_$$
        rm -rf "$stg"
        mkdir -p "$stg"
        trap 'rm -rf "$stg"' EXIT

        rsem-prepare-reference -p {threads} \
            --gtf {input.annotation} \
            --star \
            --star-sjdboverhang {params.overhang} \
            {input.genome} \
            "$stg/{params.genome_id}" > {log} 2>&1

        rm -rf {params.index_dir}
        mv "$stg" {params.index_dir}
        trap - EXIT
        """


# Escreve a proveniencia e so entao a sentinela. A partir daqui a referencia
# fica somente-leitura: um rebuild acidental falha com "Permission denied"
# em vez de apagar o indice. Para reconstruir de proposito:
#   chmod -R u+w /mnt/reference/rnaseq/<genome_id>
rule mark_done:
    input:
        index = REFP["index_flag"],
        gene_info = REFP["gene_info"],
        fasta = REFP["fasta"],
        gtf = REFP["gtf"],
    output:
        done = touch(REFP["done"]),
    params:
        manifest = REFP["manifest"],
        overhang = OVH,
        genome_url = R["genome_url"],
        annotation_url = R["annotation_url"],
        ref = REF,
        index_dir = REFP["index_dir"],
    shell:
        """
        set -euo pipefail
        m={params.manifest}
        rm -f "$m"
        date -Is                                  >> "$m"
        echo "host=$(hostname)"                   >> "$m"
        echo "STAR=$(STAR --version)"             >> "$m"
        echo "RSEM=$(rsem-calculate-expression --version | head -1)" >> "$m"
        echo "sjdbOverhang={params.overhang}"     >> "$m"
        echo "genome_url={params.genome_url}"     >> "$m"
        echo "annotation_url={params.annotation_url}" >> "$m"
        md5sum {input.fasta} {input.gtf} {input.gene_info} >> "$m"

        # congela os dados; logs/ e MANIFEST ficam graváveis para o proximo build
        chmod -R a-w {params.ref}/fasta {params.ref}/gtf {params.ref}/gene_info {params.index_dir} || true
        """
