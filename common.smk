# Caminhos de referencia compartilhados entre Snakefile e reference.smk.
# Toda a arvore de referencia e derivada de config["reference"], nunca hardcoded.

import os
import sys

from snakemake.exceptions import WorkflowError

REF_ROOT = config["reference"]["dir"]
GENOME   = config["reference"]["genome_id"]
REF      = os.path.join(REF_ROOT, GENOME)


def _resolve_sjdb_overhang(star):
    """sjdbOverhang a partir de star.read_length, com override explicito.

    O STAR quer comprimento_da_read - 1. Pedir o comprimento da read no config
    e derivar daqui evita que alguem escreva 150 onde o STAR espera 149.
    star.sjdbOverhang continua aceito e tem precedencia, para reads de
    comprimento misto (use o maior - 1) ou para reaproveitar um indice que ja
    foi construido com outro valor.
    """
    override = star.get("sjdbOverhang")
    if override is not None:
        overhang = int(override)
        origem = f"sjdbOverhang={overhang} (override explicito do config)"
        if star.get("read_length") is not None:
            print(
                f"[config] star.sjdbOverhang={overhang} tem precedencia; "
                f"star.read_length={star['read_length']} ignorado.",
                file=sys.stderr,
            )
    elif star.get("read_length") is not None:
        read_length = int(star["read_length"])
        if read_length < 2:
            raise WorkflowError(
                f"star.read_length={read_length} invalido: precisa ser >= 2."
            )
        overhang = read_length - 1
        origem = f"read_length={read_length} bp -> sjdbOverhang={overhang}"
    else:
        raise WorkflowError(
            "config star: informe star.read_length (comprimento das reads em bp).\n"
            "Opcionalmente use star.sjdbOverhang para fixar o valor na mao."
        )

    if overhang < 1:
        raise WorkflowError(f"sjdbOverhang={overhang} invalido: precisa ser >= 1.")
    return overhang, origem


OVH, OVH_ORIGEM = _resolve_sjdb_overhang(config["star"])


def _resolve_output_dirs(cfg):
    """work_dir (intermediarios, no scratch) e final_dir (entregaveis).

    Os intermediarios de RNA-seq dominam o volume: FASTQ limpos somam o mesmo
    que a entrada, e os BAMs mais ainda. O entregavel final (contagens,
    relatorios) e 1-2% disso. Por isso os dois vivem em volumes diferentes.
    """
    out = cfg.get("output", {}) or {}
    user = os.environ.get("USER") or os.environ.get("LOGNAME") or "nobody"
    project = out.get("project", "rnaseq")

    work = out.get("work_dir", "output").format(user=user, project=project)
    final = out.get("final_dir", "results").format(user=user, project=project)

    keep = str(out.get("keep_bams", "coord")).lower()
    if keep not in ("none", "coord", "all"):
        raise WorkflowError(
            f"output.keep_bams={keep!r} invalido: use none, coord ou all."
        )

    # Disco local do no (Tier 0). Vazio = tudo no work_dir, que e o
    # comportamento padrao e funciona em qualquer maquina, com ou sem SLURM.
    local = (out.get("local_scratch") or "").strip()
    if local:
        local = os.path.abspath(local.format(user=user, project=project))

    return os.path.abspath(work), os.path.abspath(final), keep, local


WORK, FINAL, KEEP_BAMS, LOCAL = _resolve_output_dirs(config)

if LOCAL:
    # EXIGE AGRUPAMENTO POR AMOSTRA. O FASTQ aparado fica no disco de UM no;
    # se trim_reads e map_reads cairem em nos diferentes, o segundo nao acha a
    # entrada do primeiro. O agrupamento vem do profile (--groups), nao do
    # Snakefile -- ver profiles/slurm/config.yaml.
    sys.stderr.write(
        f"[rnaseq] local_scratch ativo: {LOCAL}\n"
        "[rnaseq] Afeta APENAS temporarios de ferramenta (STAR _STARtmp, RSEM .temp),\n"
        "[rnaseq] que nascem e morrem dentro de um job. Saidas declaradas de regra\n"
        "[rnaseq] continuam no work_dir compartilhado, porque o Snakemake as confere\n"
        "[rnaseq] a partir do no orquestrador. Ver comentario em common.smk.\n"
    )

# Subdiretorios do work_dir. Tudo que e intermediario mora aqui.
#
# O RECORTE: o que pode ser local e o que NAO pode.
#
# Local exige que so a PROPRIA amostra consuma o arquivo. Compartilhado e
# obrigatorio para tudo que alguma regra le de TODAS as amostras:
#
#   run_multiqc       le logs/*.json de todas
#   run_deg_analysis  le counts/* de todas
#   export_results    le counts, logs e mapped/*_rg.bam de todas
#
#   clean   372 GB no lote 1  consumido so por map_reads da mesma amostra  -> LOCAL
#   mapped  520 GB no lote 1  export_results le *_rg.bam de todas          -> compartilhado
#   counts  logs  rdata       lidos por todas as regras de agregacao       -> compartilhado
#
# Mover `mapped` tambem exigiria que fix_read_group gravasse o *_rg.bam
# direto no volume compartilhado, mantendo so a saida crua do STAR local.
# Isso tiraria ~1,5 TB de trafego do array em vez dos ~744 GB de agora, e
# e a proxima refinada -- mas mexe na estrutura do Snakefile, nao so aqui.
#
# O QUE PODE SER LOCAL, E O QUE NAO PODE
#
# Descoberto em 2026-10-31 rodando 2 amostras de verdade: o Snakemake valida
# as SAIDAS DECLARADAS a partir do no ORQUESTRADOR. Um arquivo no /scratch de
# um compute node nao existe para o head, e o resultado foi:
#
#   MissingOutputException: Job 5 completed successfully,
#   but some output files are missing.
#   Removing output files of failed job map_reads since they might be corrupted
#
# O group job tinha terminado bem -- SLURM reportou COMPLETED 0:0 nos quatro
# passos -- e o Snakemake apagou 48 GB de trabalho valido por nao enxerga-lo.
# Agrupar NAO resolve: o job roda num no so, mas quem confere e o head.
#
# Portanto:
#   saida declarada de regra  -> SEMPRE compartilhada
#   temporario de ferramenta  -> pode ser local (o Snakemake nem sabe que existe)
#
# E por isso que `clean` voltou para o WORK e so `tmp` segue o LOCAL. Nao e
# perda: medido, o temporario do RSEM e ~30 GB por amostra contra ~4 GB do
# FASTQ aparado. O ganho maior estava no que nao se declara.
PESADO = LOCAL or WORK

W = {
    # Temporarios de ferramenta (STAR _STARtmp, RSEM .temp). Nao sao saida de
    # regra nenhuma: nascem e morrem dentro de um job, logo acompanham o
    # PESADO. Medido no teste de 2 amostras: ~30 GB por amostra so do RSEM,
    # mais o que o STAR usa para ordenar -- tudo isso ia para o array porque
    # as ferramentas criam o temporario AO LADO DA SAIDA quando nao se diz
    # outra coisa.
    "tmp":    f"{PESADO}/tmp",
    "clean":  f"{WORK}/clean",
    "mapped": f"{WORK}/mapped",
    "counts": f"{WORK}/counts",
    "logs":   f"{WORK}/logs",
    "rdata":  f"{WORK}/rdata",
}

# O sjdbOverhang entra no nome do indice e da sentinela para que indices com
# overhang diferente coexistam em vez de se sobrescrever.
REFP = {
    "fasta":          f"{REF}/fasta/{GENOME}.fa",
    "gtf":            f"{REF}/gtf/{GENOME}.gtf",
    "ncbi_gene_info": f"{REF}/gene_info/ncbi_gene_info.gz",
    "gene_info":      f"{REF}/gene_info/{GENOME}_gene_info.tsv",
    "index_dir":      f"{REF}/rsem_star_sjdb{OVH}",
    "index_prefix":   f"{REF}/rsem_star_sjdb{OVH}/{GENOME}",
    "index_flag":     f"{REF}/rsem_star_sjdb{OVH}/genomeParameters.txt",
    "manifest":       f"{REF}/MANIFEST.txt",
    "done":           f"{REF}/.done_sjdb{OVH}",
    "logdir":         f"{REF}/logs",
}
