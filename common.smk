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
