# Caminhos de referencia compartilhados entre Snakefile e reference.smk.
# Toda a arvore de referencia e derivada de config["reference"], nunca hardcoded.

import os

REF_ROOT = config["reference"]["dir"]
GENOME   = config["reference"]["genome_id"]
OVH      = int(config["star"]["sjdbOverhang"])
REF      = os.path.join(REF_ROOT, GENOME)

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
