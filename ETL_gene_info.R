# Extrai gene_id/symbol/sinonimos do gene_info do NCBI para uma especie.
#
# Uso: Rscript ETL_gene_info.R <gene_info.gz> <saida.tsv> [tax_id]
#
# Le o .gz diretamente (nao descomprima antes). Os defaults preservam o
# comportamento antigo para pipeline.sh, que chama o script sem argumentos.

library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

input_path  = if (length(args) >= 1) args[1] else "./reference/ncbi_gene_info.gz"
output_path = if (length(args) >= 2) args[2] else "./reference/GRCh38_p14_gene_info.tsv"
tax_id_keep = if (length(args) >= 3) as.integer(args[3]) else 9606L

gene_info = data.table::fread(input_path, select = c(1:3, 5)) %>%
  rename(tax_id = 1) %>%
  filter(tax_id == tax_id_keep) %>%
  select(gene_id = GeneID,
         gene_symbol = Symbol,
         synonyms = Synonyms)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

gene_info %>%
  write.table(.,
              file = output_path,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
