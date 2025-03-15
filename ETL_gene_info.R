library(tidyverse)

gene_info = data.table::fread("./reference/ncbi_gene_info.gz",
                              select = c(1:3,5)) %>%
  rename(tax_id = 1) %>%
  filter(tax_id == 9606) %>%
  select(gene_id = GeneID,
         gene_symbol = Symbol,
        synonyms = Synonyms)

gene_info %>%
  write.table(., 
              file = "./reference/GRCh38_p14_gene_info.tsv", 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)