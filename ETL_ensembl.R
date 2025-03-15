ensembl = data.table::fread("./reference/GRCh38_ensembl_113.gtf", skip = 5)

colnames(ensembl) = c("seqname", "source", "feature",
                      "start_pos", "end_pos", "score",
                      "strand", "frame", "attribute")

separate_info = ensembl %>% 
  rownames_to_column("index") %>% 
  select(index, attribute) %>% 
  mutate(attribute = str_squish(attribute)) %>% 
  mutate(attribute = str_replace(attribute, ";$", "")) %>% 
  mutate(attribute = str_replace_all(attribute, "\"", "")) %>% 
  mutate(attribute = str_replace_all(attribute, "; ", ";")) %>% 
  separate_rows(attribute, sep = ";") %>%
  separate(attribute, 
           into = c("key", "value"), 
           sep = " ", 
           fill = "right") %>%
  pivot_wider(names_from = key,
              values_from = value, 
              values_fill = NA,
              values_fn = 
                list(value = function(x) paste(unique(x), collapse = ",")))

ensembl = ensembl %>% 
  rownames_to_column("index") %>% 
  left_join(., separate_info, join_by(index == index)) %>% 
  select(-c(index, attribute))

ensembl %>%
  write.table(., 
              file = "./reference/ensembl_gene_annotation.tsv", 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)