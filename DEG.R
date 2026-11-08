#load libraries
library(tidyverse)
library(DESeq2)

#command args
# Lido ANTES da leitura das contagens: os caminhos vinham fixos como
# "./output/counts" e "./output/rdata", herdados do pipeline.sh, e o Snakemake
# roda a rule com work_dir ABSOLUTO -- o script nao achava as contagens e nao
# tinha onde gravar. Os dois agora chegam por argumento; o default preserva o
# comportamento antigo para quem ainda chamar o script na mao.
args = commandArgs(trailingOnly = TRUE)

counts_dir = if (length(args) >= 5) args[5] else "./output/counts"
rdata_path = if (length(args) >= 6) args[6] else "./output/rdata/DEG.RData"

#read counts
# sample_names e montado explicitamente em vez de varrer ls() depois. O ls()
# recolhia TUDO que existisse no ambiente: bastava uma variavel nova (como os
# dois caminhos acima) para ela entrar no bind_cols e quebrar a matriz de
# contagens. Com a lista explicita, acrescentar variavel aqui deixa de ser
# armadilha.
sample_names = character(0)

for (i in list.files(path = counts_dir,
                     pattern = "_ReadsPerGene.out.tab")){
  sample_name = str_replace(i, "_ReadsPerGene.out.tab", "")

  assign(paste0(sample_name),
        data.table::fread(file.path(counts_dir, i), skip = 4) %>%
          select(V1,
                 !!paste0(sample_name) := V2) %>% #Using := with !! for dynamic column naming
          arrange(V1) %>%
          column_to_rownames("V1")
        ) 
  sample_names = c(sample_names, sample_name)
}
rm(i, sample_name)

if (length(sample_names) == 0) {
  stop("Nenhum _ReadsPerGene.out.tab encontrado em: ", counts_dir)
}

# sort() reproduz a ordem alfabetica que o ls() dava antes.
raw_counts = do.call(bind_cols, mget(sort(sample_names)))

#read experiment design
exp_design = read_csv(args[1]) %>%
  column_to_rownames("sample") %>%
  mutate(condition = fct_relevel(condition, 
                                 c("control", "case")))

# Alinha as colunas das contagens a ordem do exp_design.
# As contagens sao lidas em ordem alfabetica (era o que o ls() dava), e o
# exp_design vem na ordem em que o usuario escreveu o CSV. O DESeq2 exige as
# duas na MESMA ordem e aborta quando divergem -- ou seja, o DEG so rodava se o
# CSV estivesse por acaso em ordem alfabetica. Reordenar aqui tira essa
# exigencia implicita do usuario.
faltando_no_design = setdiff(colnames(raw_counts), rownames(exp_design))
faltando_nas_contagens = setdiff(rownames(exp_design), colnames(raw_counts))
if (length(faltando_no_design) > 0 || length(faltando_nas_contagens) > 0) {
  stop("exp_design e contagens nao descrevem as mesmas amostras.\n",
       "  sem linha no design : ", paste(faltando_no_design, collapse = ", "), "\n",
       "  sem contagem        : ", paste(faltando_nas_contagens, collapse = ", "))
}
raw_counts = raw_counts[, rownames(exp_design), drop = FALSE]

#create DESeqDataSet object
dds = DESeqDataSetFromMatrix(countData = raw_counts,
                             colData = exp_design,
                             design = ~ condition)

#estimate size factors
dds = estimateSizeFactors(dds)

#remove transcript with low counts
smallest_group_size = exp_design %>% 
  group_by(condition) %>% 
  summarise(n = n()) %>% 
  summarise(min_n = min(n)) %>% 
  pull(min_n)
  
keep = rowSums(counts(dds) >= 10) >=  smallest_group_size
dds = dds[keep, ]

#export raw_counts
dds_counts = dds@assays@data$counts

#differential expression analysis
dds = DESeq(dds)

dds_results = results(dds, alpha = as.numeric(args[2]), 
                           lfcThreshold = as.numeric(args[3]), 
                           contrast = c("condition", "case", "control"))

dds_degs = data.frame(dds_results@listData,
                      row.names = dds_results@rownames)

#get entrez id
#args[4] = caminho do gene_info; default preserva o comportamento do pipeline.sh
gene_info_path = if (length(args) >= 4) args[4] else "./reference/GRCh38_p14_gene_info.tsv"
gene_info = read_delim(gene_info_path)

dds_degs  = dds_degs %>%
  rownames_to_column("gene_symbol") %>%
  left_join(., gene_info, join_by(gene_symbol == gene_symbol))

#reactome with significant DEGs
dds_signif = dds_degs %>%
  filter(padj < as.numeric(args[2])) %>%
  filter(abs(log2FoldChange) > as.numeric(args[3]))

reactome_total = ReactomePA::enrichPathway(gene = dds_signif$gene_id,
                                           pvalueCutoff = as.numeric(args[2]), 
                                           readable = TRUE)

dds_reactome = data.frame(reactome_total)

#save DESEq2 analsyis
dir.create(dirname(rdata_path), recursive = TRUE, showWarnings = FALSE)
save.image(rdata_path)
