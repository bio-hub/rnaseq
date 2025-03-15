#load libraries
library(tidyverse)
library(DESeq2)
# library(ReactomePA)
# library(biomaRt)

#read counts
for (i in list.files(path = "./output/counts", 
                     pattern = "_ReadsPerGene.out.tab")){
  sample_name = str_replace(i, "_ReadsPerGene.out.tab", "")

  assign(paste0(sample_name),
        data.table::fread(paste0("./output/counts/",i), skip = 4) %>%
          select(V1,
                 !!paste0(sample_name) := V2) %>% #Using := with !! for dynamic column naming
          arrange(V1) %>%
          column_to_rownames("V1")
        ) 
}
rm(i, sample_name)

raw_counts = eval(parse(
  text = paste0("bind_cols(", paste0(ls(), collapse = ","), ")")))

#command args
args = commandArgs(trailingOnly = TRUE)
args = c("./samples/GSE185919_exp_design.csv", "38")

#read experiment design
exp_design = read_csv(args[1]) %>%
  column_to_rownames("sample") %>%
  mutate(condition = fct_relevel(condition, 
                                 c("control", "case")))

#create DESeqDataSet object
dds = DESeqDataSetFromMatrix(countData = raw_counts,
                             colData = exp_design,
                             design = ~ condition)

#estimate size factors
dds = estimateSizeFactors(dds)

####Exploratory analysis and visualization
nrow(dds)
dds = dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)

###Differential Expression analysis
dds = DESeq(dds)

#View counts
dds_counts = dds@assays@data$counts

#DEG analysis
dds_results = results(dds, alpha = 0.05, 
                           lfcThreshold =1, 
                           contrast = c("condition", "case", "control"))

dds_degs = data.frame(dds_results@listData,
                      row.names = dds_results@rownames)

#get entrez id



# ensembl = biomaRt::useMart("ensembl", 
#                            dataset = "hsapiens_gene_ensembl")

# attributes = c("hgnc_symbol", 
#                "entrezgene_id")

# filters = c("hgnc_symbol")

# values = list(hgnc_symbol = rownames(dds_degs))

# query_entrezid = biomaRt::getBM(attributes = attributes, 
#                                 filters = filters, 
#                                 values = values, 
#                                 mart = ensembl)

#reactome
reactome_total = enrichPathway(gene=c(dep_downregulated, dep_upregulated),
                               pvalueCutoff = 0.05, 
                               readable=TRUE)

#save DESEq2 analsyis
save.image("./output/rdata/DEG.RData")
