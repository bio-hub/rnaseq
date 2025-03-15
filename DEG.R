#command args
args = commandArgs(trailingOnly = TRUE)

#libraries
base_packages = c("data.table", "ggpubr", "devtools", "BiocManager")
#bioconductor_packages = c("GenomicFeatures", "GenomicAlignments", "DESeq2", "ReactomePA")

bioconductor_packages = c("DESeq2", "ReactomePA")

if (!require(base_packages, quietly = TRUE))
    install.packages(base_packages, 
    ask= FALSE, 
    repo = "https://cran-r.c3sl.ufpr.br/")

if (!require(bioconductor_packages, quietly = TRUE))
    BiocManager::install(bioconductor_packages, 
    ask = FALSE)

#load libraries
library(tidyverse)
#library(Rsamtools)
#library(GenomicFeatures)
#library(GenomicAlignments)
#library(BiocParallel)
library(DESeq2)
library(ReactomePA)

#read counts
for (i in list.files(path = "./output/counts", pattern = "_ReadsPerGene.out.tab")){
  assign(paste0(str_replace(i, "_ReadsPerGene.out.tab", "")),
        data.table::fread(paste0("./output/counts/",i),   skip = 4)
        )
}
rm(i)

#load input file
input = read.table("./tmp/input_2R.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
input$V1 = gsub(" ", "", input$V1)
input$V2 = gsub(" ", "", input$V2)

#DESeq2 import functions
bamfiles <- Rsamtools::BamFileList(input$V3, yieldSize=1000000)

#Defining gene models
gtffile <- "/home/operator/InterOmics/db/GRCh38_genes.gtf"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
ebg <- GenomicFeatures::exonsBy(txdb, by="gene")

#Read counting step
BiocParallel::register(BiocParallel::MulticoreParam(progressbar = TRUE, workers=56))

### Quantification
se <- GenomicAlignments::summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE)
save.image("./output/rdata/DEG.RData")

###DEG Starting from SummarizedExperiment

#se$condition1 = exp_design$condition1
#se$condition2 = exp_design$condition2
#se$condition3 = exp_design$condition3

se$group=input$V2
dds <- DESeqDataSet(se, design = ~group) # Set group to compare
dds <- estimateSizeFactors(dds)

####Exploratory analysis and visualization
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10, ]
nrow(dds)

###Differential Expression analysis
dds <- DESeq(dds)

#View counts
#View(dds@assays@data$counts)

### OUTPUT####
groups = unique(input$V2)

# 1 - Lista dos DEGs
for (i in 1:length(groups)){
  for (j in 1:length(groups)){
    if(i != j){
      nam  = paste0("DEG.",groups[i],".", groups[j])
      nam2 = paste0("ls.",nam)
      assign(nam,
             eval(parse(text=paste0("results(dds, alpha = 0.05, lfcThreshold =1, contrast = c('group','", groups[i],"','",groups[j],"'))")))
      )
      assign(nam2,
             data.frame(get(nam)@listData, SYMBOL = get(nam)@rownames, stringsAsFactors = FALSE)
      )
    }
  }
}

# 2 -Gene Expression para cada individuo - Normalizada
tempa = data.frame(dds@assays@data$counts)
colnames(tempa) = sub("_Aligned.sortedByCoord.out.bam","",names(tempa))
tempa = data.frame(SYMBOL = rownames(get(nam)),tempa) #Nao remover o obejto nam.
# Gene SYMBOL to ENTREZID
library(org.Hs.eg.db)
tempb = select(x = org.Hs.eg.db,keys = tempa$SYMBOL,columns = c("SYMBOL","ENTREZID"),keytype = "SYMBOL")
source(file = "../../../../transcriptome/output/functions/fix.geneID.R")
fix.geneID(df_input = tempb,df_name_output = "tempc")

genes.count = merge(tempc,tempa, by ="SYMBOL", all = TRUE)
# DEG Final object
for (i in ls(pattern = "ls.DEG.")){
  assign(i,
         merge(tempc,get(i),by ="SYMBOL", all = TRUE))
}
rm(tempa,tempb)

for (i in ls(pattern = "ls.DEG.")){
  nam = sub("ls.DEG.","",i)
  assign(paste0("up.",nam),
         get(i)[which(get(i)$log2FoldChange  > 1 & get(i)$padj < 0.05),]
  )
  assign(paste0("down.",nam),
         get(i)[which(get(i)$log2FoldChange < -1 & get(i)$padj < 0.05),]
  )
  
}

DEG = list()
for (i in ls(pattern = "ls.DEG.")){
  nam = sub("ls.DEG.","",i)
  # eval(parse(text = paste0("DEG = list(",nam,"= list(all = ls.DEG.",nam,", up = up.",nam,", down = down.",nam,"))")))
  eval(parse(text = paste0("DEG$",nam," = list(all = ls.DEG.",nam,", up = up.",nam,", down = down.",nam,",DESeqResults = DEG.",nam,")")))
  
}
rm(list = ls(pattern = "up\\.|down\\.|DEG\\."))



#### Reactome ####
source(file = "../../../../transcriptome/output/functions/deg.enrichPathway.R")

for (i in names(DEG)){
  deg.enrichPathway(comparisons = i)
}


#save rdata
save.image("./output/rdata/DEG.RData")

