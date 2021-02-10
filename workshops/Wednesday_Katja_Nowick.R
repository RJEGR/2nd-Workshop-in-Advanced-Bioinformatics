# DESeq2, edgR, biomaRt, TopGo, Rgraphviz, wTO, CoDiNA

# ==============
## Checking and Load packages ----
# ==============
.cran_packages <- c("wTO", "CoDiNA") # "tidyverse"
.bioc_packages <- c("edgeR","DESeq2", "biomaRt", "topGO", "Rgraphviz", "wTO")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

#
setwd("~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/workshops/katja/")

readcounts=read.table("Tutorial_macaques_readcounts.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)

dim(readcounts)

coldata=read.table("Tutorial_macaques_ids.txt", header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
dim(coldata)

coldata$condition = factor(coldata$condition)


dds=DESeqDataSetFromMatrix(countData = readcounts, colData = coldata, design = ~ condition)
dds

# filtering step

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds

dds$condition = factor(dds$condition, levels = c("NC","LPS"))


# DESEQ test

dds = DESeq(dds)

result = results(dds)
result

mcols(result)$description
resultLFC = lfcShrink(dds, coef="condition_LPS_vs_NC", type="normal")

summary(result)

resultOrdered = result[order(result$pvalue),]
sum(result$padj < 0.1, na.rm=TRUE)

sum(result$padj < 0.05, na.rm=TRUE)

result005 = results(dds, alpha=0.05)

summary(result005)

# plotMA(result, ylim=c(-2,2))

plotCounts(dds, gene=which.min(result$padj), intgroup="condition")

# write.csv(as.data.frame(resultOrdered), file="DESeq2_DEgenes_condition_LPS_NC.csv")

# test now edgeR

library(edgeR)

count_edgeR_obj=DGEList(counts=readcounts, group=coldata$condition)

count_edgeR_obj

count_edgeR_obj=estimateCommonDisp(count_edgeR_obj)

count_edgeR_obj=estimateTagwiseDisp(count_edgeR_obj)

edgeR_DEgenes=exactTest(count_edgeR_obj)

topTags(edgeR_DEgenes, sort.by = "logFC")

edgeR_DEgenesTable=edgeR_DEgenes$table

head(edgeR_DEgenesTable)

signedgeR_DEgenes=edgeR_DEgenesTable[edgeR_DEgenesTable[,3]<0.05,]

edgeROrdered <- edgeR_DEgenesTable[order(edgeR_DEgenesTable$PValue),]

# write.csv(as.data.frame(edgeR_DEgenesTable), file="edgeR_DEgenes_condition_LPS_NC.csv")

head(edgeROrdered)

head(resultOrdered)

sum(edgeROrdered$PValue < 0.05, na.rm=TRUE)
sum(resultOrdered$padj < 0.05, na.rm=TRUE)

#

# dim(signedgeR_DEgenes)
# dim(as.data.frame(resultOrdered)[resultOrdered$padj < 0.05, ])

# DESeq2 with interaction between factors ----

coldata$rank = as.factor(coldata$rank)
dds_interact=DESeqDataSetFromMatrix(countData = readcounts, 
                                    colData = coldata, 
                                    design = ~ condition + rank + condition:rank)
dds_interact

dds_interact = DESeq(dds_interact)

# Network

require(wTO)
require(magrittr)

# Splitting the input file into one for control (NC) and one for treatment (LPS)
# First, collecting all NC samples:

NC = readcounts[,coldata$condition == 'NC']
dim(NC)

# Removing genes that have less than 10 counts:
  
NC = NC[rowSums(NC)> 10,]
dim(NC)

# Then, collecting all LPS samples and removing genes that have less than 10 counts:
  
LPS = readcounts[,coldata$condition == 'LPS']
dim(LPS)

LPS = LPS[rowSums(LPS)> 10,]
dim(LPS)

# We can now select only the significant genes from the tables with LPS and NC samples. But before that, we will have a look at the results from above again:

summary(result)

DE_genes = subset(row.names(result), result$padj<0.01)
NC = subset(NC, row.names(NC) %in% DE_genes)
LPS = subset(LPS, row.names(LPS) %in% DE_genes)

# We want to make a TF wTO network. Thus, we need to retrieve the information about TF genes. The list of TFs is in the file “TFs.txt”. Let’s read this file into our R session.

TFs = read.table("TFs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

dim(TFs)

TFs = TFs[,1]
length(TFs)

# We obtained 1834 TFs. However, their IDs are GeneSymbols, while our table with readcounts used Ensembl gene IDs. So, we cannot easily match these IDs to identify the TFs in our LPS and NC table.

# We will thus take a short digression into Biomart, with which we can annotate our genes. Here, we want to change all Ensembl gene IDs in the LPS and NC tables to GeneSymbols. This way, we can match them with the TF table and our resulting network will display the TF GeneSymbols at the nodes.

# retrive gene annot via BiomaRt
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
mart = useDataset("mmulatta_gene_ensembl", useMart("ensembl"))
expressedGenes=row.names(result)

# listAttributes(mart)

GeneSymbols = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", 'external_gene_name'),values=expressedGenes,mart= mart)
dim(GeneSymbols)
head(GeneSymbols)

library(plyr)
LPS$ensembl_gene_id = row.names(LPS)
LPS = join(LPS, GeneSymbols, type = 'inner', match = 'first') 

LPS = LPS[!duplicated(LPS$external_gene_name), ]
row.names(LPS) = LPS$external_gene_name
head(LPS)

LPS = LPS[, -c(1,27)]
head(LPS)

NC$ensembl_gene_id = row.names(NC)
NC = join(NC, GeneSymbols, type = 'inner', match = 'first') 

NC = NC[!duplicated(NC$external_gene_name), ]
row.names(NC) = NC$external_gene_name
head(NC)

NC = NC[, -c(1,27)]
head(NC)

# Constructing the networks ----

# Now we have gathered all information to construct the networks. wTO performs a bootstrapping analysis to evaluate how likely it is that an inferred link is real. Note, that in a real analysis, you should run 1000 bootstraps, i.e. n should be 1000. Here we will run it only with n=10 to save time.

# First, constructing the NC network
# TFs will be the nodes, correlations will be calculated for all TFs with all expressed genes. The output, network_NC, will be an object with information on nodes, signed and absolute wTO values and p-values.

# This should take less than 5 minutes.

network_NC = wTO.Complete(n = 10, k = 5,  Data = NC, 
                          method_resampling = 'Bootstrap', 
                          Overlap = TFs, method = 's', plot = F) 

network_NC = network_NC$wTO
head(network_NC)

network_NC$wTO = ifelse(network_NC$Padj_sig<0.05, network_NC$wTO_sign, 0 )
head(network_NC)

# From the output table, remove everything but the info on nodes and the wTO between two nodes.

network_NC = network_NC[,c(1:2,9)] %>% as.data.frame()
head(network_NC)

network_LPS = wTO.Complete(n = 10,  Data = LPS, 
                           method_resampling = 'Bootstrap', 
                           Overlap = TFs, method = 's', plot = F)

network_LPS =network_LPS$wTO
network_LPS$wTO = ifelse(network_LPS$Padj_sig <0.05, network_LPS$wTO_sign, 0 )

network_LPS = network_LPS[,c(1,2,9)] %>% as.data.frame()

head(network_LPS)
require(CoDiNA)

network_NC$Node.1 = as.character(network_NC$Node.1)
network_NC$Node.2 = as.character(network_NC$Node.2)

Diff_LPS_NC = MakeDiffNet(Data = list(network_NC, network_LPS),
                          Code = c('NC', 'LPS'))


Diff_LPS_NC_clean = subset(Diff_LPS_NC, 
                           Diff_LPS_NC$Score_Phi_tilde/Diff_LPS_NC$Score_internal > 1)
Diff_LPS_NC_clean

DiffNodes = ClusterNodes(Diff_LPS_NC_clean, 
                         cutoff.external = 0, 
                         cutoff.internal = 1)

barplot(table(DiffNodes$Phi_tilde))
plot(Diff_LPS_NC_clean, layout = 'layout_with_drl')
