library("tximport")
library("tidyverse")
library("DESeq2")
library("reshape2")
library("pheatmap")
library("RColorBrewer")
library("IHW")

#set working directory
setwd("C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpuratus/)

#path to the samples
dir_sp27 <- list.files(file.path("C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpuratus/sp27wtbra/"))
quant_files_sp27 <-paste0("C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpuratus/sp27wtbra/",dir_sp27,"/quant.sf")
names(quant_files_sp27) <- factor(c("sp27wt1", "sp27wt2", "sp27wt3", "sp27bra1", "sp27bra2", "sp27bra3"))

#conversion table for gene level analysis
tx2gene_sp <- read.table("tx2gene.txt")
txi_sp27 <- tximport(quant_files_sp27, tx2gene = tx2gene_sp, type="salmon", countsFromAbundance = "no")

#create sample table
sampleTable_sp27 <- data.frame(condition = factor(rep(c("control", "experiment"), each = 3)))
rownames(sampleTable_sp27) <- colnames(txi_sp27$counts)

#add experimental conditions
sampleTable_sp27$batch <- factor(c("1", "2", "3", "1", "2", "3"))
dds_sp27 <- DESeqDataSetFromTximport(txi_sp27, sampleTable_sp27, ~batch + condition)
dds_sp27$condition <- relevel(dds_sp27$condition, ref = "control")

#perform differential analysis
deseq_sp27 <- DESeq(dds_sp27)

#export and save results
res_sp27 <- results(deseq_sp27)
res27wtbra <- data.frame(res_sp27)
write.table(res27wtbra_tx2, "C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpurats/gene_level_results/res27wt_vs_bra.txt", sep="\t")

#Independent hypothesis weighting
resIHW_sp27 <- results(deseq_sp27_tx2, filterFun=ihw)

#Inspect results
summary(resIHW_sp27)

#Write results and subset results by padj value
write.table(resIHW_sp27, file="C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpurats/gene_level_results/sp27wt_bra_results_IHW.txt", sep="\t")
resIHW_sp27padj005 <- subset(resIHW_sp27, padj <= 0.05)
resIHW_sp27padj005 <- data.frame(resIHW_sp27padj005_tx2)
write.table(resIHW_sp27padj005, "C:/Users/sznuser15/Documents/PhD/Data/RNAseq_szn/S_purpurats/gene_level_results/IHW_sp27wt_vs_bra_padj005.txt", sep="\t")

#visualisation PCA
vst_sp27 <- vst(deseq_sp27)
qplot(PC1, PC2, color=condition, shape=batch, data=plotPCA(vst_sp27, intgroup = c("condition", "batch"), returnData=TRUE))

#clustering and dendrogram
sampleDists27 <- dist(t(assay(vst_sp27)))
sampleDistMatrix27 <- as.matrix(sampleDists27)
rownames(sampleDistMatrix27) <- paste(vst_sp27$condition, vst_sp27$batch, sep="-")
colnames(sampleDistMatrix27) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix27,
         clustering_distance_rows=sampleDists27,
         clustering_distance_cols=sampleDists27,
         col=colors)
