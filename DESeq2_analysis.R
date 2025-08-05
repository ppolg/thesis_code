# Running DESEQ2 like a champ to get bad results outta bad data
# (No, this isn't a waste of taxpayer money!)
#


# Load libraries
library(here)
library(ggplot2)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(pheatmap)

# Def



# Init

directory <- paste(here::here("seq/RNAseq"))

# Load data

#CDS
sampleFiles_CDS <- grep("CDS.count",list.files(directory),value=TRUE)
sampleTypes_CDS <- sub("[12].*","\\1",sampleFiles_CDS)
sampleTimes_CDS <- sub(".*-([^.]+)_.*","\\1",sampleFiles_CDS)
sampleTable_CDS <- data.frame(sampleName = sampleFiles_CDS,
                               fileName = sampleFiles_CDS,
                               type = sampleTypes_CDS,
                               time = sampleTimes_CDS)

sampleTable_CDS$type <- factor(sampleTable_CDS$type)

sampleTable_CDS30 <- sampleTable_CDS[sampleTable_CDS$time == 30,]
sampleTable_CDS60 <- sampleTable_CDS[sampleTable_CDS$time == 60,]

dds_CDS <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS,
                                       directory = directory,
                                       design = formula(~ time + type))

dds_CDS30 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS30,
                                         directory = directory,
                                         design = ~ type)

dds_CDS60 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS60,
                                         directory = directory,
                                         design = ~ type)

#tRNA
sampleFiles_tRNA <- grep("tRNA.count",list.files(directory),value=TRUE)
sampleTypes_tRNA <- sub("[12].*","\\1",sampleFiles_tRNA)
sampleTimes_tRNA <- sub(".*-([^.]+)_.*","\\1",sampleFiles_tRNA)
sampleTable_tRNA <- data.frame(sampleName = sampleFiles_tRNA,
                               fileName = sampleFiles_tRNA,
                               type = sampleTypes_tRNA,
                               time = sampleTimes_tRNA)

sampleTable_tRNA$type <- factor(sampleTable_tRNA$type)


sampleTable_tRNA30 <- sampleTable_tRNA[sampleTable_tRNA$time == 30,]
sampleTable_tRNA60 <- sampleTable_tRNA[sampleTable_tRNA$time == 60,]

dds_tRNA <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA,
                                       directory = directory,
                                       design = formula(~ time + type))

dds_tRNA30 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA30,
                                       directory = directory,
                                       design = ~ type)

dds_tRNA60 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA60,
                                         directory = directory,
                                         design = ~ type)

# Pre-filter

keep_CDS <- rowSums(counts(dds_CDS)) >= 10
keep_CDS30 <- rowSums(counts(dds_CDS30)) >= 10
keep_CDS60 <- rowSums(counts(dds_CDS60)) >= 10

dds_CDS <- dds_CDS[keep_CDS,]
dds_CDS30 <- dds_CDS30[keep_CDS30,]
dds_CDS60 <- dds_CDS60[keep_CDS60,]

keep_tRNA <- rowSums(counts(dds_tRNA)) >= 10
dds_tRNA <- dds_tRNA[keep_tRNA,]

# DE analysis

dds_CDS <- DESeq(dds_CDS)
res_CDS <- results(dds_CDS)

dds_CDS30 <- DESeq(dds_CDS30)
res_CDS30 <- results(dds_CDS30)
res_CDS30 <- res_CDS30[order(res_CDS30$pvalue),]
LFC_CDS30 <- lfcShrink(dds_CDS30, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS30 <- LFC_CDS30[order(LFC_CDS30$pvalue),]

dds_CDS60 <- DESeq(dds_CDS60)
res_CDS60 <- results(dds_CDS60)
res_CDS60 <- res_CDS60[order(res_CDS60$pvalue),]
LFC_CDS60 <- lfcShrink(dds_CDS60, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS60 <- LFC_CDS60[order(LFC_CDS60$pvalue),]



dds_tRNA <- DESeq(dds_tRNA)
res_tRNA <- results(dds_tRNA)

dds_tRNA30 <- DESeq(dds_tRNA30)
res_tRNA30 <- results(dds_tRNA30)
res_tRNA30 <- res_tRNA30[order(res_tRNA30$pvalue),]
LFC_tRNA30 <- lfcShrink(dds_tRNA30, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA30 <- LFC_tRNA30[order(LFC_tRNA30$log2FoldChange),]

dds_tRNA60 <- DESeq(dds_tRNA60)
res_tRNA60 <- results(dds_tRNA60)
res_tRNA60 <- res_tRNA60[order(res_tRNA60$pvalue),]
LFC_tRNA60 <- lfcShrink(dds_tRNA60, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA60 <- LFC_tRNA60[order(LFC_tRNA60$log2FoldChange),]

LFC_CDS60_table <- as.data.frame(LFC_CDS60)
LFC_CDS30_table <- as.data.frame(LFC_CDS30)

vsd_CDS <- vst(dds_CDS, blind=FALSE)
vsd_tRNA <- varianceStabilizingTransformation(dds_tRNA)

# Plot

plotMA(LFC_tRNA30, ylim=c(-2,2))
plotMA(LFC_tRNA60, ylim=c(-2,2))
plotMA(LFC_CDS30, ylim=c(-2,2))
plotMA(LFC_CDS60, ylim=c(-2,2))

meanSdPlot(assay(vsd_CDS))


sampleDists_CDS <- dist(t(assay(vsd_CDS)))
sampleDistMatrix_CDS <- as.matrix(sampleDists_CDS)
rownames(sampleDistMatrix_CDS) <- paste(vsd_CDS$type, vsd_CDS$time, sep="-")
colnames(sampleDistMatrix_CDS) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_CDS,
         clustering_distance_rows=sampleDists_CDS,
         clustering_distance_cols=sampleDists_CDS,
         col=colors)

plotPCA(vsd_CDS, intgroup=c("type", "time"))

# SELECTED DATASETS WITH ONLY 1 EMPTY
directory1 <- paste(here::here("RNAseq/Try1"))
directory2 <- paste(here::here("RNAseq/Try2"))

#CDS
sampleFiles_CDS1 <- grep("CDS.count",list.files(directory1),value=TRUE)
sampleTypes_CDS1 <- sub("[12].*","\\1",sampleFiles_CDS1)
sampleTimes_CDS1 <- sub(".*-([^.]+)_.*","\\1",sampleFiles_CDS1)
sampleTable_CDS1 <- data.frame(sampleName = sampleFiles_CDS1,
                              fileName = sampleFiles_CDS1,
                              type = sampleTypes_CDS1,
                              time = sampleTimes_CDS1)

sampleTable_CDS1$type <- factor(sampleTable_CDS1$type)

sampleTable_CDS301 <- sampleTable_CDS1[sampleTable_CDS1$time == 30,]
sampleTable_CDS601 <- sampleTable_CDS1[sampleTable_CDS1$time == 60,]

dds_CDS1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS1,
                                      directory = directory1,
                                      design = formula(~ time + type))

dds_CDS301 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS301,
                                        directory = directory1,
                                        design = ~ type)

dds_CDS601 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS601,
                                        directory = directory1,
                                        design = ~ type)

#tRNA
sampleFiles_tRNA1 <- grep("tRNA.count",list.files(directory1),value=TRUE)
sampleTypes_tRNA1 <- sub("[12].*","\\1",sampleFiles_tRNA1)
sampleTimes_tRNA1 <- sub(".*-([^.]+)_.*","\\1",sampleFiles_tRNA1)
sampleTable_tRNA1 <- data.frame(sampleName = sampleFiles_tRNA1,
                               fileName = sampleFiles_tRNA1,
                               type = sampleTypes_tRNA1,
                               time = sampleTimes_tRNA1)

sampleTable_tRNA1$type <- factor(sampleTable_tRNA1$type)


sampleTable_tRNA301 <- sampleTable_tRNA1[sampleTable_tRNA1$time == 30,]
sampleTable_tRNA601 <- sampleTable_tRNA1[sampleTable_tRNA1$time == 60,]

dds_tRNA1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA1,
                                       directory = directory1,
                                       design = formula(~ time + type))

dds_tRNA301 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA301,
                                         directory = directory1,
                                         design = ~ type)

dds_tRNA601 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA601,
                                         directory = directory1,
                                         design = ~ type)


dds_CDS1 <- DESeq(dds_CDS1)
res_CDS1 <- results(dds_CDS1)

dds_CDS301 <- DESeq(dds_CDS301)
res_CDS301 <- results(dds_CDS301)
res_CDS301 <- res_CDS301[order(res_CDS301$pvalue),]
LFC_CDS301 <- lfcShrink(dds_CDS301, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS301 <- LFC_CDS301[order(LFC_CDS301$pvalue),]

dds_CDS601 <- DESeq(dds_CDS601)
res_CDS601 <- results(dds_CDS601)
res_CDS601 <- res_CDS601[order(res_CDS601$pvalue),]
LFC_CDS601 <- lfcShrink(dds_CDS601, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS601 <- LFC_CDS601[order(LFC_CDS601$pvalue),]

dds_tRNA1 <- DESeq(dds_tRNA1)
res_tRNA1 <- results(dds_tRNA1)

dds_tRNA301 <- DESeq(dds_tRNA301)
res_tRNA301 <- results(dds_tRNA301)
res_tRNA301 <- res_tRNA301[order(res_tRNA301$pvalue),]
LFC_tRNA301 <- lfcShrink(dds_tRNA301, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA301 <- LFC_tRNA301[order(LFC_tRNA301$log2FoldChange),]

dds_tRNA601 <- DESeq(dds_tRNA601)
res_tRNA601 <- results(dds_tRNA601)
res_tRNA601 <- res_tRNA601[order(res_tRNA601$pvalue),]
LFC_tRNA601 <- lfcShrink(dds_tRNA601, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA601 <- LFC_tRNA601[order(LFC_tRNA601$log2FoldChange),]

LFC_CDS601_table <- as.data.frame(LFC_CDS601)
LFC_CDS301_table <- as.data.frame(LFC_CDS301)



#CDS
sampleFiles_CDS2 <- grep("CDS.count",list.files(directory2),value=TRUE)
sampleTypes_CDS2 <- sub("[12].*","\\1",sampleFiles_CDS2)
sampleTimes_CDS2 <- sub(".*-([^.]+)_.*","\\1",sampleFiles_CDS2)
sampleTable_CDS2 <- data.frame(sampleName = sampleFiles_CDS2,
                               fileName = sampleFiles_CDS2,
                               type = sampleTypes_CDS2,
                               time = sampleTimes_CDS2)

sampleTable_CDS2$type <- factor(sampleTable_CDS2$type)

sampleTable_CDS302 <- sampleTable_CDS2[sampleTable_CDS2$time == 30,]
sampleTable_CDS602 <- sampleTable_CDS2[sampleTable_CDS2$time == 60,]

dds_CDS2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS2,
                                       directory = directory2,
                                       design = formula(~ time + type))

dds_CDS302 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS302,
                                         directory = directory2,
                                         design = ~ type)

dds_CDS602 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_CDS602,
                                         directory = directory2,
                                         design = ~ type)

#tRNA
sampleFiles_tRNA2 <- grep("tRNA.count",list.files(directory2),value=TRUE)
sampleTypes_tRNA2 <- sub("[12].*","\\1",sampleFiles_tRNA2)
sampleTimes_tRNA2 <- sub(".*-([^.]+)_.*","\\1",sampleFiles_tRNA2)
sampleTable_tRNA2 <- data.frame(sampleName = sampleFiles_tRNA2,
                                fileName = sampleFiles_tRNA2,
                                type = sampleTypes_tRNA2,
                                time = sampleTimes_tRNA2)

sampleTable_tRNA2$type <- factor(sampleTable_tRNA2$type)


sampleTable_tRNA302 <- sampleTable_tRNA2[sampleTable_tRNA2$time == 30,]
sampleTable_tRNA602 <- sampleTable_tRNA2[sampleTable_tRNA2$time == 60,]

dds_tRNA2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA2,
                                        directory = directory2,
                                        design = formula(~ time + type))

dds_tRNA302 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA302,
                                          directory = directory2,
                                          design = ~ type)

dds_tRNA602 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_tRNA602,
                                          directory = directory2,
                                          design = ~ type)


dds_CDS2 <- DESeq(dds_CDS2)
res_CDS2 <- results(dds_CDS2)

dds_CDS302 <- DESeq(dds_CDS302)
res_CDS302 <- results(dds_CDS302)
res_CDS302 <- res_CDS302[order(res_CDS302$pvalue),]
LFC_CDS302 <- lfcShrink(dds_CDS302, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS302 <- LFC_CDS302[order(LFC_CDS302$pvalue),]

dds_CDS602 <- DESeq(dds_CDS602)
res_CDS602 <- results(dds_CDS602)
res_CDS602 <- res_CDS602[order(res_CDS602$pvalue),]
LFC_CDS602 <- lfcShrink(dds_CDS602, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_CDS602 <- LFC_CDS602[order(LFC_CDS602$pvalue),]

dds_tRNA2 <- DESeq(dds_tRNA2)
res_tRNA2 <- results(dds_tRNA2)

dds_tRNA302 <- DESeq(dds_tRNA302)
res_tRNA302 <- results(dds_tRNA302)
res_tRNA302 <- res_tRNA302[order(res_tRNA302$pvalue),]
LFC_tRNA302 <- lfcShrink(dds_tRNA302, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA302 <- LFC_tRNA302[order(LFC_tRNA302$log2FoldChange),]

dds_tRNA602 <- DESeq(dds_tRNA602)
res_tRNA602 <- results(dds_tRNA602)
res_tRNA602 <- res_tRNA602[order(res_tRNA602$pvalue),]
LFC_tRNA602 <- lfcShrink(dds_tRNA602, coef="type_RpfB_vs_Empty", type="apeglm")
LFC_tRNA602 <- LFC_tRNA602[order(LFC_tRNA602$log2FoldChange),]

LFC_CDS602_table <- as.data.frame(LFC_CDS602)
LFC_CDS302_table <- as.data.frame(LFC_CDS302)





