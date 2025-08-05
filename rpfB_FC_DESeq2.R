
# Deseq2 for featurecounts table
# For RNAseq of rpfB overexpress in Msm

########
# Func #
########

theme_mycopore <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_blank(),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_text(angle=90,vjust = 3, size = rel(1)),
           axis.title.x = element_text(vjust = -0.1, size = rel(1)),
           axis.text = element_text(face="bold",size = rel(1.3)), 
           axis.line = element_line(colour="black",size=1),
           axis.ticks = element_line(colour="black",size = 1),
           axis.ticks.length = unit(.35, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "top",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.x = unit(0.3, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(12,6,6,6),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}


###########################################################################

########
# Init #
########

invisible(library(here))
invisible(library(DESeq2))
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))


gff <- get_gff_input("msmeg") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";product=|;identified=|;this=", 2)[,1], "Name=", 2)[,2],) %>%
  dplyr::mutate(product = str_split_fixed(str_split_fixed(attributes, ";uniProt_AC", 2)[,1], "product=", 2)[,2],)

gff2 <- get_gff_input("msmeg3") %>%
  filter(type == "gene") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "Name=", 2)[,2],) %>%
  dplyr::mutate(locus_name2 = sub(".*locus_tag=", "\\1",attributes)) %>%
  mutate(weird_locus = str_extract(attributes,"MSMEG_RS\\d{5}"))

gff2$realID <- apply(gff2[10:11],1, function(x) x[which.min(nchar(x))])


########
# Main #
########

##### Get files from featurecounts csv list #####

directory <- paste(here::here("Brindhaseqcounts"))

sampleFiles <- grep(".counts",list.files(directory),value=TRUE)
names <- sub(".counts","\\1",sampleFiles)
sampleTypes <- sub("_[02]h_.*","\\1",sampleFiles)
sampleTimes <- regmatches(sampleFiles,regexpr("[02]h",sampleFiles))
sampleTable <- tibble(sampleName = names,
                      fileName = sampleFiles,
                      time = sampleTimes,
                      type = sampleTypes,
                      reps = rep(c(1,2,3), times=length(sampleFiles)/3))

sampleTable$type <- factor(sampleTable$type)

# 2h here selects out outliers aswell!
sampleTable_2h <- sampleTable %>% filter(time == "2h", reps != 1)


sampleTable_trunc <- sampleTable %>% filter(sampleName %nin% c("304-1_S1"))
sampleTable_compare <- sampleTable_trunc %>% filter(type == "rpfB")
sampleTable_0h <- sampleTable_trunc %>% filter(time == "0h")


##### Full #####
# Loop to make matrix
samples <- data.table(GeneID = gff2$realID)
columns <- data.table(sample = sampleTable$sampleName,
                  type = sampleTable$type,
                  time = sampleTable$time,
                  reps = sampleTable$reps)

for(i in 1:nrow(sampleTable)){
  
  # get the file
  n <- read.csv(paste0(directory,"/",sampleTable$fileName[i]))[8] %>%
    dplyr::rename(!!sampleTable$sampleName[i] := 1)
  n[,1] = ifelse(n[,1]<10,0,n[,1])
  samples <- cbind(samples,n)
}

cts <- as.matrix(samples, rownames="GeneID")
coldata <- as.matrix(columns, rownames="sample")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ type+time)

dds <- DESeq(dds)
k <- resultsNames(dds)
res_dds <- results(dds)
#LFC_type <- lfcShrink(dds, coef="type_304.SW_vs_304", type="apeglm")
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("type","time"))
PCAdata <- plotPCA(vsd, intgroup=c("type", "time"), returnData=T)

fig_PCA <- ggplot(PCAdata, aes(size = 1.2)) +
  geom_point(aes(x = PC1, y= PC2, colour = group)) +
  theme_mycopore() +
  ggtitle("PCA of full dataset")+
  xlab("PC1: 68% variance") +
  ylab("PC2: 14% variance") +
  scale_colour_manual(values = c(colours$red,colours$orange,colours$green,colours$azure)) +
  theme(legend.title = element_blank()) +
  guides(size = "none")

print(fig_PCA)

##### For comparing 2h #####
samples_2h <- data.table(GeneID = gff2$realID)
columns_2h <- data.table(sample = sampleTable_2h$sampleName,
                      type = sampleTable_2h$type,
                      time = sampleTable_2h$time,
                      reps = sampleTable_2h$reps)

for(i in 1:nrow(sampleTable_2h)){

  # get the file
  n <- read.csv(paste0(directory,"/",sampleTable_2h$fileName[i]))[8] %>%
    dplyr::rename(!!sampleTable_2h$sampleName[i] := 1)
  n[,1] = ifelse(n[,1]<10,0,n[,1])
  samples_2h <- cbind(samples_2h,n)
}

cts_2h <- as.matrix(samples_2h, rownames="GeneID")
coldata_2h <- as.matrix(columns_2h, rownames="sample")

dds_2h <- DESeqDataSetFromMatrix(countData = cts_2h,
                              colData = coldata_2h,
                              design = ~ type)

dds_2h <- DESeq(dds_2h)
res_dds_2h <- results(dds_2h)
LFC_2h <- lfcShrink(dds_2h, coef="type_rpfB_vs_Empty", type="apeglm")
LFC_2h <- LFC_2h[order(LFC_2h$pvalue),]
vsd_2h <- vst(dds_2h, blind=FALSE)
LFC_2h_table <- as.data.frame(LFC_2h) %>%
  dplyr::arrange(padj) %>%
  filter(padj < 0.01) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  rownames_to_column()

h <- gff$product[match(x = LFC_2h_table$rowname, table=gff$locus_name)]
LFC_2h_table$prod <- h


fwrite(LFC_2h_table, file=paste0(directory,"/lfc_table.csv"), row.names = T)

##### Removed outliers #####
# Loop to make matrix
# samples_trunc <- data.table(GeneID = gff2$realID)
# columns_trunc <- data.table(sample = sampleTable_trunc$sampleName,
#                       type = sampleTable_trunc$type,
#                       reps = sampleTable_trunc$reps)
# 
# for(i in 1:nrow(sampleTable_trunc)){
#   
#   # get the file
#   n <- read.csv(paste0(directory,"/",sampleTable_trunc$fileName[i]))[8] %>%
#     dplyr::rename(!!sampleTable_trunc$sampleName[i] := 1)
#   n[,1] = ifelse(n[,1]<10,0,n[,1])
#   samples_trunc <- cbind(samples_trunc,n)
# }
# 
# cts_trunc <- as.matrix(samples_trunc, rownames="GeneID")
# coldata_trunc <- as.matrix(columns_trunc, rownames="sample")
# 
# dds_trunc <- DESeqDataSetFromMatrix(countData = cts_trunc,
#                               colData = coldata_trunc,
#                               design = ~ type)
# 
# dds_trunc <- DESeq(dds_trunc)
# res_dds_trunc <- results(dds_trunc)
# LFC_type_trunc <- lfcShrink(dds_trunc, coef="type_304.SW_vs_304", type="apeglm")
# 
# trunc_result_table <- as.data.frame(LFC_type_trunc) %>%
#   rownames_to_column("gene")
# 
# fwrite(trunc_result_table, file=paste0(directory,"/lfc_table.csv"), row.names = T)
# 
# 
# vsd_trunc <- vst(dds_trunc, blind=FALSE)
# plotPCA(vsd_trunc, intgroup=c("type"))
# PCAdata_trunc <- plotPCA(vsd_trunc, intgroup=c("type"), returnData=T)

##### Trunc, 0h vs 2h for changes #####
# Loop to make matrix
# samples_compare <- data.table(GeneID = gff2$realID)
# columns_compare <- data.table(sample = sampleTable_compare$sampleName,
#                             type = sampleTable_compare$type,
#                             time = sampleTable_compare$time,
#                             reps = sampleTable_compare$reps)
# 
# for(i in 1:nrow(sampleTable_compare)){
#   
#   # get the file
#   n <- read.csv(paste0(directory,"/",sampleTable_compare$fileName[i]))[8] %>%
#     dplyr::rename(!!sampleTable_compare$sampleName[i] := 1)
#   n[,1] = ifelse(n[,1]<10,0,n[,1])
#   samples_compare <- cbind(samples_compare,n)
# }
# 
# cts_compare <- as.matrix(samples_compare, rownames="GeneID")
# coldata_compare <- as.matrix(columns_compare, rownames="sample")
# 
# dds_compare <- DESeqDataSetFromMatrix(countData = cts_compare,
#                                     colData = coldata_compare,
#                                     design = ~ time)
# 
# dds_compare <- DESeq(dds_compare)
# res_dds_compare <- results(dds_compare)
# LFC_compare <- lfcShrink(dds_compare, coef="time_2h_vs_0h", type="apeglm")
# LFC_compare <- LFC_compare[order(LFC_compare$pvalue),]
# LFC_compare_table <- as.data.frame(LFC_compare)

##### Leaky plasmid test of 0h #####

# samples_0h <- data.table(GeneID = gff2$realID)
# columns_0h <- data.table(sample = sampleTable_0h$sampleName,
#                          type = sampleTable_0h$type,
#                          time = sampleTable_0h$time,
#                          reps = sampleTable_0h$reps)
# 
# for(i in 1:nrow(sampleTable_0h)){
#   
#   # get the file
#   n <- read.csv(paste0(directory,"/",sampleTable_0h$fileName[i]))[8] %>%
#     dplyr::rename(!!sampleTable_0h$sampleName[i] := 1)
#   n[,1] = ifelse(n[,1]<10,0,n[,1])
#   samples_0h <- cbind(samples_0h,n)
# }
# 
# cts_0h <- as.matrix(samples_0h, rownames="GeneID")
# coldata_0h <- as.matrix(columns_0h, rownames="sample")
# 
# dds_0h <- DESeqDataSetFromMatrix(countData = cts_0h,
#                                  colData = coldata_0h,
#                                  design = ~ type)
# 
# dds_0h <- DESeq(dds_0h)
# res_dds_0h <- results(dds_0h)
# LFC_0h <- lfcShrink(dds_0h, coef="type_rpfB_vs_Empty", type="apeglm")
# LFC_0h <- LFC_0h[order(LFC_0h$pvalue),]
# vsd_0h <- vst(dds_0h, blind=FALSE)
# LFC_0h_table <- as.data.frame(LFC_0h)
# 
# fwrite(LFC_0h_table, file=paste0(directory,"/0h_comp.csv"), row.names = T)



##### MA/volcano plots #####

cap_limit <- 2 #+/- limit on y axis

###### rpfB 0 vs 2 ######
ma_rpfB <- as.data.frame(results(dds_2h, contrast = c("type","Empty","rpfB"), cooksCutoff = F))

ma_rpfB <- ma_rpfB %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit,
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma_rpfB) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma_rpfB$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma_rpfB$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma_rpfB$log2FoldChange >= cap_limit, 24,
                            ifelse(ma_rpfB$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6","1e7")) +
  scale_y_continuous(limits=c(-cap_limit,cap_limit)) +
  theme_mycopore() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot - empty vs overexp")

fig_volcano <- ggplot(ma_rpfB) +
  geom_point(aes(x=log2FoldChange,
                 y=-log10(padj),
                 color=ifelse(padj<0.05,"grey60",colours$alex_R1),
                 fill=ifelse(padj<0.05,"grey60",colours$alex_R1))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
  theme_mycopore() +
  scale_color_manual(values = c("grey60",colours$alex_R3)) +
  theme(legend.position = "none") +
  xlab("log2-fold enrichment") +
  ylab("-log10 p(adjusted)") +
  ggtitle(paste0("RNAseq volcano plot"))

print(fig_volcano)

###### rpfB 0 vs 2 ######
# ma_rpfB <- as.data.frame(results(dds_compare, contrast = c("time","0h","2h"), cooksCutoff = F))
# 
# ma_rpfB <- ma_rpfB %>%
#   dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
#                                              sign(log2FoldChange) * cap_limit, log2FoldChange))
# 
# ggplot(ma_rpfB) +
#   geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
#              size = 0.8,
#              color = ifelse(ma_rpfB$padj < 0.05, colours$alex_R2, "grey60"),
#              fill = ifelse(ma_rpfB$padj < 0.05, colours$alex_R2, "grey60"),
#              shape = ifelse(ma_rpfB$log2FoldChange >= cap_limit, 24, 
#                             ifelse(ma_rpfB$log2FoldChange <= -1*cap_limit, 25, 21)),
#              alpha = 0.8) +
#   scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6","1e7")) +
#   scale_y_continuous(limits=c(-4,4)) +
#   theme_mycopore() +
#   geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
#   xlab("Mean of counts") +
#   ylab("log2-fold change") +
#   ggtitle("MA plot - rpfB - 0h and 2h")

###### 0h ######
# ma0 <- as.data.frame(results(dds_0h, contrast = c("type","Empty","rpfB"), cooksCutoff = F))
# 
# ma0 <- ma0 %>%
#   dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
#                                              sign(log2FoldChange) * cap_limit, log2FoldChange))
# 
# ggplot(ma0) +
#   geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
#              size = 0.8,
#              color = ifelse(ma0$padj < 0.05, colours$alex_R2, "grey60"),
#              fill = ifelse(ma0$padj < 0.05, colours$alex_R2, "grey60"),
#              shape = ifelse(ma0$log2FoldChange >= cap_limit, 24, 
#                             ifelse(ma0$log2FoldChange <= -1*cap_limit, 25, 21)),
#              alpha = 0.8) +
#   scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6","1e7")) +
#   scale_y_continuous(limits=c(-4,4)) +
#   theme_mycopore() +
#   geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
#   xlab("Mean of counts") +
#   ylab("log2-fold change") +
#   ggtitle("MA plot - empty vs rpfB, 0h")


##### So are any targets significantly changed? #####
RNA_targets <- read.csv(paste(here::here("Brindhaseqcounts/targetRNA3.csv"))) %>%
  mutate(real_name = ifelse(grepl("\\(MSMEG",Target),
                            sub(".\\(.*","\\1",Target),
                            gff2$realID[match(Target,gff2$weird_locus)]))

hits <- RNA_targets %>%
  filter(real_name %in% LFC_2h_table$rowname) %>%
  mutate(foldchange = LFC_2h_table$log2FoldChange[match(real_name,LFC_2h_table$rowname)],
         change_p = LFC_2h_table$padj[match(real_name,LFC_2h_table$rowname)])

hits_sig <- hits %>%
  filter(change_p < 0.05) %>%
  select(12,10,13,14,4,5) %>%
  dplyr::rename("locus" = 1,
                "RNAseq-change" = 3,
                "RNAseq-padjust" = 4,
                "sRNA-pvalue" = 5,
                "sRNA-probability" = 6)


fwrite(hits_sig, file=paste0(directory,"/sRNA_hits.csv"), row.names = T)



