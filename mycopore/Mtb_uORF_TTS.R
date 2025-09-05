###########################################################################

# Additional data analysis for D'Halluin et al., 2023

#   • Check cooccurrence of CondTTS and potential uORF

###########################################################################

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","reshape2","DESeq2")
invisible(lapply(packages, require, character.only = TRUE))

# ggplot2 theme for figures
theme_alex <- function(base_size=10) {
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

########
# Main #
########

#### uORF ####

# read list of predicted peptides in leaders
uORF_new <- read.csv(here::here("NUGA_Data/New_peptide.csv")) %>%
  dplyr::rename(start = NTG.start,
                end = NTG.end,
                strand = Strand,
                Rv_name = X5.UTR,
                type = Category) %>%
  # correct annotations
  mutate(Rv_name = ifelse(strand == "-" & str_sub(Rv_name, -1) != "c",
                          paste0(Rv_name,"c",sep=""),Rv_name)) %>%
  arrange(desc(strand),start)

# Remove duplicates
uORF_new <- unique(data.table(uORF_new), by = c('start','end'))

# Sort by leader
uORF_new$group = NA
uORF_new$group[1] = as.numeric(1)
for (row in 2:nrow(uORF_new)) {
  if (uORF_new$Rv_name[row] == uORF_new$Rv_name[row-1]) {
    uORF_new$group[row] <- uORF_new$group[row-1] + 1
  }
  else {
    uORF_new$group[row] <- as.numeric(1)
  }
}

# Set new names for easier identification
uORF_new <- uORF_new %>%
  dplyr::mutate(locus_name = paste0(Rv_name,"_uORF",group))

# Filter out so that only one entry per leader
uORF_leaders <- uORF_new %>% filter(group == 1) # 383 5' leaders that likely have uORF

#### UTR_TTS ####

# Read  list of conditional TTS
UTR_TTS <- read.csv(here::here("NUGA_data/condTTS_table.csv")) %>%
  filter(Class == 'O') %>%
  dplyr::rename(start = Position.Start,
                end = Position.End,
                strand = Strand,
                Rv_name = Locus,
                locus_name = Gene.name) %>%
  select(ID,start,end,strand,Rv_name,locus_name)

# The overlap between translated leaders and conditional TTS
TTS_in_transleader <- UTR_TTS %>%
  filter(Rv_name %in% uORF_leaders$Rv_name) #73

# Write the file
fwrite(TTS_in_transleader, file = here::here("NUGA_Out/TTS_uORF_cooc.csv"))

# Hypergeometric test for significance of enrichment
pval <- 1 - phyper(q = nrow(TTS_in_transleader),
               m = nrow(UTR_TTS),
               n = (1434 - nrow(UTR_TTS)),
               k = nrow(uORF_leaders)) #2.22044e-16

ratio1 <- 123/1434 #0.08574
ratio2 <- 73/383 #0.19060

# Ratio of ratios
ratio_all <- ratio2/ratio1 # 2.22212-fold increase between the two

