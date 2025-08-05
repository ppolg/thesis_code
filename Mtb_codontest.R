###########################################################################

# Additional data analysis for D'Halluin et al.: post-publication follow-up

#   • Tine's theory on wholesome wobbles 
#         - (CGA codon matching stop codon with wobble)

###########################################################################

########
# Libs #
########

invisible(library(here))
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))

########
# Func #
########

# Read gff and optionally filter for type
get_gff_input <- function(gff_file,filter_type="none") {
  input_gff <- paste(here::here("NUGA_data/"), gff_file, ".gff", sep = "") 
  gff <- read.gff(input_gff)
  
  # Filter if argument passed, to get only rRNA, tRNA, CDS etc...
  if(filter_type != "none"){ 
    gff %>% dplyr::filter(type == filter_type)
  }
  else{invisible(gff)}
}

# Function to read the bed(esque) files the TermSeq peak caller made
read_alex_bed <- function(file,strand = "+") {
  path <- paste(here::here("NUGA_data/"), file, ".bed", sep = "")
  bed <- as.tibble(read.table(path)) %>%
    select(2:3) %>%
    dplyr::rename(start = 1,
                  end = 2) %>%
    mutate(ID = row_number(),
           main = floor((start+end)/2)) %>%
    select(ID,main,start,end)
}

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

########
# Init #
########

colours <- hash::hash(
  red = "#A3280A",
  orange = "#E3812B",
  brown = "#8A5122",
  yellow = "#E0D253",
  grey = "#858482",
  green = "#195928",
  lime = "#83e649",
  lightlime = "#c5f081",
  lightlime2 = "#bed698",
  blue = "#4A749E",
  darkblue = "#18099c",
  azure = "#092780",
  purple = "#612882",
  pink = "#d64daf",
  darkgold = "#755807",
  alex_R1 = "#F0BFC1",
  alex_R2 = "#D46859",
  alex_R3 = "#871E12",
  alex_B1 = "#D7E0F7",
  alex_B2 = "#87A9D9",
  alex_B3 = "#5577AA",
  alex_B1_edge = "#BBD0FA"
)

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

# Read Mtb H37Rv FASTA (from NCBI)
mtb_fasta <- readDNAStringSet(here("NUGA_data/Mtb.fasta"))

mtb_fasta_string_plus <- toString(mtb_fasta)
mtb_fasta_string_minus <- toString(reverseComplement(mtb_fasta))

# Codon table for Mtb with rarities
mtb_codons <- read.csv(here::here("seq/R/Data/Mtb_codons.csv"))
mtb_codons$triplet <- str_replace_all(mtb_codons$triplet,"U","T")
# Read codon test of NUGA
codon_test <- read.csv(here::here("NUGA_OUT/codontests.csv"))
# Read test all
codon_test_all <- read.csv(here::here("NUGA_OUT/codons_all.csv"))
# Read codon test of all u-A
uORF_test_all <- read.csv(here::here("NUGA_Out/uORF_test_all.csv"))
uORF_test_tere <- uORF_test_all %>%
  filter(n_overlap_main == 4)

# Arginine NUGAs
arg_test <- codon_test %>%
  filter(aa10 == "R")

arg_sum <- arg_test %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  add_row(codon10 = "CGC", n = 0) %>%
  mutate(k = n/nrow(arg_test))


h <- mtb_codons$fraction[match(x = arg_sum$codon10, table=mtb_codons$triplet)]
arg_sum <- arg_sum %>%
  add_column(h = h) %>%
  dplyr::rename("NUGA last codon" = k,
         "All codons" = h)

arg_melt <- reshape2::melt(arg_sum, id.vars="codon10", measure.vars=c("NUGA last codon","All codons"))

ggplot(arg_melt, aes(x=codon10, y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_B3,colours$lightlime2)) +
  scale_x_discrete(expand=c(0.12,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6), breaks = c(0,0.2,0.4,0.6)) +
  xlab("Codon sequence") +
  ylab("Codon frequency") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("Arginine codon frequency")
                                                                                                                                                                        
# NGA last codon frequency
nga_test <- codon_test %>%
  filter(codon10  %in% c("CGA","GGA","TGA","AGA")) #TGA is stop. 91 out of 515

nga_test_all <- codon_test_all %>%
  filter(codon10  %in% c("CGA","GGA","TGA","AGA")) #TGA is stop. 91 out of 515

nga_sum <- nga_test %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  add_row(codon10 = "TGA", n = 0) %>%
  mutate(k = (n/nrow(codon_test))*1000)

h <- mtb_codons$per_thousand[match(x = nga_sum$codon10, table=mtb_codons$triplet)]
nga_sum <- nga_sum %>%
  add_column(h = h) %>%
  dplyr::rename("NUGA last codon" = k,
                "All codons" = h)

nga_melt <- reshape2::melt(nga_sum, id.vars="codon10", measure.vars=c("NUGA last codon","All codons"))

ggplot(nga_melt, aes(x=codon10, y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_B3,colours$lightlime2)) +
  scale_x_discrete(expand=c(0.3,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,130), breaks = c(0,50,100)) +
  xlab("Codon sequence") +
  ylab("Codons per thousand") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("NGA codon frequency")

# All last codon frequency

last_sum <- codon_test %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  mutate(k = (n/nrow(codon_test))*1000)

j <- mtb_codons %>%
  filter(triplet %nin% last_sum$codon10)


not_in_nuga <- mtb_codons %>%
  filter(triplet %nin% last_sum$codon10) %>%
  select(1) %>%
  rename(triplet = "codon10") %>%
  mutate(n = 0, k = 0)

last_sum <- rbind(last_sum,not_in_nuga) %>%
  arrange(codon10)

h <- mtb_codons$per_thousand[match(x = last_sum$codon10, table=mtb_codons$triplet)]
last_sum <- last_sum %>%
  add_column(h = h) %>%
  dplyr::rename("NUGA last codon" = k,
                "All codons" = h) %>%
  select(1,3,4)

last_melt <- reshape2::melt(last_sum, id.vars="codon10", measure.vars=c("NUGA last codon","All codons"))

ggplot(last_melt, aes(x=fct_rev(codon10), y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_B3,colours$lightlime2)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,130), breaks = c(0,50,100)) +
  xlab("Codon sequence") +
  ylab("Codons pet thousand") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  ggtitle("All frequency")

##### Check the last 10 codons of u-A! #####
# Grab + and -
uORF_test_plus <- uORF_test_tere %>%
  dplyr::filter(strand == "+") %>%
  select(locus_name,strand,start,end,main_5) %>%
  rowwise() %>%
  mutate(codon10 = toString(subseq(mtb_fasta,end-5,end-3)))

uORF_test_minus <- uORF_test_tere %>%
  dplyr::filter(strand == "-") %>%
  select(locus_name,strand,start,end,main_5) %>%
  rowwise() %>%
  mutate(codon10 = toString(reverseComplement(subseq(mtb_fasta,start+3,start+5))))

# Merge
uORF_test <- rbind(uORF_test_plus,uORF_test_minus)

# Get last n codons of overlap, STOP included

uORF_aa = mtb_codons$amino_acid[match(x = uORF_test$codon10, table=mtb_codons$triplet)]

uORF_test <- cbind(uORF_test,uORF_aa)

# Get aa
#for(i in 1:n_codons){ 
#  k <- paste0("codon",i)
#  m <- paste0("aa",i)
#  h <- mtb_codons$amino_acid[match(x = codon_test[[rlang::as_name(k)]], table=mtb_codons$triplet)]
#  codon_test <- codon_test %>%
#    add_column(!!m := h)
#}

#thing
uORF_last_sum <- uORF_test %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  mutate(k = (n/nrow(uORF_test))*1000)



uORF_not_in_nuga <- mtb_codons %>%
  filter(triplet %nin% uORF_last_sum$codon10) %>%
  select(1) %>%
  rename(triplet = "codon10") %>%
  mutate(n = 0, k = 0)

uORF_last_sum <- rbind(uORF_last_sum,uORF_not_in_nuga) %>%
  arrange(codon10)

h <- mtb_codons$per_thousand[match(x = uORF_last_sum$codon10, table=mtb_codons$triplet)]
uORF_last_sum <- uORF_last_sum %>%
  add_column(h = h) %>%
  dplyr::rename("uORF last" = k,
                "All codons" = h) %>%
  select(1,3,4)

uORF_last_melt <- reshape2::melt(uORF_last_sum, id.vars="codon10", measure.vars=c("uORF last","All codons"))

ggplot(uORF_last_melt, aes(x=fct_rev(codon10), y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_R3,colours$purple)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,130), breaks = c(0,50,100)) +
  xlab("Codon sequence") +
  ylab("Codons per thousand") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  ggtitle("uORF codon")



##### Composite plots #####

# Composite - Arginine codon frequency - all codons, all last, NUGA(annot) last

arg_test_all <- codon_test_all %>%
  filter(aa10 == "R")

arg_sum2 <- arg_test_all %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  mutate(k = n/nrow(arg_test_all))

arg_sum_final<- arg_sum %>%
  add_column("Last codons" = arg_sum2$k) %>%
  select(1,3,5,4)

arg_final_melt <- reshape2::melt(arg_sum_final, id.vars="codon10", measure.vars=c("NUGA last codon","Last codons","All codons"))

ggplot(arg_final_melt, aes(x=codon10, y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_R3,colours$alex_R2,"grey70")) +
  scale_x_discrete(expand=c(0.12,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6), breaks = c(0,0.2,0.4,0.6)) +
  xlab("Codon sequence") +
  ylab("Codon frequency") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("Arginine codon frequencies")


# Composite - NGA codon frequency - all codons, all last, NUGA(annot) last

nga_sum2 <- nga_test_all %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  add_row(codon10 = "TGA", n = 0) %>%
  mutate(k = (n/nrow(codon_test_all))*1000)

nga_sum_final<- nga_sum %>%
  add_column("Last codons" = nga_sum2$k) %>%
  select(1,3,5,4)

nga_final_melt <- reshape2::melt(nga_sum_final, id.vars="codon10", measure.vars=c("NUGA last codon","Last codons","All codons"))

ggplot(nga_final_melt, aes(x=codon10, y = value, fill=variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  theme_alex() +
  scale_fill_manual(values = c(colours$alex_R3,colours$alex_R2,"grey70")) +
  scale_x_discrete(expand=c(0.3,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,130), breaks = c(0,50,100)) +
  xlab("Codon sequence") +
  ylab("Codons per thousand") +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("NGA codon frequencies")



#### RER ####
# Evaluating "enrichment" - relative enrichment ratio score (RER)?

# • Show the ratio of NUGA overlap last codon against codon freq
#   • This already normalises to G/C content.
#   • Other things to normalise to:
#     • Log2: otherwise ratio is biased for lower occupancy
#     • No CUGA - we expect everything to increase due to not_in_nuga
#       • Define compared to the baseline??
#       • Expected increase due to all CUGA = sum(freq_CUGA)/sum(freq_non_CUGA)
#       • so baseline is this divided by 48.

# RER is: nuga_freq/baseline_freq, logged(?)

# For now, just a simple ratio

last_sum2 <- codon_test_all %>%
  group_by(codon10) %>%
  summarise(n = length(codon10)) %>%
  mutate(k = (n/nrow(codon_test_all))*1000)

j <- mtb_codons %>%
  filter(triplet %nin% last_sum2$codon10)


not_in_last <- mtb_codons %>%
  filter(triplet %nin% last_sum2$codon10) %>%
  select(1) %>%
  rename(triplet = "codon10") %>%
  mutate(n = 0, k = 0)

last_sum2 <- rbind(last_sum2,not_in_last) %>%
  arrange(codon10)

last_sum_final <- last_sum %>%
  add_column("Last codons" = last_sum2$k) %>%
  dplyr::rename(nuga = 2,
                bg = 3,
                last = 4) %>%
  dplyr::mutate(ratio = nuga/last,
                ratio_bg = nuga/bg,
                ratio_last = last/bg)

