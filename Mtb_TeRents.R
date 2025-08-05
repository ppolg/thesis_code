###########################################################################

# See nt distribution of overlapping ORFs

###########################################################################

########
# Libs #
########

packages <- c("insect","hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","stringi")
invisible(lapply(packages, require, character.only = TRUE))

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

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

########
# Init #
########


# For classifying nucleotides as pyrimidine/purine
replacements <- c("A"="R","G"="R","T"="Y","U"="Y","C"="Y" )

# Range of plot colours
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
                                     size = rel(2.8), hjust = 0.5,vjust=3),
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

#### Main ####

# Load exported overlaps
overlaps_all <- read.csv(here::here("NUGA_Out/table_overlaps_all.csv")) %>%
  filter(Overlaps.ORF.downstream == T) %>%
  arrange(Downstream.overlap.length) %>%
  dplyr::rename("locus" = 1,
                "product" = 2,
                "start" = 3,
                "end" = 4,
                "strand" = 5,
                "overlaps_prev" = 6,
                "n_overlap_prev" = 7,
                "seq_overlap_prev" = 8,
                "overlaps_next" = 9,
                "n_overlap_next" = 10,
                "seq_overlap_next" = 11,
                "has_RBS" = 12)

overlaps_nts <- overlaps_all %>%
  select(c(1,3:5,9:12)) %>%
  mutate(n_A = str_count(seq_overlap_next, "A"),
         n_T = str_count(seq_overlap_next, "T"),
         n_C = str_count(seq_overlap_next, "C"),
         n_G = str_count(seq_overlap_next, "G"))


overlap_nt_dist <- overlaps_nts %>%
  group_by(n_overlap_next) %>%
  summarise(n = length(locus),
            n_A = sum(n_A),
            n_T = sum(n_T),
            n_C = sum(n_C),
            n_G = sum(n_G)) %>%
  mutate(n_nts = n * n_overlap_next,
         "A" = n_A/n_nts,
         "T" = n_T/n_nts,
         "C" = n_C/n_nts,
         "G" = n_G/n_nts)

nt_melt <- reshape2::melt(overlap_nt_dist,id.vars=c("n_overlap_next"), measure.vars=c("A","T","C","G"))

ggplot(nt_melt, aes(x=n_overlap_next, y = value, fill=variable)) +
  geom_bar(stat="identity") +
  theme_alex() +
  scale_x_continuous(expand=c(0,0), limits = c(0,25), breaks = c(1,4,7,10,13,16,19,22,25)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("relative nt frequency") +
  xlab("overlap length (nt)") +
  scale_fill_manual(values = c("green3","red2","blue2","gold1")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("nt distribution of A-A overlaps")


