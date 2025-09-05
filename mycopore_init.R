
###########################################################################

# Global init for Mycopore scripts                                     
# Loads libraries (except "here", which is loaded in Init of other scripts) 
# Sets ggplot graphical settings for figures into a function           
# Defines all colour palettes for ggplot globally                      
# Defines any recurrent functions for ease of use     

###########################################################################


#----------------#
# LOAD LIBRARIES #
#________________#

packages <- c("hash","ggeconodist", "tidyverse", "here", "ggthemes", "gganimate",
              "writexl","data.table", "ggExtra", "Rsamtools", "GenomicAlignments",
              "UpSetR","seqTools", "Rsubread", "ape", "DT", "ggpubr", "ggridges",
              "ggsci","LncFinder","CoverageView", "gghalves", "pryr", "fst",
              "R.utils", "readxl", "ggseqlogo", "RColorBrewer","ggpattern","insect","RRHO", "waldo")

invisible(lapply(packages, require, character.only = TRUE))


#----------------#
# GGPLOT COLOURS #
#________________#

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
  alex_B1_edge = "#BBD0FA",
  nt_A = "#1ab317",
  nt_C = "#253edb",
  nt_G = "#ede02d",
  nt_T = "#d1171d"
)

#-----------------#
# GGPLOT SETTINGS #
#_________________#

# Global settings for plotting
theme_mycopore <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Helvetica")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_rect(colour = "black", size=rel(1.4)),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_text(angle=90,vjust = 3, size = rel(1.2)),
           axis.title.x = element_text(vjust = -0.2, size = rel(1.2)),
           axis.text = element_text(size = rel(1.3)), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           axis.ticks.length = unit(.2, "cm"),
           panel.grid.major = element_line(colour="grey80"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "bottom",
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


theme_heatmap <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2), hjust = 0.5),
           text = element_text(color = "black"),
           aspect.ratio = 1/1,
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_blank(),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_blank(),
           axis.title.x = element_text(vjust = -0.1, size = rel(1)),
           axis.text.x = element_text(face="bold",size = rel(1.3)), 
           axis.text.y = element_blank(),
           axis.line.x = element_line(colour="black",size=1),
           axis.line.y = element_blank(),
           axis.ticks.x = element_line(colour="black",size = 1),
           axis.ticks.y = element_blank(),
           axis.ticks.length = unit(.2, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "right",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "vertical",
           legend.key.size= unit(0.6, "cm"),
           legend.spacing.x = unit(0.2, "cm"),
           legend.text = element_text(color = "black", size = rel(0.9)),
           legend.title = element_blank(),
           plot.margin=unit(c(12,6,6,6),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}

# Guide legend function
guide_legend_heatmap <- function(){
  guide_colourbar(title = "Count",
                  title.position = "top",
                  label.theme = element_text(size = 10,
                                             family = "Arial",
                                             colour = "black",
                                             face = "bold"),
                  ticks = F,
                  nbin = 300,
                  barwidth = 1,
                  barheight = 7)
}

#------------------#
# SHARED FUNCTIONS #
#__________________#

# Read and tag guppy table
get_guppy_input <- function(table_file,name="RUN") {
  input_table <- paste(here::here("seq/R/data/"), table_file, ".txt", sep = "")
  table <- fread(input_table) %>%
    dplyr::mutate(run = name) # Add name of run as passed argument
}

# Read gff and optionally filter for type
get_gff_input <- function(gff_file,filter_type="none") {
  input_gff <- paste(here::here("seq/R/Data/"), gff_file, ".gff", sep = "") 
  gff <- read.gff(input_gff)
  
  # Filter if argument passed, to get only rRNA, tRNA, CDS etc...
  if(filter_type != "none"){ 
     gff %>% dplyr::filter(type == filter_type)
  }
  else{invisible(gff)}
}

# Read BAM file
#
# REVISE!!!!!!!!!!!

# o.O
get_bam_input <- function(bam_file) {
  input_bam <- paste(here::here("seq/R/data/"), bam_file, ".bam", sep = "")
  allReads <- readGAlignments(input_bam, use.names = T, param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), tag=c("NM"), what=c("mapq", "flag")))
  allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
    mutate(minion_read_name = names(allReads)) 
  
  # > get read sequences
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  res0 <- scanBam(input_bam,param = param)[[1]] # always list-of-lists
  allReads_sequence <- res0[["seq"]]                 # query widths
  allReads_sequence_table <- as.list(as.character(allReads_sequence))
  allReads_table$sequence <- unlist(allReads_sequence_table)
  allReads_table$n_char   <- nchar(allReads_table$sequence[1:length(allReads_table$sequence)])
  left  <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = allReads_table$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  allReads_table$soft_l <- as_tibble(cigarOpTable(left))$S
  allReads_table$hard_l <- as_tibble(cigarOpTable(left))$H
  allReads_table$soft_r <- as_tibble(cigarOpTable(right))$S
  allReads_table$hard_r <- as_tibble(cigarOpTable(right))$H
  return(allReads_table)
  
}

# nin
`%nin%` = Negate(`%in%`)

# get mode
get_mode <- function(x){
  return(names(sort(table(x), decreasing = T, na.last = T)[1]))
}


