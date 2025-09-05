###########################################################################

# Additional data analysis for D'Halluin et al., 2023

#   • Checking connections between "type" of peaks/frequency

###########################################################################

########
# Libs #
########

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","reshape2")
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

#### Get gff annotation ####
mtb_gff <- get_gff_input("Mtb") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
  mutate(Rv_name = str_split_fixed(str_split_fixed(attributes, ";Name=", 2)[,1], "Locus=", 2)[,2],) %>%
  mutate(product_name = str_split_fixed(str_split_fixed(attributes, ";Comments=", 2)[,1], "Product=", 2)[,2],) %>%
  mutate(function_group = str_split_fixed(str_split_fixed(attributes, ";Protein Data Bank", 2)[,1], "Functional_Category=", 2)[,2],) %>%
  mutate(npeaks = 0) %>%
  dplyr::filter(locus_name != "oxyR'") %>% #Misannotated pseudogene
  select(type,start,end,strand,attributes,locus_name, Rv_name, product_name, function_group, npeaks)

##### Split by strand #####
gff_plus <- mtb_gff %>%
  dplyr::filter(strand == "+") %>%
  arrange(start)

gff_minus <- mtb_gff %>%
  dplyr::filter(strand == "-") %>%
  arrange(end)

#### Get Termseq peaks ####

term_minus <- read_alex_bed("Peaks_ExpoMin") %>%
  arrange(main)
term_plus <- read_alex_bed("Peaks_ExpoPlus") %>%
  arrange(main)

#### Assign peaks to gff ####
##### Plus #####

for(gene in 1:nrow(gff_plus)){
  # Iterate through all genes
  n <- 0
  for(peak in 1:nrow(term_plus)){
    # Iterate through peaks
    k <- ifelse(term_plus$main[peak] > gff_plus$start[gene] & 
                  term_plus$main[peak] <= gff_plus$end[gene],1,0)
    n <- n+k
  }
  gff_plus$npeaks[gene] <- n
}

##### Minus #####
for(gene in 1:nrow(gff_minus)){
  # Iterate through all genes
  n <- 0
  for(peak in 1:nrow(term_minus)){
    # Iterate through peaks
    k <- ifelse(term_minus$main[peak] > gff_minus$start[gene] & 
                  term_minus$main[peak] <= gff_minus$end[gene],1,0)
    n <- n+k
  }
  gff_minus$npeaks[gene] <- n
}

##### Combine #####
gff_all <- rbind(gff_plus,gff_minus) %>%
  arrange(start) %>%
  mutate(length = end-start+1) %>%
  mutate(peakscore = (npeaks/length)*100) %>%
  arrange(desc(peakscore))

#### Plot ####
ggplot(gff_all,aes(peakscore)) +
  geom_histogram(binwidth = 0.01,fill=colours$purple) +
  theme_alex() +
  scale_x_continuous(expand=c(0,0), limits = c(0,5)) +
  scale_y_log10(expand=c(0,0)) +
  xlab("peak score (% of nt that is a TTS)") +
  ylab("log10 count") +
  ggtitle("Distribution of peak frequencies")

#### Export table ####
fwrite(gff_all, file = here::here("NUGA_Out/peakscore.csv"))





