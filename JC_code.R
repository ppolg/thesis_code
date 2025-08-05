
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

# For checking RBS:
RBS_seq <- "RRRRR"   #  The strand of purines
RBS_range <- c(5,19) #  The range in which we look for purines

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
theme_alex <- function(base_size=20) {
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

# Read Mtb H37Rv FASTA (from NCBI)
mtb_fasta <- readDNAStringSet(here("NUGA_data/Mtb.fasta"))

mtb_fasta_string_plus <- toString(mtb_fasta)
mtb_fasta_string_minus <- toString(reverseComplement(mtb_fasta))

# Read Mtb H37Rv GFF (From Mycobrowser release 4)
mtb_gff <- get_gff_input("Mtb") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
  mutate(Rv_name = str_split_fixed(str_split_fixed(attributes, ";Name=", 2)[,1], "Locus=", 2)[,2],) %>%
  mutate(product_name = str_split_fixed(str_split_fixed(attributes, ";Comments=", 2)[,1], "Product=", 2)[,2],) %>%
  mutate(function_group = str_split_fixed(str_split_fixed(attributes, ";Protein Data Bank", 2)[,1], "Functional_Category=", 2)[,2],) %>%
  dplyr::filter(start != 1960667) %>% #Two ncRNAs with same TSS messing up my code
  dplyr::filter(locus_name != "oxyR'") %>% #Misannotated pseudogene
  filter(type == "CDS") %>%
  select(type,start,end,strand,attributes,locus_name, Rv_name, product_name, function_group)


##### PLUS STRAND #####

gff_plus <- mtb_gff %>%
  dplyr::filter(strand == "+") %>%
  dplyr::mutate(
    prev_start = dplyr::lag(start, order_by=start),
    prev_end = dplyr::lag(end, order_by=start),
    next_start = dplyr::lead(start, order_by=start),
    next_end = dplyr::lead(end, order_by=start),
    overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
    overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(start, .by_group = TRUE)

# Grab upstream for purines

gff_plus <- gff_plus %>%
  dplyr::mutate(seq_upstream = toString(subseq(mtb_fasta,start-RBS_range[2],start-RBS_range[1])))

gff_plus$seq_upstream2 <- str_replace_all(gff_plus$seq_upstream,replacements)

gff_plus <- gff_plus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream2,RBS_seq),TRUE,FALSE))

gff_plus <- gff_plus %>%
  filter(start > 5) %>%
  dplyr::mutate(nugug_test = toString(subseq(mtb_fasta,start-2,start+4)),
                has_augug = ifelse(str_detect(nugug_test,"ATGTG"),T,F),
                has_gugug = ifelse(str_detect(nugug_test,"GTGTG"),T,F))
# Check AUGUG/GUGUG
# Get overlap

plus_overlaps <- gff_plus %>%
  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
  dplyr::mutate(
    n_overlap_prev = ifelse(overlap_prev == TRUE,prev_end - start +1,NA),
    n_overlap_next = ifelse(overlap_next == TRUE,end - next_start +1,NA),
    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(subseq(mtb_fasta,start,prev_end)), NA),
    seq_overlap_next = ifelse(overlap_next == TRUE, toString(subseq(mtb_fasta,next_start,end)), NA)
  )

# Group into overlapping operons
#   i.e. if multiple overlaps, still same group

plus_overlaps$group = NA
plus_overlaps$group[1] = as.numeric(1)
for (row in 2:nrow(plus_overlaps)) {
  if (plus_overlaps$overlap_prev[row] == TRUE) {
    plus_overlaps$group[row] <- plus_overlaps$group[row-1]
  }
  else {
    plus_overlaps$group[row] <- plus_overlaps$group[row-1] + 1
  }
}





##### MINUS STRAND #####

gff_minus <- mtb_gff %>%
  dplyr::filter(strand == "-") %>%
  dplyr::mutate(
    prev_start = dplyr::lead(start, order_by=start),
    prev_end = dplyr::lead(end, order_by=start),
    next_start = dplyr::lag(start, order_by=start),
    next_end = dplyr::lag(end, order_by=start),
    overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
    overlap_next = ifelse(start <= next_end & start > next_start,TRUE,FALSE )) %>%
  group_by(start) %>%
  dplyr::arrange(desc(start))

# Grab range upstream

gff_minus <- gff_minus %>%
  dplyr::mutate(seq_upstream = toString(reverseComplement(subseq(mtb_fasta,end+RBS_range[1],end+RBS_range[2]))))

gff_minus$seq_upstream2 <- str_replace_all(gff_minus$seq_upstream,replacements)

gff_minus <- gff_minus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream2,RBS_seq),TRUE,FALSE))

# Check AUGUG/GUGUG
gff_minus <- gff_minus %>%
  dplyr::mutate(nugug_test = toString(reverseComplement(subseq(mtb_fasta,end-4,end+2))),
                has_augug = ifelse(str_detect(nugug_test,"ATGTG"),T,F),
                has_gugug = ifelse(str_detect(nugug_test,"GTGTG"),T,F))

# Get overlap

minus_overlaps <- gff_minus %>%
  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
  dplyr::mutate(
    n_overlap_prev = ifelse(overlap_prev == TRUE,end - prev_start +1,NA),
    n_overlap_next = ifelse(overlap_next == TRUE,next_end - start +1,NA),
    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(reverseComplement(subseq(mtb_fasta,prev_start,end))), NA),
    seq_overlap_next = ifelse(overlap_next == TRUE, toString(reverseComplement(subseq(mtb_fasta,start,next_end))), NA)
  )

# Group into overlapping operons
#   i.e. if multiple overlaps, still same group

minus_overlaps$group = NA
minus_overlaps$group[1] = as.numeric(plus_overlaps$group[nrow(plus_overlaps)]+1)
for (row in 2:nrow(minus_overlaps)) {
  if (minus_overlaps$overlap_prev[row] == TRUE) {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1]
  }
  else {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1] + 1
  }
}



##### COMBINED #####

#Get all gff together
gff_annot <- rbind(gff_plus,gff_minus) %>%
  mutate(start_codon = ifelse(strand == "+",toString(subseq(mtb_fasta,start,start+2)),toString(reverseComplement(subseq(mtb_fasta,end-2,end)))))



#Get all gff together
gff_annot_ol <- rbind(plus_overlaps,minus_overlaps) %>%
  mutate(start_codon = ifelse(strand == "+",toString(subseq(mtb_fasta,start,start+2)),toString(reverseComplement(subseq(mtb_fasta,end-2,end)))))

gff_nugug <- gff_annot %>%
  dplyr::filter(has_augug == T | has_gugug == T)

gff_nugug_ol <- gff_annot_ol %>%
  dplyr::filter(has_augug == T | has_gugug == T)

# Get uORf_ol
uORF_ol <- read.csv(here::here("NUGA_out/table_overlaps_uA.csv"))

gff_nugug_uORF_ol <- gff_nugug %>%
  filter(Rv_name %in% uORF_ol$Annotated.ORF.downstream)

# Get TTS_RT
TTS_RT <- read.csv(here::here("NUGA_Data/TTS_RT.csv"))

TTS_NUGUG <- TTS_RT %>%
  filter(Locus %in% gff_nugug$Rv_name)

gff_NUGUG_ol_TTS <- gff_nugug_ol %>%
  filter(Rv_name %in% TTS_NUGUG$Locus)
