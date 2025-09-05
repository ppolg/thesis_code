###########################################################################

# Additional data analysis for D'Halluin et al., 2023

# Checking the overlaps of Mtb genes:
#   • Read in FASTA, GFF
#   • Search ORFs for which overlap
#   • Filter out 4-nt
#   • Check for RBSs as 5 consecutive purine 5-15 nt upstream of TSS
#   • Make necessary tables/figures
#   • uORF analysis

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

#### Main ####

# Get lists for Rho-dependents (based in TTS-score) to check co-occurrence
#fdr_3 <- read.csv(here::here("NUGA_Data/TTS_3h.csv")) %>% filter(Class %in% c("I","O"))
#fdr_45 <- read.csv(here::here("NUGA_Data/TTS_45h.csv")) %>% filter(Class %in% c("I","O"))
#fdr_6 <- read.csv(here::here("NUGA_Data/TTS_6h.csv")) %>% filter(Class %in% c("I","O"))

# New - get list of all TTS and their TTS/RT Score

tts_threshold <- 1.1
rt_treshold <- 1.1
q_treshold <- 0.05

tts_rt <- read.csv(here::here("NUGA_Data/TTS_RT.csv")) %>%
  dplyr::rename("Name" = 5,
                "Type" = 6,
                "TTS3" = 13,
                "TTS45" = 14,
                "TTS6" = 15,
                "p_TTS3" = 16,
                "p_TTS45" = 17,
                "p_TTS6" = 18,
                "RT3" = 19,
                "RT45" = 20,
                "RT6" = 21,
                "p_RT3" = 22,
                "p_RT45" = 23,
                "p_RT6" = 24) %>%
  select(1:7,13:24) %>%
  filter(Class %in% c("I","O"))

fdr_3 <- tts_rt %>%
  filter((TTS3 >= tts_threshold &
           RT3 >= rt_treshold) |
           (p_TTS3 < q_treshold &
           p_RT3 < q_treshold))

fdr_45 <- tts_rt %>%
  filter((TTS45 >= tts_threshold &
            RT45 >= rt_treshold) |
           (p_TTS45 < q_treshold &
              p_RT45 < q_treshold))

fdr_6 <- tts_rt %>%
  filter((TTS6 >= tts_threshold &
            RT6 >= rt_treshold) |
           (p_TTS6 < q_treshold &
              p_RT6 < q_treshold))

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

gff_plus$seq_upstream <- str_replace_all(gff_plus$seq_upstream,replacements)

gff_plus <- gff_plus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

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

gff_minus$seq_upstream <- str_replace_all(gff_minus$seq_upstream,replacements)

gff_minus <- gff_minus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

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

# Get all overlaps together
overlaps_all <- rbind(plus_overlaps,minus_overlaps)

# Get remainder of n_overlap_next for colours(TM)
overlaps_all <- overlaps_all %>%
  mutate(frame = ifelse(is.na(n_overlap_next),0,n_overlap_next%%3))


# overlaps of 4 nt
overlaps_4 <- overlaps_all %>%
  dplyr::filter(n_overlap_prev==4)

# All TeRe
overlaps_tere <- overlaps_all %>%
  dplyr::filter(n_overlap_prev %in% c(1,4))

# Everything non-TeRe (overlap or not)
gff_not_tere <- gff_annot %>%
  filter(Rv_name %nin% overlaps_tere$Rv_name)

# For the pie chart of 4nt overlaps
overlaps_4_count <- overlaps_4 %>%
  dplyr::select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

# For the pie chart of 4nt overlaps, no RBS
overlaps_4_norbs_count <- overlaps_4 %>%
  dplyr::filter(has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

# For the pie chart of all start codons

gff_weirds <- gff_annot %>% # Select all that are non-NTG start
  filter(start_codon %nin% c("ATG","GTG","TTG","CTG"))

gff_noweirds <- gff_annot %>% # Select everyhting that IS NTG start
  filter(start %nin% gff_weirds$start)

start_codon_count <- gff_noweirds %>%
  dplyr::select(type,start,end,strand,locus_name,start_codon) %>%
  ungroup() %>%
  arrange(start_codon) %>%
  group_by(start_codon) %>%
  summarise(n = length(start_codon)) %>%
  add_row(start_codon = "other", n=nrow(gff_weirds)) %>%
  mutate(start_codon = factor(start_codon, levels=c("ATG","GTG","TTG","other")))


# Categorise overlap function - all overlaps
overlaps_func <- overlaps_all %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset= "overlap")

# Categorise overlap function - 4 nt
overlaps_func_4 <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4")

# Categorise overlap function - 1+4 nt
overlaps_func_tere <- overlaps_tere %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_tere")


# Categorise overlap function - 4 nt, no RBS
overlaps_func_4_noRBS <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "" & has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4_noRBS")

# Categorise overlap function - tere, no RBS
overlaps_func_tere_noRBS <- overlaps_tere %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "" & has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_tere_noRBS")

# Categorise overlap function - every ORF, ever
total_func <- mtb_gff %>%
  dplyr::filter( function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset="all")

# Combine
func_all <- rbind(overlaps_func,overlaps_func_4, overlaps_func_4_noRBS, total_func)
func_all2 <- rbind(overlaps_func,overlaps_func_tere, overlaps_func_tere_noRBS, total_func)

#### uORF ####

uORF_new <- read.csv(here::here("NUGA_Data/New_peptide.csv")) %>%
  dplyr::rename(start = NTG.start,
                end = NTG.end,
                strand = Strand,
                Rv_name = X5.UTR,
                type = Category) %>%
  mutate(Rv_name = ifelse(strand == "-" & str_sub(Rv_name, -1) != "c",
                           paste0(Rv_name,"c",sep=""),Rv_name)) %>%
  arrange(desc(strand),start)

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

uORF_new <- uORF_new %>%
  dplyr::mutate(locus_name = paste0(Rv_name,"_uORF",group)) %>%
  select(start,end,Rv_name,locus_name,strand,type)

##### Plus #####
uORF_plus <- uORF_new %>%
  filter(strand == "+") %>%
  mutate(prev_start = dplyr::lag(start, order_by=start),
         prev_end = dplyr::lag(end, order_by=start),
         next_start = dplyr::lead(start, order_by=start),
         next_end = dplyr::lead(end, order_by=start),
         overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
         overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(start, .by_group = TRUE)

# TEST for new overlap type
uORF_plus <- uORF_plus %>% dplyr::ungroup()
mains_list <- distinct(uORF_plus,Rv_name)
uORF_overlap_test_plus = tibble(strand = "+", uORF1 = 0, uORF2= 0, overlap_start = 0, overlap_end = 0)

for(row in 1:nrow(mains_list)){
  workspace <- uORF_plus %>%
    dplyr::filter(Rv_name == mains_list$Rv_name[row])
  for(i in 1:nrow(workspace)){
    for(k in i:nrow(workspace)){
      if(workspace$end[i] >= workspace$start[k] &
         workspace$end[i] < workspace$end[k]){
        add <- tibble(strand = "+",uORF1 = workspace$locus_name[i], uORF2 = workspace$locus_name[k],
                      overlap_start = workspace$start[k], overlap_end = workspace$end[i])
        uORF_overlap_test_plus <- rbind(uORF_overlap_test_plus,add)
      }
    }
  }
}


uORF_overlap_test_plus <- uORF_overlap_test_plus %>%
  filter(overlap_start != 0) %>%
  dplyr::mutate(n_overlap = overlap_end - overlap_start +1,
                seq_overlap = str_sub(mtb_fasta_string_plus,overlap_start,overlap_end))

# Overlap with main ORF
h <- gff_annot$start[match(x = uORF_plus$Rv_name, table=gff_annot$Rv_name)]
uORF_plus <- uORF_plus %>%
  add_column(main_5 = !!h)

h <- gff_annot$end[match(x = uORF_plus$Rv_name, table=gff_annot$Rv_name)]
uORF_plus <- uORF_plus %>%
  add_column(main_3 = !!h)

uORF_plus <- uORF_plus %>%
  mutate(overlap_main = ifelse(main_5 <= end & main_5 > start,TRUE,FALSE),
         n_overlap_main = ifelse(overlap_main == TRUE,end - main_5 +1,NA),
         #seq_overlap_main = ifelse(overlap_main == TRUE, toString(subseq(mtb_fasta,start=main_5,width=5)), NA)
  )

uORF_plus$seq_overlap_main <- NA # Get overlap seq, even though dplyr is not working for some reason :(
for(row in 1:nrow(uORF_plus)) {
  if(isTRUE(uORF_plus$overlap_main[row])){
    uORF_plus$seq_overlap_main[row] = toString(subseq(mtb_fasta,start=uORF_plus$main_5[row],end=uORF_plus$end[row]))
  }
}

##### Minus #####

uORF_minus <- uORF_new %>%
  filter(strand == "-") %>%
  arrange(end) %>%
  mutate(prev_start = dplyr::lead(start, order_by=end),
         prev_end = dplyr::lead(end, order_by=end),
         next_start = dplyr::lag(start, order_by=end),
         next_end = dplyr::lag(end, order_by=end),
         overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
         overlap_next = ifelse(start <= next_end & end > next_end,TRUE,FALSE )) %>%
  group_by(end) %>%
  arrange(desc(end))


# rename to fix previous issues
uORF_minus$group = NA
uORF_minus$group[1] = as.numeric(1)
for (row in 2:nrow(uORF_minus)) {
  if (uORF_minus$Rv_name[row] == uORF_minus$Rv_name[row-1]) {
    uORF_minus$group[row] <- uORF_minus$group[row-1] + 1
  }
  else {
    uORF_minus$group[row] <- as.numeric(1)
  }
}

uORF_minus <- uORF_minus %>%
  dplyr::mutate(locus_name = paste0(Rv_name,"_uORF",group)) %>%
  select(-group)

# Test for new overlap type
uORF_minus <- uORF_minus %>% dplyr::ungroup()
mains_list <- distinct(uORF_minus,Rv_name)
uORF_overlap_test_minus = tibble(strand="-",uORF1 = 0, uORF2= 0, overlap_start = 0, overlap_end = 0)

for(row in 1:nrow(mains_list)){
  workspace <- uORF_minus %>%
    dplyr::filter(Rv_name == mains_list$Rv_name[row])
  for(i in 1:nrow(workspace)){
    for(k in i:nrow(workspace)){
      if(workspace$start[i] <= workspace$end[k] &
         workspace$end[i] > workspace$end[k]){
        add <- tibble(strand="-",uORF1 = workspace$locus_name[i], uORF2 = workspace$locus_name[k],
                      overlap_start = workspace$start[i], overlap_end = workspace$end[k])
        uORF_overlap_test_minus <- rbind(uORF_overlap_test_minus,add)
      }
    }
  }
}

uORF_overlap_test_minus <- uORF_overlap_test_minus %>%
  filter(overlap_start != 0) %>%
  dplyr::mutate(n_overlap = overlap_end - overlap_start +1,
                seq_overlap = rc(str_sub(mtb_fasta_string_plus,overlap_start,overlap_end)))


# Overlap with main ORF
h <- gff_annot$end[match(x = uORF_minus$Rv_name, table=gff_annot$Rv_name)]
uORF_minus <- uORF_minus %>%
  add_column(main_5 = !!h)

h <- gff_annot$start[match(x = uORF_minus$Rv_name, table=gff_annot$Rv_name)]
uORF_minus <- uORF_minus %>%
  add_column(main_3 = !!h)


uORF_minus <- uORF_minus %>%
  mutate(overlap_main = ifelse(main_5 >= start & main_5 < end,TRUE,FALSE),
         n_overlap_main = ifelse(overlap_main == TRUE,main_5 - start +1,NA))


uORF_minus$seq_overlap_main <- NA
for(row in 1:nrow(uORF_minus)) {
  if(isTRUE(uORF_minus$overlap_main[row])){
    uORF_minus$seq_overlap_main[row] = toString(reverseComplement(subseq(mtb_fasta,start=uORF_minus$start[row],end=uORF_minus$main_5[row])))
  }
}

##### Combine #####
uORF_all <- rbind(uORF_plus,uORF_minus) %>%
  arrange(desc(strand),start)

uORF_test <- uORF_all %>%
  filter((strand=="+" & end %in% gff_plus$end) | (strand == "-" & start %in% gff_minus$start))

uORF_test2 <- uORF_all %>%
  filter((strand == "+" & start > main_3 | strand=="-" & end < main_3))

uORF_test3 <- uORF_all %>%
  filter(Rv_name %nin% gff_annot$Rv_name)

uORF_overlap_main <- uORF_all %>%
  filter(overlap_main == T) %>% 
  mutate(frame = ifelse(is.na(n_overlap_main),0,n_overlap_main%%3)) %>%
  filter(frame %in% c(1,2)) %>%
  arrange(n_overlap_main)

uORF_overlap_self <- rbind(uORF_overlap_test_plus,uORF_overlap_test_minus) %>%
  mutate(frame = ifelse(is.na(n_overlap),0,n_overlap%%3)) %>%
  filter(frame %in% c(1,2)) %>%
  arrange(n_overlap)

#### Testing the 4th nt in regards to NUGA ####

nt4_test <- gff_annot %>% #ALL NUGA overlaps
  dplyr::filter(overlap_prev == FALSE) %>%
  mutate(nt4 = toString(subseq(mtb_fasta,start+3,start+3)))

nt4_test_noSD <- nt4_test %>%
  filter(has_RBS == F)

nt4_has_uORF <- nt4_test %>%
  filter(Rv_name %in% uORF_all$Rv_name)

nt4_noSD_uORF <- nt4_test %>%
  filter(Rv_name %in% uORF_all$Rv_name & has_RBS == F)

nt4_sum <- nt4_test %>%
  group_by(nt4) %>%
  summarise(n=length(nt4)) %>%
  mutate(dataset = "all",
         n_mod = n/sum(n))

nt4_noSD_sum <- nt4_test_noSD %>%
  group_by(nt4) %>%
  summarise(n=length(nt4)) %>%
  mutate(dataset = "noSD",
         n_mod = n/sum(n))

nt4_uORF_sum <- nt4_has_uORF %>%
  group_by(nt4) %>%
  summarise(n=length(nt4)) %>%
  mutate(dataset = "has_uORF",
         n_mod = n/sum(n))

nt4_noSD_uORF_sum <- nt4_noSD_uORF %>%
  group_by(nt4) %>%
  summarise(n=length(nt4)) %>%
  mutate(dataset = "noSD_uORF",
         n_mod = n/sum(n))

nt4_all <- rbind(nt4_sum,nt4_noSD_sum,nt4_uORF_sum,nt4_noSD_uORF_sum)

#### TTS and TeRe ####

# How many NUGA has a conditional terminator? In all datasets?

uORF_tere <- uORF_overlap_main %>%
  filter(n_overlap_main %in% c(1,4))

uORF_non <- uORF_all %>%
  filter(locus_name %nin% uORF_tere$locus_name)

cooc_counts <- hash::hash(
  # annotated_annotated
  a_a_3 = overlaps_tere %>% filter(locus_name %in% fdr_3$Locus),
  n_a_a_3 = nrow(overlaps_tere %>% filter(locus_name %in% fdr_3$Locus)),
  a_a_45 = overlaps_tere %>% filter(locus_name %in% fdr_45$Locus),
  n_a_a_45 = nrow(overlaps_tere %>% filter(locus_name %in% fdr_45$Locus)),
  a_a_6 = overlaps_tere %>% filter(locus_name %in% fdr_6$Locus),
  n_a_a_6 = nrow(overlaps_tere %>% filter(locus_name %in% fdr_6$Locus)),
  # All the annotated non-NUGA ORFs
  a_n_3 = gff_not_tere %>% filter(Rv_name %in% fdr_3$Locus),
  n_a_n_3 = nrow(gff_not_tere %>% filter(Rv_name %in% fdr_3$Locus)),
  a_n_45 = gff_not_tere%>% filter(Rv_name %in% fdr_45$Locus),
  n_a_n_45 = nrow(gff_not_tere %>% filter(Rv_name %in% fdr_45$Locus)),
  a_n_6 = gff_not_tere %>% filter(Rv_name %in% fdr_6$Locus),
  n_a_n_6 = nrow(gff_not_tere %>% filter(Rv_name %in% fdr_6$Locus)),
  # uORF_annotated
  u_a_3 = uORF_tere %>% filter(Rv_name %in% fdr_3$Locus),
  n_u_a_3 = nrow(uORF_tere %>% filter(Rv_name %in% fdr_3)),
  u_a_45 = uORF_tere %>% filter(Rv_name %in% fdr_45$Locus),
  n_u_a_45 = nrow(uORF_tere %>% filter(Rv_name %in% fdr_45$Locus)),
  u_a_6 = uORF_tere %>% filter(Rv_name %in% fdr_6$Locus),
  n_u_a_6 = nrow(uORF_tere %>% filter(Rv_name %in% fdr_6$Locus)),
  # All the uORFs that do NOT NUGA
  u_n_3 = uORF_non %>% filter(Rv_name %in% fdr_3$Locus),
  n_u_n_3 = nrow(uORF_non %>% filter(Rv_name %in% fdr_3$Locus)),
  u_n_45 = uORF_non%>% filter(Rv_name %in% fdr_45$Locus),
  n_u_n_45 = nrow(uORF_non %>% filter(Rv_name %in% fdr_45$Locus)),
  u_n_6 = uORF_non %>% filter(Rv_name %in% fdr_6$Locus),
  n_u_n_6 = nrow(uORF_non %>% filter(Rv_name %in% fdr_6$Locus))
)

cooc_ratios <- tibble(name = c("a_a_3", "a_a_45","a_a_6",
                               "a_n_3","a_n_45","a_n_6",
                               "u_a_3","u_a_45","u_a_6",
                               "u_n_3","u_n_45","u_n_6"),
                      type = c("aa","aa","aa","an","an","an",
                               "ua","ua","ua","un","un","un"),
                      time = c("3h","4.5h","6h","3h","4.5h","6h",
                               "3h","4.5h","6h","3h","4.5h","6h"),
                      all_NUGA = c(nrow(overlaps_tere),nrow(overlaps_tere),nrow(overlaps_tere),
                                   nrow(gff_not_tere), nrow(gff_not_tere), nrow(gff_not_tere),
                                   nrow(uORF_tere), nrow(uORF_tere), nrow(uORF_tere),
                                   nrow(uORF_non), nrow(uORF_non), nrow(uORF_non)),
                      NUGA_with_TTS = c(cooc_counts$n_a_a_3,cooc_counts$n_a_a_45,
                                        cooc_counts$n_a_a_6,cooc_counts$n_a_n_3,
                                        cooc_counts$n_a_n_45,cooc_counts$n_a_n_6,
                                        cooc_counts$n_u_a_3,cooc_counts$n_u_a_45,
                                        cooc_counts$n_u_a_6,cooc_counts$n_u_n_3, 
                                        cooc_counts$n_u_n_45,cooc_counts$n_u_n_6))

cooc_ratios <- cooc_ratios %>% 
  mutate(NUGA_non_TTS = all_NUGA - NUGA_with_TTS,
         ratio = NUGA_with_TTS/all_NUGA)

#### Plot ####

# Overlap length distribution - annotated
ggplot(overlaps_all,aes(x = n_overlap_next)) +
  geom_histogram(data=subset(overlaps_all,frame==1),binwidth = 1, color = "black", fill = colours$alex_R3, size=1) +
  geom_histogram(data=subset(overlaps_all,frame==2),binwidth = 1, color = "black", fill = colours$alex_R1, size=1) +
  geom_histogram(data=subset(overlaps_all,frame%nin%c(1,2)),binwidth = 1, color = "black", fill = colours$alex_R2, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Count") +
  ggtitle("Annotated-annotated (count)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,600), breaks=c(0,150,300,450,600)) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# A-A: Percentage/frequency
overlaps_percent <- overlaps_all %>%
  group_by(frame,n_overlap_next) %>%
  summarize(n=length(n_overlap_next)) %>%
  ungroup() %>%
  mutate(n=n/sum(n)) %>%
  filter(frame != 0)

ggplot(overlaps_percent,aes(x = n_overlap_next, y=n*100)) +
  geom_col(data=subset(overlaps_percent,frame==1), width=1, color = "black", fill = colours$alex_B3, size=1) +
  geom_col(data=subset(overlaps_percent,frame==2), width=1,color = "black", fill = colours$alex_B1, size=1) +
  geom_col(data=subset(overlaps_percent,frame%nin%c(1,2)), color = "black", fill = colours$alex_B2, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Frequency") +
  ggtitle("Annotated-annotated overlaps") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,50), breaks=c(0,25,50), labels=c(0.00,0.25,0.5) ) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# Overlap length distribution - uORF-annotated
ggplot(uORF_overlap_main,aes(x = n_overlap_main)) +
  geom_histogram(data=subset(uORF_overlap_main,frame==1),binwidth = 1, color = "black", fill = colours$alex_B3, size=1) +
  geom_histogram(data=subset(uORF_overlap_main,frame==2),binwidth = 1, color = "black", fill = colours$alex_B1, size=1) +
  geom_histogram(data=subset(uORF_overlap_main,frame%nin%c(1,2)),binwidth = 1, color = "black", fill = colours$alex_B2, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Count") +
  ggtitle("uORF-annotated (count)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# u-A frequency
uORF_main_percent<- uORF_overlap_main %>%
  group_by(frame,n_overlap_main) %>%
  summarize(n=length(n_overlap_main)) %>%
  ungroup() %>%
  mutate(n=n/sum(n)) %>%
  filter(frame != 0)

ggplot(uORF_main_percent,aes(x = n_overlap_main, y=n*100)) +
  geom_col(data=subset(uORF_main_percent,frame==1), width=1, color = "black", fill = colours$green, size=1) +
  geom_col(data=subset(uORF_main_percent,frame==2), width=1,color = "black", fill = colours$lightlime2, size=1) +
  geom_col(data=subset(uORF_main_percent,frame%nin%c(1,2)), color = "black", fill = colours$lightlime, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Frequency") +
  ggtitle("uORF-annotated overlaps") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,50), breaks=c(0,25,50), labels=c(0.00,0.25,0.5) ) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# Overlap length distribution - uORF-uORF
ggplot(uORF_overlap_self,aes(x = n_overlap)) +
  geom_histogram(data=subset(uORF_overlap_self,frame==1),binwidth = 1, color = "black", fill = colours$green, size=1) +
  geom_histogram(data=subset(uORF_overlap_self,frame==2),binwidth = 1, color = "black", fill = colours$lightlime2, size=1) +
  geom_histogram(data=subset(uORF_overlap_self,frame%nin%c(1,2)),binwidth = 1, color = "black", fill = colours$lightlime, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Count") +
  ggtitle("uORF-uORF (count)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,20)) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

#u-u frequency
uORF_self_percent<- uORF_overlap_self %>%
  group_by(frame,n_overlap) %>%
  summarize(n=length(n_overlap)) %>%
  ungroup() %>%
  mutate(n=n/sum(n)) %>%
  filter(frame != 0)

ggplot(uORF_self_percent,aes(x = n_overlap, y=n*100)) +
  geom_col(data=subset(uORF_self_percent,frame==1), width=1, color = "black", fill = colours$alex_R3, size=1) +
  geom_col(data=subset(uORF_self_percent,frame==2), width=1,color = "black", fill = colours$alex_R1, size=1) +
  geom_col(data=subset(uORF_self_percent,frame%nin%c(1,2)), color = "black", fill = colours$lightlime, size=1) +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Frequency") +
  ggtitle("uORF-uORF overlaps") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,50), breaks=c(0,25,50), labels=c(0.00,0.25,0.5) ) +  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# Function categories distribution
ggplot(func_all, aes(x=function_group, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_alex() +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits=c(0,32),expand=c(0,0),breaks=c(0,10,20,30), labels = c("0%","10%","20%","30%")) +
  scale_fill_manual(values = c(colours$alex_B3,colours$alex_B2,colours$alex_B1,"grey80"), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

# Function categories distribution - merged TeRe
ggplot(func_all2, aes(x=function_group, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_alex() +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits=c(0,32),expand=c(0,0),breaks=c(0,10,20,30), labels = c("0%","10%","20%","30%")) +
  scale_fill_manual(values = c(colours$alex_B3,colours$alex_B2,colours$alex_B1,"grey80"), 
                    labels=c("All","Overlap","TeRe overlap","TeRe overlap, no SD"))

# nt4 distribution
ggplot(nt4_all, aes(x=nt4, y=n_mod*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity",width = 0.7) +
  theme_alex() +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits=c(0,40),expand=c(0,0),breaks=c(0,10,20,30,40), labels = c("0%","10%","20%","30%","40%")) +
  scale_fill_manual(values = c(colours$alex_R3,colours$alex_R2,colours$alex_R1,"grey80"),
                    labels=c("All","ORFs with uORF", "With no SD","With uORF no SD"))

# Cooc distribution
ggplot(cooc_ratios, aes(x=time, y=ratio*100, fill=type)) +
  geom_bar(position="dodge",stat="identity", width = 0.85) +
  theme_alex() +
  xlab("") +
  ylab("Coocurrence with Rho-TTS") +
  ggtitle("") +
  theme(legend.title = element_blank(),
        aspect.ratio = 1/1,
        legend.position = "bottom") +
  scale_x_discrete() +
  scale_y_continuous(limits=c(0,40),expand=c(0,0),breaks=c(0,10,20,30,40), labels = c("0%","10%","20%","30%","40%")) +
  scale_fill_manual(values = c(colours$blue,colours$alex_B1,colours$green,colours$lightlime2),
                               labels=c('Annotated-annotated',"Annotated-nonoverlap","uORF-annotated","uORF-nonoverlap"))


#### Stat/Out ####


##### Hypergeo #####

# Hypergeo for functional groups
func_annot <- as.data.frame(func_all) %>% filter(dataset == "all") %>% select(function_group,n)
func_overlap <- as.data.frame(func_all) %>% filter(dataset == "overlap") %>% select(function_group,n)

func_4nt <- func_all %>% filter(dataset == "overlap_4") %>% select(function_group,n)
func_norbs <- func_all %>% filter(dataset == "overlap_4_noRBS") %>% select(function_group,n)
func_4nt <- func_4nt %>% add_row (function_group = "unknown", n=0, .before = 9)
func_norbs <- func_norbs %>% add_row (function_group = "unknown", n=0, .before = 9)

# NEW for merged TeRe group!
func_tere <- func_all2 %>% filter(dataset == "overlap_tere") %>% select(function_group,n)
func_tere_norbs <- func_all2 %>% filter(dataset == "overlap_tere_noRBS") %>% select(function_group,n)
func_tere <- func_tere %>% add_row (function_group = "unknown", n=0, .before = 9)
func_tere_norbs <- func_tere_norbs %>% add_row (function_group = "unknown", n=0, .before = 9)

func_pvals <- func_annot %>% select(function_group)

p1 = c()
for(row in 1:nrow(func_overlap)){
  n = 1.0-phyper(func_overlap$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_overlap$n))
  p1 <- base::append(p1,n)
}

p2 = c()
for(row in 1:nrow(func_4nt)){
  n = 1.0-phyper(func_4nt$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_4nt$n))
  p2 <- base::append(p2,n)
}

p3 = c()
for(row in 1:nrow(func_norbs)){
  n = 1.0-phyper(func_norbs$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_norbs$n))
  p3 <- base::append(p3,n)
}

p4 = c()
for(row in 1:nrow(func_tere)){
  n = 1.0-phyper(func_tere$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_tere$n))
  p4 <- base::append(p4,n)
}

p5 = c()
for(row in 1:nrow(func_tere_norbs)){
  n = 1.0-phyper(func_tere_norbs$n[row]-1,func_annot$n[row],sum(func_annot$n)-func_annot$n[row],sum(func_tere_norbs$n))
  p5 <- base::append(p5,n)
}


##### FDR #####
# Both BH and BY, because twice as many numbers look twice as good!
p1_fdr <- p.adjust(p1,method = "BH")
p2_fdr <- p.adjust(p2,method = "BH")
p3_fdr <- p.adjust(p3,method = "BH")
p4_fdr <- p.adjust(p4,method = "BH")
p5_fdr <- p.adjust(p5,method = "BH")

p1_by <- p.adjust(p1,method = "BY")
p2_by <- p.adjust(p2,method = "BY")
p3_by <- p.adjust(p3,method = "BY")

func_pvals <- func_annot %>% select(function_group) %>%
  mutate(p_overlap := p1,
         p_overlap_4nt := p2,
         p_overlap_4bt_norbs := p3,
         p_overlap_BH := p1_fdr,
         p_overlap_4nt_BH := p2_fdr,
         p_overlap_4bt_norbs_BH := p3_fdr,
         p_overlap_BY := p1_by,
         p_overlap_4nt_BY := p2_by,
         p_overlap_4bt_norbs_BY := p3_by)

func_pvals_new <- func_annot %>% select(function_group) %>%
  mutate("padj between overlapping and all" := p1_fdr,
         "padj between TeRe overlap and all":= p4_fdr,
         "padj between TeRe (no RBS) and all" := p5_fdr)

##### Fisher #####

###### AA-UA ######
fish3_aaua <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[1],cooc_ratios$NUGA_non_TTS[7]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[1],cooc_ratios$NUGA_with_TTS[7]),
  row.names = c("A-A","u-A")
)

p_fish_3_aaua <- fisher.test(fish3_aaua)$p.value

fish45_aaua <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[2],cooc_ratios$NUGA_non_TTS[8]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[2],cooc_ratios$NUGA_with_TTS[8]),
  row.names = c("A-A","u-A")
)

p_fish_45_aaua <- fisher.test(fish45_aaua)$p.value

fish6_aaua <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[3],cooc_ratios$NUGA_non_TTS[9]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[3],cooc_ratios$NUGA_with_TTS[9]),
  row.names = c("A-A","u-A")
)

p_fish_6_aaua <- fisher.test(fish6_aaua)$p.value

###### AN-UN ######
fish3_anun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[4],cooc_ratios$NUGA_non_TTS[10]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[4],cooc_ratios$NUGA_with_TTS[10]),
  row.names = c("A-N","u-N")
)

p_fish_3_anun <- fisher.test(fish3_anun)$p.value

fish45_anun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[5],cooc_ratios$NUGA_non_TTS[11]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[5],cooc_ratios$NUGA_with_TTS[11]),
  row.names = c("A-N","u-N")
)

p_fish_45_anun <- fisher.test(fish45_anun)$p.value

fish6_anun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[6],cooc_ratios$NUGA_non_TTS[12]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[6],cooc_ratios$NUGA_with_TTS[12]),
  row.names = c("A-N","u-N")
)

p_fish_6_anun <- fisher.test(fish6_anun)$p.value

###### AA-AN ######
fish3_aaan <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[1],cooc_ratios$NUGA_non_TTS[4]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[1],cooc_ratios$NUGA_with_TTS[4]),
  row.names = c("A-A","A-N")
)

p_fish_3_aaan <- fisher.test(fish3_aaan)$p.value

fish45_aaan <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[2],cooc_ratios$NUGA_non_TTS[5]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[2],cooc_ratios$NUGA_with_TTS[5]),
  row.names = c("A-A","A-N")
)

p_fish_45_aaan <- fisher.test(fish45_aaan)$p.value

fish6_aaan <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[3],cooc_ratios$NUGA_non_TTS[6]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[3],cooc_ratios$NUGA_with_TTS[6]),
  row.names = c("A-A","A-N")
)

p_fish_6_aaan <- fisher.test(fish6_aaan)$p.value


###### UA-UN ######

fish3_uaun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[7],cooc_ratios$NUGA_non_TTS[10]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[7],cooc_ratios$NUGA_with_TTS[10]),
  row.names = c("u-A","u-N")
)

p_fish_3_uaun <- fisher.test(fish3_uaun)$p.value

fish45_uaun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[8],cooc_ratios$NUGA_non_TTS[11]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[8],cooc_ratios$NUGA_with_TTS[11]),
  row.names = c("u-A","u-N")
)

p_fish_45_uaun <- fisher.test(fish45_uaun)$p.value

fish6_uaun <- data.frame(
  no_TTS = c(cooc_ratios$NUGA_non_TTS[9],cooc_ratios$NUGA_non_TTS[12]),
  has_TTS = c(cooc_ratios$NUGA_with_TTS[9],cooc_ratios$NUGA_with_TTS[12]),
  row.names = c("u-A","u-N")
)

p_fish_6_uaun <- fisher.test(fish6_uaun)$p.value

#### PRINT ####

# Write lists of overlaps - all
print_overlaps <- overlaps_all %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,n_overlap_prev,seq_overlap_prev,
         overlap_next,n_overlap_next,seq_overlap_next,has_RBS) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlaps ORF downstream" = overlap_next,
                "Downstream overlap length" = n_overlap_next,
                "Sequence overlapping downstream" = seq_overlap_next,
                "Overlaps ORF upstream" = overlap_prev,
                "Upstream overlap length" = n_overlap_prev,
                "Sequence overlapping upstream" = seq_overlap_prev,
                "RBS present" = has_RBS)
fwrite(print_overlaps, file = here::here("NUGA_Out/table_overlaps_all.csv"))

# Write lists of overlaps - tere
print_overlaps_tere <- overlaps_all %>%
  filter(n_overlap_prev %in% c(1,4) | n_overlap_next %in% c(1,4)) %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,n_overlap_prev,
         seq_overlap_prev,overlap_next,n_overlap_next,seq_overlap_next,has_RBS) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlaps ORF downstream" = overlap_next,
                "Downstream overlap length" = n_overlap_next,
                "Sequence overlapping downstream" = seq_overlap_next,
                "Overlaps ORF upstream" = overlap_prev,
                "Upstream overlap length" = n_overlap_prev,
                "Sequence overlapping upstream" = seq_overlap_prev,
                "RBS present" = has_RBS)

fwrite(print_overlaps_tere, file = here::here("NUGA_Out/table_overlaps_tere.csv"))

# Write uORF lists for overlaps

print_uORF_self <- uORF_overlap_self %>%
  select(uORF1,uORF2,strand,overlap_start,overlap_end,n_overlap,seq_overlap) %>%
  dplyr::rename("First uORF" = uORF1,
                "Second uORF" = uORF2,
                "Strand" = strand,
                "First nt of overlap" = overlap_start,
                "Last nt of overlap" = overlap_end,
                "Length of overlap(nt)" = n_overlap,
                "Overlap_sequence" = seq_overlap)

print_uORF_main <- uORF_overlap_main %>%
  mutate(overlap_start = ifelse(strand == "+",main_5,start),
         overlap_end = ifelse(strand == "+", end,main_5)) %>%
  select(locus_name,Rv_name,strand,overlap_start,overlap_end,n_overlap_main,seq_overlap_main) %>%
  dplyr::rename("uORF" = locus_name,
                "Annotated ORF downstream" = Rv_name,
                "Strand" = strand,
                "First nt of overlap" = overlap_start,
                "Last nt of overlap" = overlap_end,
                "Overlap length" = n_overlap_main,
                "Overlap sequence" = seq_overlap_main)

fwrite(print_uORF_self, file = here::here("NUGA_Out/table_overlaps_uu.csv"))
fwrite(print_uORF_main, file = here::here("NUGA_Out/table_overlaps_uA.csv"))


# Write tables for pie charts to make later

fwrite(start_codon_count, file = here::here("NUGA_Out/pie_start_codons.csv"))
fwrite(overlaps_4_count, file = here::here("NUGA_Out/pie_4nt.csv"))
fwrite(overlaps_4_norbs_count, file = here::here("NUGA_Out/pie_noRBS.csv"))

# Write table for p values
fwrite(func_pvals, file = here::here("NUGA_Out/p_values.csv"))
fwrite(func_pvals_new, file = here::here("NUGA_Out/5E_p_values.csv"))

# Write cooc
cooc_print <- cooc_ratios %>%
  dplyr::select(type,time,all_NUGA,NUGA_with_TTS,ratio) %>%
  mutate(type = ifelse(type=="aa"
                         ,"Annotated-annotated",
                       ifelse(type=="an", "Annotated-nothing",
                              ifelse(type=="ua","uORF-annotated","uORF-nothing")))) %>%
  dplyr::rename("Overlap type" = type,
                "Timepoint" = time,
                "All NUGA overlaps" = all_NUGA,
                "NUGA overlaps with Rho-dependent TTS" = NUGA_with_TTS,
                "Ratio" = ratio)

fwrite(cooc_print, file = here::here("NUGA_Out/coocurrence.csv"))

# Write cooc p-values

cooc_pval <- data.frame(
  "comparison" = c("AA-UA","AN-UN","AA-AN","UA-UN"),
  "sample_1" = c("Annotated-annotated","Annotated-nonoverlap",
                 "Annotated-annotated","uORF-annotated"),
  "sample_2" = c("uORF-annotated","uORF-nonoverlap",
                 "Annotated-nonoverlap","uORF-nonoverlap"),
  "h3" = c(p_fish_3_aaua,p_fish_3_anun,p_fish_3_aaan,p_fish_3_uaun),
  "h45" = c(p_fish_45_aaua,p_fish_45_anun,p_fish_45_aaan,p_fish_45_uaun),
  "h6" = c(p_fish_6_aaua,p_fish_6_anun,p_fish_6_aaan,p_fish_6_uaun)
)

fwrite(cooc_pval, file = here::here("NUGA_Out/NUGA-pvalues.csv"))


# Write all
fwrite(overlaps_tere, file = here::here("NUGA_Out/NUGA_URAUG_overlaps_annotated.csv"))
fwrite(uORF_overlap_main, file = here::here("NUGA_Out/uORF_test_all.csv"))
