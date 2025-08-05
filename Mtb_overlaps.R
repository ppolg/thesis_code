###########################################################################

# Reads in gff, fasta
# Pick ORF start/stop, get sequence
# Search for ORFs with:
#             • TSS is before previous TTS
#             • TTS is after previous TTS
# For each overlap, check/filter by length
# For filtered, select out NUGA 4-codon overlap
  
###########################################################################

########
# Func #
########

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


# Getting the next ORF for overlaps in uORFs

###########################################################################

########
# Init #
########

invisible(library(here))
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))

# Mtb FASTA read
mtb_fasta <- readDNAStringSet(here("seq/R/Data/Mtb.fasta"))

mtb_fasta_string_plus <- toString(mtb_fasta)
mtb_fasta_string_minus <- toString(reverseComplement(mtb_fasta))

# Mtb GFF read
mtb_gff <- get_gff_input("Mtb") %>%
  dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
  mutate(Rv_name = str_split_fixed(str_split_fixed(attributes, ";Name=", 2)[,1], "Locus=", 2)[,2],) %>%
  mutate(product_name = str_split_fixed(str_split_fixed(attributes, ";Comments=", 2)[,1], "Product=", 2)[,2],) %>%
  mutate(function_group = str_split_fixed(str_split_fixed(attributes, ";Protein Data Bank", 2)[,1], "Functional_Category=", 2)[,2],) %>%
  select(type,start,end,strand,attributes,locus_name, Rv_name, product_name, function_group)


# For classifying nucleotides as pyrimidine/purine
replacements <- c("A"="R","G"="R","T"="Y","U"="Y","C"="Y" )

# For checking RBS, the strand of purines
RBS_seq <- "RRRRR"

# For checking proline doublets
proline_seq <- "PP"

# Codon table for Mtb with rarities
mtb_codons <- read.csv(here::here("seq/R/Data/Mtb_codons.csv"))

# Leaderless predictions by Cortes et al. (2013)
ll_pred <- read.csv(here::here("seq/R/Data/leaderless_pred.csv")) %>%
  mutate_all(., list(~na_if(.,""))) %>%
  mutate(Leaderless = replace_na(Leaderless, "N"))

# Replace U with T to allow for direct reading from the coding strand DNA
mtb_codons$triplet <- str_replace_all(mtb_codons$triplet,"U","T")


# Tag rare, acidic, basic codons
mtb_codons <- mtb_codons %>%
  dplyr::mutate(is_rare = ifelse(per_thousand < 5,TRUE,FALSE),
                is_start = ifelse(triplet %in% c("ATG","GTG","TTG"),TRUE,FALSE),
                is_acidic = ifelse(amino_acid %in% c("D","E"),TRUE,FALSE),
                is_basic = ifelse(amino_acid %in% c("R","K","H"),TRUE,FALSE))


# For checking codons:
n_codons = 11

aa_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("aa",i)
  aa_positions=c(aa_positions,k)
}

rare_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("is_rare",i)
  rare_positions=c(rare_positions,k)
}

acid_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("is_acidic",i)
  acid_positions=c(acid_positions,k)
}

# "Name" of last codon for plotting later
last_aa <- paste0("aa",n_codons-1)
last_rare <- paste0("is_rare",n_codons-1)
last_acid <- paste0("is_acidic",n_codons-1)
###########################################################################

########
# Main #
########


# Get main GFF from Mycobrowser, tag ll from Cortes
h <- ll_pred$Leaderless[match(x = mtb_gff$Rv_name, table=ll_pred$RvNumber)]
mtb_gff <- mtb_gff %>%
  add_column(is_ll = !!h) %>%
  dplyr::filter(locus_name != "oxyR'") %>%
  filter(type == "CDS")

# Get uORF from Alex
uORF_gff <- read.csv(here::here("seq/R/Data/Alex_uORF.csv")) %>%
  dplyr::mutate(type = "uORF",
                attributes = NA,
                product_name = "uORF",
                function_group = "uORF") %>%
  mutate(is_ll = ifelse(Distance.between.SD.and.Start.codon == "Leaderless","L","N"),
         has_RBS = ifelse(is_ll == "L",F,T)) %>%
  arrange(desc(strand),start)

uORF_gff$group = NA
uORF_gff$group[1] = as.numeric(1)
for (row in 2:nrow(uORF_gff)) {
  if (uORF_gff$Rv_name[row] == uORF_gff$Rv_name[row-1]) {
    uORF_gff$group[row] <- uORF_gff$group[row-1] + 1
  }
  else {
    uORF_gff$group[row] <- as.numeric(1)
  }
}

uORF_gff <- uORF_gff %>%
  dplyr::mutate(locus_name = paste0(Rv_name,"_uORF",group)) %>%
  select(type,start,end,strand,attributes,
         locus_name,Rv_name,product_name,
         function_group,is_ll,has_RBS)

# Get 109 leaders from Alex

leaders_alex <- read.csv(here::here("seq/R/Data/RIBOSEQ_109.txt"),sep = "\t") %>%
  dplyr::distinct(UTR.start,UTR.end,strand,locus) %>%
  dplyr::rename(start = "UTR.start",
         end = "UTR.end")

# Get ORF data from Wade RiboRET:
wade_gff <- read.csv(here::here("seq/R/Data/Wade_ORF.csv")) %>%
  mutate_all(., list(~na_if(.,""))) %>%
  dplyr::rename(start = Start.Coordinate,
                end = Stop.Coordinate,
                strand = Strand,
                classification = Classification) %>%
  mutate(type = "Wade_pred") %>%
  select(type,start,end,strand,classification)

wade_gff <- wade_gff %>%
  mutate(buffer = start,
         start = ifelse(strand == "-",end,start),
         end = ifelse(strand =="-",buffer,end)) %>%
  select(type,start,end,strand,classification)

###########################################################################

# Annotated Dataset
# Get overlaps +/-
# Combine


#-------------#
# PLUS STRAND #
#-------------#

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

# Grab 5-19 nt upstream

gff_plus <- gff_plus %>%
  dplyr::filter(start != 1960667) %>% #Two ncRNAs with same TSS messing up my code
  dplyr::mutate(seq_upstream = toString(subseq(mtb_fasta,start-18,start-5)))

gff_plus$seq_upstream <- str_replace_all(gff_plus$seq_upstream,replacements)

gff_plus <- gff_plus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Grab 5-19nt upstream of next
gff_plus <- gff_plus %>%
  dplyr::mutate(seq_next_upstream = toString(subseq(mtb_fasta,next_start-18,next_start-5)))

gff_plus$seq_next_upstream <- str_replace_all(gff_plus$seq_next_upstream,replacements)

gff_plus <- gff_plus %>%
  dplyr::mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))

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



#--------------#
# MINUS STRAND #
#--------------#

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

# Grab 5-19 nt upstream

gff_minus <- gff_minus %>%
  dplyr::mutate(seq_upstream = toString(reverseComplement(subseq(mtb_fasta,end+5,end+19))))

gff_minus$seq_upstream <- str_replace_all(gff_minus$seq_upstream,replacements)

gff_minus <- gff_minus %>%
  dplyr::mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Grab 5-19 nt upstream of next

gff_minus <- gff_minus %>%
  dplyr::mutate(seq_next_upstream = ifelse(Rv_name != "Rv0008c",toString(reverseComplement(subseq(mtb_fasta,next_end+5,next_end+19))),NA))

gff_minus$seq_next_upstream <- str_replace_all(gff_minus$seq_next_upstream,replacements)

gff_minus <- gff_minus %>%
  dplyr::mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))


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
minus_overlaps$group[1] = as.numeric(plus_overlaps$group[nrow(plus_overlaps)])
for (row in 2:nrow(minus_overlaps)) {
  if (minus_overlaps$overlap_prev[row] == TRUE) {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1]
  }
  else {
    minus_overlaps$group[row] <- minus_overlaps$group[row-1] + 1
  }
}


#----------#
# COMBINED #
#----------#


gff_annot <- rbind(gff_plus,gff_minus) %>%
  mutate(start_codon = ifelse(strand == "+",toString(subseq(mtb_fasta,start,start+2)),toString(reverseComplement(subseq(mtb_fasta,end-2,end)))))
overlaps_all <- rbind(plus_overlaps,minus_overlaps)

gff_weirds <- gff_annot %>%
  filter(start_codon %nin% c("ATG","GTG","TTG","CTG"))

gff_noweirds <- gff_annot %>%
  filter(start %nin% gff_weirds$start)

overlaps_25nt <- overlaps_all %>%
  filter(n_overlap_next <26 & !is.na(n_overlap_next))
#  filter(!is.na(is_ll))

overlaps_3 <- overlaps_all %>%
  filter(overlap_next == T)


# overlaps of 4 nt
overlaps_4 <- overlaps_all %>%
  dplyr::filter(n_overlap_prev==4)

overlaps_4_3 <-overlaps_all %>%
  dplyr::filter(n_overlap_next==4)

# Only non-leaderless (for RBS%)
gff_non_ll <- gff_annot %>%
  dplyr::filter(is_ll == "N",)

###########################################################################

# Alex's Predictions
# Overlaps of uORF between each other
# Overlaps of uORF with actual ORFs



#Split +/-, get start of main annotated ORF from Alex

########
# Plus #
########

uORF_plus <- uORF_gff %>%
  filter(strand == "+") %>%
  mutate(start = start+1,
         end = end+3,
         prev_start = dplyr::lag(start, order_by=start),
         prev_end = dplyr::lag(end, order_by=start),
         next_start = dplyr::lead(start, order_by=start),
         next_end = dplyr::lead(end, order_by=start),
         overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
         overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  mutate(start = ifelse(is_ll == "L",start-1,start)) %>%
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

#uORF_overlaps_plus <- uORF_plus %>%
#  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
#  dplyr::mutate(
#    n_overlap_prev = ifelse(overlap_prev == TRUE,prev_end - start +1,NA),
#    n_overlap_next = ifelse(overlap_next == TRUE,end - next_start +1,NA),
#    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(subseq(mtb_fasta,start,prev_end)), NA),
#    seq_overlap_next = ifelse(overlap_next == TRUE, toString(subseq(mtb_fasta,next_start,end)), NA)
#  )

#########
# Minus #
#########

uORF_minus <- uORF_gff %>%
  filter(strand == "-") %>%
  arrange(end) %>%
  mutate(start = start-2,
         prev_start = dplyr::lead(start, order_by=end),
         prev_end = dplyr::lead(end, order_by=end),
         next_start = dplyr::lag(start, order_by=end),
         next_end = dplyr::lag(end, order_by=end),
         overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
         overlap_next = ifelse(start <= next_end & end > next_end,TRUE,FALSE )) %>%
  mutate(start = ifelse(is_ll == "L",start-1,start)) %>%
  group_by(end) %>%
  arrange(desc(end))


# rename cuz I'm bad
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

# TEST for new overlap type
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

uORF_minus <- uORF_minus %>%
  mutate(overlap_main = ifelse(main_5 >= start & main_5 < end,TRUE,FALSE),
         n_overlap_main = ifelse(overlap_main == TRUE,main_5 - start +1,NA))


uORF_minus$seq_overlap_main <- NA
for(row in 1:nrow(uORF_minus)) {
  if(isTRUE(uORF_minus$overlap_main[row])){
    uORF_minus$seq_overlap_main[row] = toString(reverseComplement(subseq(mtb_fasta,start=uORF_minus$start[row],end=uORF_minus$main_5[row])))
  }
}

#uORF_overlaps_minus <- uORF_minus %>%
#  dplyr::filter(overlap_next == TRUE | overlap_prev == TRUE ) %>%
#  dplyr::mutate(
#    n_overlap_prev = ifelse(overlap_prev == TRUE,end - prev_start +1,NA),
#    n_overlap_next = ifelse(overlap_next == TRUE,next_end - start +1,NA),
#    seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(reverseComplement(subseq(mtb_fasta,prev_start,end))), NA),
#    seq_overlap_next = ifelse(overlap_next == TRUE, toString(reverseComplement(subseq(mtb_fasta,start,next_end))), NA)
#  )

###########
# Combine #
###########

uORF_all <- rbind(uORF_plus,uORF_minus) %>%
  arrange(desc(strand),start)

#For testing - misannotations, yo!
uORF_all <- uORF_all %>%
  mutate(codons = ((end-start+1)/3)-1)
  #       start_codod = ifelse(strand == "+",toString(subseq(mtb_fasta,start,start+3)),toString(reverseComplement(subseq(mtb_fasta,start,start+3)))))



uORF_test <- uORF_all %>%
  filter(start %in% gff_annot$start | end %in% gff_annot$end) 

uORF_overlap_test_all <- rbind(uORF_overlap_test_plus,uORF_overlap_test_minus)

#uORF_overlaps_all <- rbind(uORF_overlaps_plus,uORF_overlaps_minus) %>%
#  arrange(desc(strand),start)
###########################################################################

# Wade's Data
# Compare with gff from mycobrowser, and with Alex's data
#     • See how much of an overlap there is
#     • Select all actually novel
#     • Test overlaps


# Get all new
wade_test <- wade_gff %>%
  filter(start %in% gff_annot$start & end %in% gff_annot$end)

wade_new <- wade_gff %>%
  filter(start %nin% gff_annot$start & end %nin% gff_annot$end)

wade_new2 <- wade_gff %>%
  filter(classification == "Novel")

wade_newdiff <- wade_new2 %>%
  filter(start %nin% wade_new$start)

#wade_test2 <- wade_gff %>%
#  filter(classification == "Annotated")

#wade_test3 <- wade_test2 %>%
#  filter(start %nin% wade_test$start)

#wade_test4 <- wade_test %>%
#  filter(start %nin% wade_test2$start)



# Get the wades in Alex's 109 leaders
wade_in_alex <- tibble("start" = 0,"end" = 0,"strand" = NA,"classification" = NA, "locus" = NA)

for(waderow in 1:nrow(wade_gff)){
  for(alexrow in 1:nrow(leaders_alex)){
    if(wade_gff$start[waderow] >= leaders_alex$start[alexrow] & 
       wade_gff$end[waderow] <= leaders_alex$end[alexrow]){
      add <- tibble("start" = wade_gff$start[waderow],
                    "end" = wade_gff$end[waderow],
                    "strand" = wade_gff$strand[waderow],
                    "classification" = wade_gff$classification[waderow],
                    "locus" = leaders_alex$locus[alexrow])
      wade_in_alex <- rbind(wade_in_alex,add)
    }
  }
}

wade_in_alex <- wade_in_alex %>%
  dplyr::filter(start != 0)

###########################################################################
# Overlap length distribution

# Annot-Annot
ggplot(overlaps_all,aes(n_overlap_next)) +
  geom_histogram(binwidth = 1, color="grey40", fill=colours$blue) +
  theme_mycopore() +
  xlab("Overlap length (nt)") +
  ylab("Count (log10)") +
  ggtitle("Length of overlaps (nt)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_log10(expand = c(0,0), breaks = c(0,10,100,1000), limits=c(1,1e3)) +
  scale_x_continuous(limits = c(0,70), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))
  
# uORF_main
uORF_all %>% dplyr::filter(overlap_main == T) %>%
  ggplot(aes(n_overlap_main)) +
  geom_histogram(binwidth = 1, color="grey40", fill=colours$blue) +
  theme_mycopore() +
  xlab("Overlap length (nt)") +
  ylab("Count") +
  ggtitle("Length of uORF overlaps (nt)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,10), limits=c(0,20)) +
  scale_x_continuous(limits = c(0,65), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))

# uORF_uORF

###########################################################################
# NTGA distribution of 4-overlaps (different to overall start codon or not?)

overlaps_4_count <- overlaps_4 %>%
  dplyr::select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

ggplot(overlaps_4_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
  geom_bar(width = 1,stat = "identity") +
  geom_text(aes(x=1.3, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_manual(values = c("grey20","grey20","grey20")) +
  scale_fill_manual(values = rev(c(colours$orange,colours$red,colours$brown))) +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  ggtitle("4nt-overlap sequence") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))

# NTGA distribution of non_RBS 4-overlaps

overlaps_4_norbs_count <- overlaps_4 %>%
  dplyr::filter(has_RBS == FALSE) %>%
  select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(seq_overlap_prev) %>%
  group_by(seq_overlap_prev) %>%
  summarise(n = length(seq_overlap_prev))

ggplot(overlaps_4_norbs_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
  geom_bar(width = 1,stat = "identity") +
  geom_text(aes(x=1.3, y = c(236,72,5), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_manual(values = c("grey20","grey20","grey20")) +
  scale_fill_manual(values = rev(c(colours$orange,colours$red,colours$brown))) +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  ggtitle("4nt-overlap sequence (no RBS)") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))

# NTG distribution of all start codons
start_codon_count <- gff_noweirds %>%
  dplyr::select(type,start,end,strand,locus_name,start_codon) %>%
  ungroup() %>%
  arrange(start_codon) %>%
  group_by(start_codon) %>%
  summarise(n = length(start_codon)) %>%
  add_row(start_codon = "other", n=nrow(gff_weirds)) %>%
  mutate(start_codon = factor(start_codon, levels=c("ATG","GTG","TTG","other")))

ggplot(start_codon_count,aes(x="",y=n, colour=start_codon, fill=start_codon))+
  geom_bar(width = 1,stat = "identity") +
  geom_text(aes(x=1.3, y = c(3000,1000,190,1), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_manual(values = c("grey20","grey20","grey20","grey20")) +
  scale_fill_manual(values = rev(c(colours$yellow,colours$orange,colours$red,colours$brown))) +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  ggtitle("All start codons") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))


###########################################################################
#CDS vs ncRNA etc.

overlaps_type <- overlaps_all %>%
  dplyr::filter(overlap_next == TRUE) %>%
  select(type,start,end,strand,locus_name,seq_overlap_next) %>%
  ungroup() %>%
  arrange(type) %>%
  group_by(type) %>%
  summarise(n = length(type))


###########################################################################
# Function group distribution

# Overlaps per function_group:

overlaps_func <- overlaps_all %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset= "overlap")

ggplot(overlaps_func,aes(x="",y=n, colour=function_group, fill=function_group))+
  geom_bar(width = 1,stat = "identity") +
#  geom_text(aes(x=1.25, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_grey(start=0.4, end=0.4)+
  scale_fill_brewer(palette="Set3") +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank())

# Only 4-laps
overlaps_func_4 <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4")

ggplot(overlaps_func_4,aes(x="",y=n, colour=function_group, fill=function_group))+
  geom_bar(width = 1,stat = "identity") +
  #  geom_text(aes(x=1.25, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_grey(start=0.4, end=0.4)+
  scale_fill_brewer(palette="Set3") +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank())

# 4lap-noRBS
overlaps_func_4_noRBS <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE & function_group != "" & has_RBS == FALSE & is_ll == "N") %>%
  select(type,start,end,strand,locus_name, function_group,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset = "overlap_4_noRBS")

ggplot(overlaps_func_4_noRBS,aes(x="",y=n, colour=function_group, fill=function_group))+
  geom_bar(width = 1,stat = "identity") +
  #  geom_text(aes(x=1.25, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_grey(start=0.4, end=0.4)+
  scale_fill_brewer(palette="Set3") +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank())

#Everything (for control)
total_func <- mtb_gff %>%
  dplyr::filter( function_group != "") %>%
  select(type,start,end,strand,locus_name, function_group) %>%
  ungroup() %>%
  arrange(function_group) %>%
  group_by(function_group) %>%
  summarise(n = length(function_group)) %>%
  mutate(n_norm = n/sum(n), dataset="all")

ggplot(total_func,aes(x="",y=n, colour=function_group, fill=function_group))+
  geom_bar(width = 1,stat = "identity") +
  #  geom_text(aes(x=1.25, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_mycopore() +
  scale_colour_grey(start=0.4, end=0.4)+
  scale_fill_brewer(palette="Set3") +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank())


# Get normalised distribution of func_groups, box_plot

func_all <- rbind(overlaps_func,overlaps_func_4, overlaps_func_4_noRBS, total_func)

ggplot(func_all, aes(x=function_group, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_mycopore() +
  coord_flip() +
  xlab("") +
  ylab("") +
  ggtitle("ORF functions") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks=c(10,20,30), labels = c("10%","20%","30%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

###########################################################################
# RBS for overlapped sites (all and 4-only?)

overlaps_rbs <- overlaps_all %>%
  dplyr::filter(overlap_prev == TRUE) %>%
  select(type,start,end,strand,locus_name, has_RBS,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(has_RBS) %>%
  group_by(has_RBS) %>%
  summarise(n = length(has_RBS)) %>%
  mutate(n_norm = n/sum(n), dataset= "overlap")

overlaps_4_rbs <- overlaps_4 %>%
  dplyr::filter(overlap_prev == TRUE) %>%
  select(type,start,end,strand,locus_name, has_RBS,seq_overlap_prev) %>%
  ungroup() %>%
  arrange(has_RBS) %>%
  group_by(has_RBS) %>%
  summarise(n = length(has_RBS)) %>%
  mutate(n_norm = n/sum(n), dataset= "overlap_4")

total_rbs <- gff_annot %>%
  select(type,start,end,strand,locus_name, has_RBS) %>%
  ungroup() %>%
  arrange(has_RBS) %>%
  group_by(has_RBS) %>%
  summarise(n = length(has_RBS)) %>%
  mutate(n_norm = n/sum(n), dataset= "all")

non_ll_rbs <- gff_non_ll %>%
  select(type,start,end,strand,locus_name, has_RBS) %>%
  ungroup() %>%
  arrange(has_RBS) %>%
  group_by(has_RBS) %>%
  summarise(n = length(has_RBS)) %>%
  mutate(n_norm = n/sum(n), dataset= "non_ll")

#Combine, plot
rbs_all <- rbind(overlaps_rbs,overlaps_4_rbs,non_ll_rbs, total_rbs) %>%
  dplyr::filter (has_RBS == TRUE)

ggplot(rbs_all, aes(x=has_RBS, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_mycopore() +
  xlab("") +
  ylab("") +
  ggtitle("Percentage of ORFs with RBS") +
  scale_y_continuous(breaks=c(10,20,30,40), labels = c("10%","20%","30%","40%")) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank()) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue),
                    labels = c("All","Non-leaderless","Overlaps","4nt-overlaps"))
  

# Overlap length vs RBS count

lengths_vs_RBS <- overlaps_all %>%
  dplyr::filter(overlap_prev == TRUE) %>%
  select(type,start,end,strand,locus_name, has_RBS,seq_overlap_prev,n_overlap_prev) %>%
  ungroup() %>%
  arrange(n_overlap_prev) %>%
  group_by(n_overlap_prev) %>%
  summarise(RBS = (sum(has_RBS == TRUE))/length(has_RBS))

ggplot(lengths_vs_RBS, aes(x=n_overlap_prev,y=RBS*100)) +
  geom_point(size=2) +
  scale_x_continuous(limits = c(0,72), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,70)) +
  theme_mycopore() +
  xlab("Overlap length (nt)") +
  ylab("Percentage with RBS") +
  ggtitle("RBS% against overlap length") +
  scale_y_continuous(breaks=c(0,25,50,75,100), labels = c("0%","25%","50%","75%","100%"))
  


###########################################################################
# How often does a non-overlap_prev, non-SD ORF has and NTGA motif at start?
# for finding unannotated

nt4a_test <- gff_annot %>% #Non-SD NUGA overlaps
  dplyr::filter(overlap_prev == FALSE & has_RBS == FALSE) %>%
  mutate(nt4 = toString(subseq(mtb_fasta,start+3,start+3))) %>%
  filter (nt4 == "A")

nt4a_test2 <- gff_annot %>% #ALL NUGA overlaps
  dplyr::filter(overlap_prev == FALSE) %>%
  mutate(nt4 = toString(subseq(mtb_fasta,start+3,start+3))) %>%
  filter (nt4 == "A")


###########################################################################
# Check last 10(?) codons of ORF ovelapped_next
#     • Get sequence as 10 codons
#           -> count rarity out of 10
#
#     • Translate
#           -> count acidic/basic/proline?


# Grab + and -
codon_test_plus <- gff_annot %>%
  dplyr::filter(strand == "+" & type == "CDS") %>%
  select(locus_name,strand,start,end, overlap_next, next_has_RBS)

codon_test_minus <- gff_annot %>%
  dplyr::filter(strand == "-" & type == "CDS") %>%
  select(locus_name,strand,start,end, overlap_next, next_has_RBS)


# Get last n codons of overlap, STOP included

for(i in 1:n_codons){    # + strand
  k <- paste0("codon",i)
  codon_test_plus <- codon_test_plus %>%
    dplyr::mutate(!!k := toString(subseq(mtb_fasta,end-3*(n_codons-i+1)+1,end-3*(n_codons-i+1)+3)))
}


for(i in 1:n_codons){    # - strand
  k <- paste0("codon",i)
  codon_test_minus <- codon_test_minus %>%
    dplyr::mutate(!!k := toString(reverseComplement(subseq(mtb_fasta,start+3*(n_codons-i),start+3*(n_codons-i)+2))))
}

# Merge
codon_test <- rbind(codon_test_plus,codon_test_minus)

# Get aa
for(i in 1:n_codons){ 
  k <- paste0("codon",i)
  m <- paste0("aa",i)
  h <- mtb_codons$amino_acid[match(x = codon_test[[rlang::as_name(k)]], table=mtb_codons$triplet)]
  codon_test <- codon_test %>%
    add_column(!!m := h)
}

# Get codon rarity
for(i in 1:n_codons){ 
  k <- paste0("codon",i)
  m <- paste0("is_rare",i)
  h <- mtb_codons$is_rare[match(x = codon_test[[rlang::as_name(k)]], table=mtb_codons$triplet)]
  codon_test <- codon_test %>%
    add_column(!!m := h)
}

codon_test <- codon_test %>%
  rowwise() %>%
  mutate(n_rare = sum(c_across(starts_with("is_rare")))-1) # -1 for stop codon



# Get codon acidity
for(i in 1:n_codons){ 
  k <- paste0("codon",i)
  m <- paste0("is_acidic",i)
  h <- mtb_codons$is_acidic[match(x = codon_test[[rlang::as_name(k)]], table=mtb_codons$triplet)]
  codon_test <- codon_test %>%
    add_column(!!m := h)
}

codon_test <- codon_test %>%
  rowwise() %>%
  mutate(n_acidic = sum(c_across(starts_with("is_acidic"))))

# Get codon basicity
for(i in 1:n_codons){ 
  k <- paste0("codon",i)
  m <- paste0("is_basic",i)
  h <- mtb_codons$is_basic[match(x = codon_test[[rlang::as_name(k)]], table=mtb_codons$triplet)]
  codon_test <- codon_test %>%
    add_column(!!m := h)
}

codon_test <- codon_test %>%
  rowwise() %>%
  mutate(n_basic = sum(c_across(starts_with("is_basic"))))

#Get proline count
codon_test <- codon_test %>%
  rowwise() %>%
  mutate(n_proline = sum(c_across(!!aa_positions) == "P"))

# Get aa string
codon_test <- codon_test %>%
  unite(col = "aa_full",!!aa_positions,sep="",remove = FALSE, na.rm=TRUE)

# Check aa string for double prolines
codon_test <- codon_test %>%
  dplyr::mutate(has_P_stretch = ifelse(str_detect(aa_full,proline_seq),TRUE,FALSE))

# Select overlap only
codon_test_overlap <- codon_test %>%
  dplyr::filter(locus_name %in% overlaps_all$locus_name & overlap_next == TRUE)

# Select 4-overlap only
codon_test_overlap_4 <- codon_test_overlap %>%
  dplyr::filter(locus_name %in% overlaps_4_3$locus_name)

# Select 4-overlap, no-RBS only
codon_test_overlap_4_noRBS <- codon_test_overlap_4 %>%
  dplyr::filter(next_has_RBS == FALSE)

########
# RARE #
########

# Check rare codon counts for plotting
rare_codons_all <- codon_test %>%
  dplyr::filter(n_rare >= 0) %>%
  ungroup() %>%
  arrange(n_rare) %>%
  group_by(n_rare) %>%
  summarise(count = length(n_rare)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "all")

rare_codons_overlap <- codon_test_overlap %>%
  dplyr::filter(n_rare >= 0) %>%
  ungroup() %>%
  arrange(n_rare) %>%
  group_by(n_rare) %>%
  summarise(count = length(n_rare)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")

rare_codons_overlap_4 <- codon_test_overlap_4 %>%
  dplyr::filter(n_rare >= 0) %>%
  ungroup() %>%
  arrange(n_rare) %>%
  group_by(n_rare) %>%
  summarise(count = length(n_rare)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")

rare_codons_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  dplyr::filter(n_rare >= 0) %>%
  ungroup() %>%
  arrange(n_rare) %>%
  group_by(n_rare) %>%
  summarise(count = length(n_rare)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")

rare_codons_sum = rbind(rare_codons_all,rare_codons_overlap,rare_codons_overlap_4,rare_codons_overlap_4_noRBS)

ggplot(rare_codons_sum, aes(x=n_rare, y=count_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_mycopore() +
  xlab("Number of rare codons") +
  ylab("Relative abundance") +
  ggtitle("Number of rare codons in last 10") +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(breaks = c(0,20,40), labels = c("0%","20%","40%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

# For rare codon per position
# Get rare composition for all
rare_test <- codon_test %>%
  select(starts_with("is_rare")) %>%
  gather(position,rare) %>%
  group_by(position,rare) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-rare),funs(. / nrow(codon_test))) %>%
  select(rare,!!rare_positions) %>%
  filter(rare == TRUE) %>%
  mutate(category = "all")

rare_test_overlap <- codon_test_overlap %>%
  select(starts_with("is_rare")) %>%
  gather(position,rare) %>%
  group_by(position,rare) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap))) %>%
  select(rare,!!rare_positions) %>%
  filter(rare == TRUE) %>%
  mutate(category = "overlap")

rare_test_overlap_4 <- codon_test_overlap_4 %>%
  select(starts_with("is_rare")) %>%
  gather(position,rare) %>%
  group_by(position,rare) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4))) %>%
  select(rare,!!rare_positions) %>%
  filter(rare == TRUE) %>%
  mutate(category = "overlap_4")

rare_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  select(starts_with("is_rare")) %>%
  gather(position,rare) %>%
  group_by(position,rare) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
  select(rare,!!rare_positions) %>%
  filter(rare == TRUE) %>%
  mutate(category = "overlap_4_noRBS")

rare_test_all <- rbind(rare_test,rare_test_overlap,rare_test_overlap_4,rare_test_overlap_4_noRBS) %>%
  select(category, !!rare_positions)

rare_test_all %>%
  gather(pos,value,is_rare1:!!last_rare) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_mycopore() +
  xlab("Codon position") +
  ylab("Percentage of rare codons") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,5,10,15,20), labels = c("0%","5%","10%","15%","20%")) +
  ggtitle("Rare codons in last 10 codons") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))


########
# ACID #
########

# Check acidic codons for plotting
acid_codons_all <- codon_test %>%
  ungroup() %>%
  arrange(n_acidic) %>%
  group_by(n_acidic) %>%
  summarise(count = length(n_acidic)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "all")

acid_codons_overlap <- codon_test_overlap %>%
  ungroup() %>%
  arrange(n_acidic) %>%
  group_by(n_acidic) %>%
  summarise(count = length(n_acidic)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")

acid_codons_overlap_4 <- codon_test_overlap_4 %>%
  ungroup() %>%
  arrange(n_acidic) %>%
  group_by(n_acidic) %>%
  summarise(count = length(n_acidic)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")

acid_codons_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  ungroup() %>%
  arrange(n_acidic) %>%
  group_by(n_acidic) %>%
  summarise(count = length(n_acidic)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")

acid_codons_sum = rbind(acid_codons_all,acid_codons_overlap,acid_codons_overlap_4,acid_codons_overlap_4_noRBS)

ggplot(acid_codons_sum, aes(x=n_acidic, y=count_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_mycopore() +
  xlab("Number of acidic amino acids") +
  ylab("Relative abundance") +
  ggtitle("Number of acidic codons in last 10") +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9)) +
  scale_y_continuous(breaks = c(0,10,20,30,40), labels = c("0%","10%","20%","30%","40%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))


# For acid codon per position
# Get acid composition for all
acid_test <- codon_test %>%
  select(starts_with("is_acidic")) %>%
  gather(position,acid) %>%
  group_by(position,acid) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-acid),funs(. / nrow(codon_test))) %>%
  select(acid,!!acid_positions) %>%
  filter(acid == TRUE) %>%
  mutate(category = "all")

acid_test_overlap <- codon_test_overlap %>%
  select(starts_with("is_acid")) %>%
  gather(position,acid) %>%
  group_by(position,acid) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap))) %>%
  select(acid,!!acid_positions) %>%
  filter(acid == TRUE) %>%
  mutate(category = "overlap")

acid_test_overlap_4 <- codon_test_overlap_4 %>%
  select(starts_with("is_acid")) %>%
  gather(position,acid) %>%
  group_by(position,acid) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap_4))) %>%
  select(acid,!!acid_positions) %>%
  filter(acid == TRUE) %>%
  mutate(category = "overlap_4")

acid_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  select(starts_with("is_acid")) %>%
  gather(position,acid) %>%
  group_by(position,acid) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
  select(acid,!!acid_positions) %>%
  filter(acid == TRUE) %>%
  mutate(category = "overlap_4_noRBS")

acid_test_all <- rbind(acid_test,acid_test_overlap,acid_test_overlap_4,acid_test_overlap_4_noRBS) %>%
  select(category, !!acid_positions)

acid_test_all %>%
  gather(pos,value,is_acidic1:!!last_acid) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_mycopore() +
  xlab("Codon position") +
  ylab("Percentage of acidic codons") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,5,10,15,20),labels = c("0%","5%","10%","15%","20%")) +
  ggtitle("Acidic codons in last 10 codons") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))


###########
# PROLINE #
###########

proline_codons_all <- codon_test %>%
  ungroup() %>%
  arrange(n_proline) %>%
  group_by(n_proline) %>%
  summarise(count = length(n_proline)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "all")

proline_codons_overlap <- codon_test_overlap %>%
  ungroup() %>%
  arrange(n_proline) %>%
  group_by(n_proline) %>%
  summarise(count = length(n_proline)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")

proline_codons_overlap_4 <- codon_test_overlap_4 %>%
  ungroup() %>%
  arrange(n_proline) %>%
  group_by(n_proline) %>%
  summarise(count = length(n_proline)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")

proline_codons_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  ungroup() %>%
  arrange(n_proline) %>%
  group_by(n_proline) %>%
  summarise(count = length(n_proline)) %>%
  mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")

proline_codons_sum = rbind(proline_codons_all,proline_codons_overlap,proline_codons_overlap_4,proline_codons_overlap_4_noRBS)

ggplot(proline_codons_sum, aes(x=n_proline, y=count_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_mycopore() +
  xlab("Number of prolines") +
  ylab("Relative abundance") +
  ggtitle("Number of prolines in last 10") +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9)) +
  scale_y_continuous(breaks = c(0,10,20,30,40,50), labels = c("0%","10%","20%","30%","40%","50%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

# Do they have proline double?
proline_stretch_all <- codon_test %>%
  ungroup() %>%
  arrange(has_P_stretch) %>%
  group_by(has_P_stretch) %>%
  summarise(count = length(has_P_stretch)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "all")

proline_stretch_overlap <- codon_test_overlap %>%
  ungroup() %>%
  arrange(has_P_stretch) %>%
  group_by(has_P_stretch) %>%
  summarise(count = length(has_P_stretch)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "overlap")

proline_stretch_overlap_4 <- codon_test_overlap_4 %>%
  ungroup() %>%
  arrange(has_P_stretch) %>%
  group_by(has_P_stretch) %>%
  summarise(count = length(has_P_stretch)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "overlap_4")

proline_stretch_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  ungroup() %>%
  arrange(has_P_stretch) %>%
  group_by(has_P_stretch) %>%
  summarise(count = length(has_P_stretch)) %>%
  mutate(count_norm = count/nrow(codon_test), dataset = "overlap_4_noRBS")

proline_stretch_sum <- rbind(proline_stretch_all,
                             proline_stretch_overlap,
                             proline_stretch_overlap_4,
                             proline_stretch_overlap_4_noRBS) %>%
  dplyr::filter(has_P_stretch == T)

ggplot(proline_stretch_sum, aes(x=has_P_stretch, y=count_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_mycopore() +
  xlab("") +
  ylab("") +
  ggtitle("Proline doublets in last 10 aa") +
  scale_y_continuous(breaks=c(1,2,3,4,5), labels = c("1%","2%","3%","4%","5%")) +
  theme(legend.title = element_blank(),
        axis.text.y = element_blank()) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue),
                    labels = c("All","Overlap","4nt-overlap","4nt-overlap, no RBS"))


########
# N_AA #
########

# Get aa composition at all
aa_test_all <- codon_test %>%
  select(starts_with("aa")) %>%
  gather(position,amino) %>%
  group_by(position,amino) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-amino),funs(. / nrow(codon_test))) %>%
  select(amino,!!aa_positions)

# Get aa composition at overlap
aa_test_overlap <- codon_test_overlap %>%
  select(starts_with("aa")) %>%
  gather(position,amino) %>%
  group_by(position,amino) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap))) %>%
  select(amino,!!aa_positions)

# Get aa composition at 4nt-overlap
aa_test_overlap_4 <- codon_test_overlap_4 %>%
  select(starts_with("aa")) %>%
  gather(position,amino) %>%
  group_by(position,amino) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4))) %>%
  select(amino,!!aa_positions)

# Get aa composition at 4nt-overlap, no RBS
aa_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
  select(starts_with("aa")) %>%
  gather(position,amino) %>%
  group_by(position,amino) %>%
  summarise(no = n()) %>%
  spread(position,no) %>%
  mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
  select(amino,!!aa_positions)

# Long format and plot
aa_test_all %>%
  gather(pos,value,aa1:!!last_aa) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=amino)) +
  geom_bar(stat="identity") +
  theme_mycopore() +
  xlab("Amino acid position") +
  ylab("Percentage occurrence") +
  ggtitle("Last 10 aa of all ORFs") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,25,50,75,100), labels = c("0%","25%","50%","75%","100%")) +
  theme(legend.title = element_blank())

aa_test_overlap %>%
  gather(pos,value,aa1:!!last_aa) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=amino)) +
  geom_bar(stat="identity") +
  theme_mycopore() +
  xlab("Amino acid position") +
  ylab("Percentage occurrence") +
  ggtitle("Last 10 aa of all overlapping ORFs") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,25,50,75,100), labels = c("0%","25%","50%","75%","100%")) +
  theme(legend.title = element_blank())

aa_test_overlap_4 %>%
  gather(pos,value,aa1:!!last_aa) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=amino)) +
  geom_bar(stat="identity") +
  theme_mycopore() +
  xlab("Amino acid position") +
  ylab("Percentage occurrence") +
  ggtitle("Last 10 aa of 4nt-overlaps") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,25,50,75,100), labels = c("0%","25%","50%","75%","100%")) +
  theme(legend.title = element_blank())

aa_test_overlap_4_noRBS %>%
  gather(pos,value,aa1:!!last_aa) %>%
  ggplot(aes(x = fct_inorder(pos), y=value*100, fill=amino)) +
  geom_bar(stat="identity") +
  theme_mycopore() +
  xlab("Amino acid position") +
  ylab("Percentage occurrence") +
  ggtitle("Last 10 aa of 4nt, no RBS") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_y_continuous(breaks = c(0,25,50,75,100), labels = c("0%","25%","50%","75%","100%")) +
  theme(legend.title = element_blank())

##########
# STR_AA #
##########

#All all
ggplot() +
  geom_logo(codon_test$aa_full) +
  theme_mycopore() +
  theme(legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Amino acid position") +
  ylab("Bits") +
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  ggtitle("Last 10 aa of all ORFs")

#Overlap all
ggplot() +
  geom_logo(codon_test_overlap$aa_full) +
  theme_mycopore() +
  theme(legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Amino acid position") +
  ylab("Bits") +
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  ggtitle("Last 10 aa of all overlaps")


#Overlap4
ggplot() +
  geom_logo(codon_test_overlap_4$aa_full) +
  theme_mycopore() +
  theme(legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Amino acid position") +
  ylab("Bits") +
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  ggtitle("Last 10 aa of 4nt overlaps")

#Overlap4, noRRBS
ggplot() +
  geom_logo(codon_test_overlap_4_noRBS$aa_full) +
  theme_mycopore() +
  theme(legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank()) +
  xlab("Amino acid position") +
  ylab("Bits") +
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  ggtitle("Last 10 aa of 4nt non-RBS overlaps")


###########################################################################
# Make printable tables and print them

# • overlaps (all,4)
# • nt4a
# • uORF (ovelaps_main, overlaps_self, the 27 club)

# Overlaps all
print_overlaps <- overlaps_all %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,n_overlap_prev,seq_overlap_prev,
                overlap_next,n_overlap_next,seq_overlap_next,group) %>%
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
                "Group" = group)

# Overlaps 4 
print_overlaps_4 <- overlaps_all %>%
  filter(n_overlap_prev == 4 | n_overlap_next == 4) %>%
  select(locus_name,product_name,start,end,strand,overlap_prev,seq_overlap_prev,
         overlap_next,seq_overlap_next,group) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlaps ORF downstream" = overlap_next,
                "Sequence overlapping downstream" = seq_overlap_next,
                "Overlaps ORF upstream" = overlap_prev,
                "Sequence overlapping upstream" = seq_overlap_prev,
                "Group" = group)


# Everything that is 4th nt A, non-overlap, non-SD (e.g. could be unannotated upstream)
print_nt4a <- nt4a_test %>%
  select(locus_name,product_name,start,end,strand) %>%
  dplyr::rename("Locus name" = locus_name,
                "Product" = product_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand)

# uORFs - all that overlaps main
print_uORF_main <- uORF_all %>%
  filter(overlap_main == T) %>%
  select(Rv_name,start,end,strand,n_overlap_main,seq_overlap_main) %>%
  dplyr::rename("Annotated ORF downstream" = Rv_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand,
                "Overlap length" = n_overlap_main,
                "Overlap sequence" = seq_overlap_main)

# uORFs - all that overlap each other
print_uORF_overlaps <- uORF_overlap_test_all %>%
  dplyr::rename("uORF 1" = uORF1,
                "uORF 2" = uORF2,
                "Overlap start" = overlap_start,
                "Overlap End" = overlap_end,
                "Strand" = strand,
                "Overlap length" = n_overlap,
                "Sequence overlapping" = seq_overlap)

# uORFs - the 27 club
print_uORF_test <- uORF_test %>%
  select(Rv_name,start,end,strand) %>%
  dplyr::rename("Annotated ORF downstream" = Rv_name,
                "Start" = start,
                "End" = end,
                "Strand" = strand)

tinetest <- print_uORF_main %>%
  dplyr::filter(`Overlap length` == 4)

tine_asked <- gff_annot %>%
  dplyr::filter(Rv_name %in% tinetest$`Annotated ORF downstream`,
                has_RBS == F)

###########################################################################
# Alex figs

#Length distribution (fig a)
ggplot(overlaps_all,aes(n_overlap_next)) +
  geom_histogram(binwidth = 1, color="black", size=1, fill="grey80") +
  theme_alex() +
  xlab("Overlap length (nt)") +
  ylab("Count (log10)") +
  scale_linetype_manual(values = c(4,1)) +
  scale_y_log10(expand = c(0,0), breaks = c(1,10,100,1000), limits=c(1,1e3)) +
  scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))

# Pie charts (fig b) - 4nt
ggplot(overlaps_4_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
  geom_bar(width = 1,stat = "identity",size=1) +
  geom_text(aes(x=1.3, y = c(380,123,11), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_alex() +
  scale_colour_manual(values = c("black","black","black")) +
  scale_fill_manual(values = rev(c("grey40","grey60","grey80"))) +
  scale_linetype_manual(values = c(4,1)) +
  ggtitle("all 4nt") +
  xlab("") +
  ylab("") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)),
        axis.ticks = element_blank(), axis.line = element_blank())

# Pie charts (fig b) - noRBS
ggplot(overlaps_4_norbs_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
  geom_bar(width = 1,stat = "identity", size=1) +
  geom_text(aes(x=1.3, y = c(231,72,5), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_alex() +
  scale_colour_manual(values = c("black","black","black")) +
  scale_fill_manual(values = rev(c("grey40","grey60","grey80"))) +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  ggtitle("4nt-noRBS") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)),
        axis.ticks = element_blank(), axis.line = element_blank())

# Pie Charts (fig b) - start codons
ggplot(start_codon_count,aes(x="",y=n, colour=start_codon, fill=start_codon))+
  geom_bar(width = 1,stat = "identity", size=1) +
#  geom_text(aes(x=1.3, y = c(3000,1000,170,1), label = n), size=5, show.legend = FALSE) +
  coord_polar("y", start = 0) +
  theme_alex() +
  scale_colour_manual(values = c("black","black","black","black")) +
  scale_fill_manual(values = rev(c("grey20","grey40","grey60","grey80"))) +
  scale_linetype_manual(values = c(4,1)) +
  xlab("") +
  ylab("") +
  ggtitle("All start codons") +
  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)),
        axis.ticks = element_blank(), axis.line = element_blank())

# Categories (fig c)
ggplot(func_all, aes(x=function_group, y=n_norm*100, fill=dataset)) +
  geom_bar(position="dodge",stat="identity") +
  theme_alex() +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme(legend.title = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits=c(0,32),expand=c(0,0),breaks=c(0,10,20,30), labels = c("0%","10%","20%","30%")) +
  scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue), 
                    labels=c("All","Overlap","4-nt overlap","4-nt overlap, no SD"))

# Pattern
#ggplot(overlaps_4_norbs_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev,
#                                  pattern=seq_overlap_prev))+
#  geom_bar_pattern(width = 1,stat = "identity", size=1) +
#  coord_polar("y", start = 0) +
#  theme_alex() +
#  scale_colour_manual(values = c("black","black","black")) +
#  scale_fill_manual(values = rev(c("grey40","grey60","grey80"))) +
#  scale_linetype_manual(values = c(4,1)) +
#  scale_pattern_manual(values = c ("none","circle","none")) +
#  xlab("") +
#  ylab("") +
#  ggtitle("4nt-noRBS") +
#  scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
#  theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)),
#        axis.ticks = element_blank(), axis.line = element_blank())

###########################################################################

#########
# STATS #
#########

# RRHO on function groups

func_annot <- as.data.frame(func_all) %>% filter(dataset == "all") %>% select(function_group,n)
func_overlap <- as.data.frame(func_all) %>% filter(dataset == "overlap") %>% select(function_group,n)
func_4nt <- func_all %>% filter(dataset == "overlap_4") %>% select(function_group,n)
func_norbs <- func_all %>% filter(dataset == "overlap_4_noRBS") %>% select(function_group,n)

func_4nt <- func_4nt %>% add_row (function_group = "unknown", n=0, .before = 9)
func_norbs <- func_norbs %>% add_row (function_group = "unknown", n=0, .before = 9)



#RRHO_1 <- RRHO(func_annot,func_overlap,BY=F,alternative = 'enrichment',stepsize = 1)
#pval_RRHO_1 <- pvalRRHO(RRHO_1, 50)

#RRHO_2 <- RRHO(func_annot,func_4nt,BY=TRUE,alternative = 'enrichment',stepsize = 1)
#pval_RRHO_2 <- pvalRRHO(RRHO_2, 50)

#RRHO_3 <- RRHO(func_annot,func_norbs,BY=TRUE,alternative = 'enrichment',stepsize = 1)
#pval_RRHO_3 <- pvalRRHO(RRHO_3, 50)



# Indivicual hypergeometric pvals

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

p1_fdr <- p.adjust(p1,method = "BH")
p2_fdr <- p.adjust(p2,method = "BH")
p3_fdr <- p.adjust(p3,method = "BH")

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


# Hypergeo for pie
p_atg <- 1-phyper(289-1,2446,4030-2446,514)
p_gtg <- 1-phyper(204-1,1363,4030-1363,514)
p_ttg <- 1-phyper(22-1,187,4030-187,514)

p_pie_fdr <- p.adjust(c(p_atg,p_gtg,p_ttg),method="BY")

# Hypergeo for rare
rare10_all <- nrow(codon_test %>% filter(is_rare10 == T)) #339
rare10_overlap <- nrow(codon_test_overlap %>% filter(is_rare10 == T)) #102, ***
rare10_4nt <- nrow(codon_test_overlap_4 %>% filter(is_rare10 == T)) #92, ***
rare10_norbs <- nrow(codon_test_overlap_4_noRBS %>% filter(is_rare10 == T)) #54 ***

#R10
r10_all <- nrow(codon_test %>% filter(aa10 == "R")) #673
r10_4nt <- nrow(codon_test_overlap_4 %>% filter(aa10 == "R")) #121 ***

#P10
p10_all <- nrow(codon_test %>% filter(aa10 == "P")) #293
p10_4nt <- nrow(codon_test_overlap_4 %>% filter(aa10 == "P")) #74 ***

#A10
a10_all <- nrow(codon_test %>% filter(aa10 == "A")) # 446
a10_4nt <- nrow(codon_test_overlap_4 %>% filter(aa10 == "A")) #54

#Q10
q10_all <- nrow(codon_test %>% filter(aa10 == "Q")) #174
q10_4nt <- nrow(codon_test_overlap_4 %>% filter(aa10 == "Q")) #51 ***

#G7
g7_all <- nrow(codon_test %>% filter(aa7 == "G")) #445
g7_4nt <- nrow(codon_test_overlap_4 %>% filter(aa7 == "G")) #116 ***

#R7
r7_all <- nrow(codon_test %>% filter(aa7 == "R")) #449
r7_4nt <- nrow(codon_test_overlap_4 %>% filter(aa7 == "R")) #66

#V7
v7_all <- nrow(codon_test %>% filter(aa7 == "V")) #292
v7_4nt <- nrow(codon_test_overlap_4 %>% filter(aa7 == "V")) #54 *

#G8
g8_all <- nrow(codon_test %>% filter(aa8 == "G")) # 390
g8_4nt <- nrow(codon_test_overlap_4 %>% filter(aa8 == "G")) #91 ***

###########################################################################

##########
# WRITE! #
##########

# Overlaps table save into tsv
fwrite(print_overlaps, file = here::here("seq/R/mycopore_redux/Out/overlaps_table.csv"))

# 4-overlaps table save into tsv
fwrite(print_overlaps_4, file = here::here("seq/R/mycopore_redux/Out/overlaps_4_table.csv"))

# nt4a save into tsv
fwrite(print_nt4a, file = here::here("seq/R/mycopore_redux/Out/NTGA_start_nonoverlap.csv"))

# uORF_2main
fwrite(print_uORF_main, file = here::here("seq/R/mycopore_redux/Out/uORF_main_overlap.csv"))

# uORF_selfs
fwrite(print_uORF_overlaps, file = here::here("seq/R/mycopore_redux/Out/uORF_self_overlap.csv"))

# uORF_27
fwrite(print_uORF_test, file = here::here("seq/R/mycopore_redux/Out/uORF_or_maybe_not.csv"))

# p values
fwrite(func_pvals, file = here::here("seq/R/mycopore_redux/Out/p_values.csv"))

# Codon test overlap 4 for new test by Tine (TM)
fwrite(codon_test_overlap_4, file = here::here("NUGA_Out/codontests.csv"))
fwrite(codon_test, file = here::here("NUGA_Out/codons_all.csv"))

