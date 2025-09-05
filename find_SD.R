# PhD Final Stretch: Comparisons between organisms

#   • Attempts to better find SD from gff
#   • Will mostly use NUGA_multispecies code I think 
###########################################################################

########
# Libs #
########

invisible(library(here))

# Contains: themes, get_input, colours hash, nin, packages to load
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))

########
# Vars #
########

# Global folder for gff/fasta/codon table
read_path <- paste0(here(),"/NUGA_Data/species/")
#read_path <- paste0(here(),"/seq/R/Data/")
save_path <- paste0(here(),"/NUGA_Out/species/")
# Species folder path (this is output folder of NUGA_multispecies.R)
species_path <- paste0(here(),"/NUGA_Out/species/")

# For classifying nucleotides as pyrimidine/purine
replacements <- c("A"="R","G"="R","T"="Y","U"="Y","C"="Y" )

# For checking RBS, the strand of purines & lengths
RBS_seq <- "RRRRR"
RBS_range_min = 5
RBS_range_max = 15
search_range = 10 # n of nt at 3' to check for aSD
search_mismatch = 1
SD_mismatch = 1
aSD_consensus = "ACCTCC" #CCUCC as per e coli mutagenesis, could change it to ACCTCC (Huber)

# What counts as mycobacteria
myco_list <- c("Mtb","Msm","Mabs","bovis","leprae","avium","marinum")

########
# Main #
########

# from NUGA_compare
read_full <- function(species, do.GC = T){
  
  if(!dir.exists(file.path(species_path,species))){
    print("you twat you gave a wrong argument!")
    return(NULL)
  }
  
  else{
    full_list <- readRDS(paste0(species_path,species,"/",species,"_full.rds"))
    
    # Calc GC content?
    if(do.GC == T){
      num_g <- str_count(toString(full_list[["files"]][["fasta"]]), "G")
      num_c <- str_count(toString(full_list[["files"]][["fasta"]]), "C")
      
      full <- str_length(toString(full_list[["files"]][["fasta"]]))
      
      gc_percent <- ((num_g + num_c)/full)*100
      full_list[["files"]][["GC"]] <- gc_percent
    }
    
    return(full_list)
  }
  
  
  
}

k <- read_full("Mtb")

gff <- k$files$gff
fasta <- k$files$fasta

# get 16S sequences out - NCBI gff
gff_rrn <- gff %>% 
  filter(type == "rRNA") %>%
  mutate(length = end-start) %>%
  filter(length <2000 & length > 1000) %>%
  group_by(start) %>% # THIS IS NEEDED FOR SUBSEQ TO WORK FOR SOME BIZARRE REASON
  mutate(end_seq = ifelse(strand == "+",
                          toString(subseq(fasta,end-search_range,end)),
                          toString(reverseComplement(subseq(fasta,start,start+search_range))))) %>%
  mutate(test = toString(matchPattern(aSD_consensus,end_seq,max.mismatch = search_mismatch)))


SD_seq <- as.character(reverseComplement(DNAString(get_mode(gff_rrn$test))))

# Annotate as per NUGA_multispecies

# +
gff_plus <- gff %>%
  filter(strand == "+") %>%
  mutate(
    prev_start = lag(start, order_by=start),
    prev_end = lag(end, order_by=start),
    next_start = lead(start, order_by=start),
    next_end = lead(end, order_by=start),
    overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
    overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(start, .by_group = TRUE)

# Grab defined range of nt upstream
gff_plus <- gff_plus %>%
  mutate(seq_upstream = toString(subseq(fasta,start-nchar(RBS_seq)-RBS_range_max+1,start-RBS_range_min)))
gff_plus$seq_upstream_replace <- str_replace_all(gff_plus$seq_upstream,replacements)
#gff_plus <- gff_plus %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Grab defined range of nt upstream of next
gff_plus <- gff_plus %>% mutate(seq_next_upstream = toString(subseq(fasta,next_start-nchar(RBS_seq)-RBS_range_max+1,next_start-RBS_range_min)))
#gff_plus$seq_next_upstream <- str_replace_all(gff_plus$seq_next_upstream,replacements)
#gff_plus <- gff_plus %>% mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))

# Get start codon
gff_plus <- gff_plus %>%
  mutate(start_codon = toString(subseq(fasta,start,start+2)))


# -
gff_minus <- gff %>%
  filter(strand == "-") %>%
  mutate(
    prev_start = lead(start, order_by=start),
    prev_end = lead(end, order_by=start),
    next_start = lag(start, order_by=start),
    next_end = lag(end, order_by=start),
    overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
    overlap_next = ifelse(start <= next_end & start > next_start,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(desc(start))

# Grab defined range of nt upstream
gff_minus <- gff_minus %>% 
  filter(end + 30 < fasta@ranges@width) %>%
  mutate(seq_upstream = toString(reverseComplement(subseq(fasta,end+RBS_range_min,end+nchar(RBS_seq)+RBS_range_max-1))))
gff_minus$seq_upstream_replace <- str_replace_all(gff_minus$seq_upstream,replacements)
#gff_minus <- gff_minus %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))


# Grab defined range of nt upstream of next
gff_minus <- gff_minus %>%
  dplyr::mutate(seq_next_upstream = ifelse(!is.na(next_start),toString(reverseComplement(subseq(fasta,next_end+RBS_range_min,next_end+nchar(RBS_seq)+RBS_range_max-1))),NA))
#gff_minus$seq_next_upstream <- str_replace_all(gff_minus$seq_next_upstream,replacements)
#gff_minus <- gff_minus %>% dplyr::mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))

# Get start codon via mismatched antiSD
gff_minus <- gff_minus %>%
  mutate(start_codon = toString(reverseComplement(subseq(fasta,end-2,end))))

gff_annot <- rbind(gff_plus,gff_minus)

gff_annot <- gff_annot %>%
  mutate(testRBS = toString(matchPattern(SD_seq,seq_upstream,max.mismatch = SD_mismatch, with.indels = T)))


# Via free2bind
aSD_3first <- IRanges::reverse(aSD_consensus)


# get the free2bind for aSD vs upstream seq
gff_annot$deltaG = 0
for(i in 1:nrow(gff_annot)){
  msg <- paste("perl", "/Users/polgpeter7/free2bind/free_align.pl", aSD_3first, gff_annot$seq_upstream[i])
  kk <- system(msg, intern = T)
  # get deltaG
  n <- str_match(kk, "(?<=Delta-G).*")
  n <- n[!is.na(n)]
  nn <- as.numeric(str_match(n, "[-0-9.]+$"))
  # get start of overlap
  l <- str_match(kk, "(?<=seq2 binding).*")
  l <- l[!is.na(l)]
  ll <- as.numeric(str_match(l, "[0-9]+$"))
  # get length of overlap
  j <- str_match(kk, "(?<=Length of bound).*")
  j <- j[!is.na(j)]
  jj <- as.numeric(str_match(j, "[0-9]+$"))
  # calc RBS end as: 0-RBS_range_max-RBS-range_min+start of overlap + length of overlap
  free2bind_dist <- RBS_range_max + RBS_range_min - ll - jj
  gff_annot$deltaG[i] <- nn
  gff_annot$RBS_dist[i] <- free2bind_dist
}

# Check number of purines and number of consecutive purines
gff_annot <- gff_annot %>%
  mutate(n_R = str_count(seq_upstream_replace, "R"))

consec_R <- str_extract_all(gff_annot$seq_upstream_replace, "(R)*") %>%
  map_chr(~.x[which.max(nchar(.x))])

gff_annot$n_consec_R <- nchar(consec_R)
gff_annot <- gff_annot %>%
  mutate(boolRBS = ifelse(testRBS == "", F,T),
         boolNonC = ifelse(type == "CDS", F, T))

# testplot - conseq vs deltaG
ggplot(gff_annot, aes(y=n_consec_R, x=deltaG, colour=boolRBS, shape=boolNonC)) +
  geom_point(size = 3) +
  theme_mycopore() +
  scale_x_reverse() +
  scale_colour_manual(values = c(colours$lightlime2, colours$pink))

# testplot - nR vs deltaG
ggplot(gff_annot, aes(y=n_R, x=deltaG, colour=boolRBS, shape=boolNonC)) +
  geom_point(size = 3) +
  theme_mycopore() +
  scale_x_reverse() +
  scale_colour_manual(values = c(colours$lightlime, colours$purple))

# testplot - dist vs deltaG - not great
ggplot(gff_annot, aes(group=RBS_dist, x=RBS_dist, y=deltaG)) +
  geom_boxplot() +
  theme_mycopore()+
  scale_x_reverse() +
  scale_y_reverse()

#### random controls ####
lastnt <- as.numeric(fasta@ranges@width)
lengthRange <- 1000
gff_random <- tibble("type" = "CDS",
                     "start" = 1000,
                     "end" = 2000,
                     "strand" = "+",
                     "locus_name" = "test") %>%
  ungroup()

for(i in 1:nrow(gff_annot)){
  # randomise length
  geneLength <- floor(runif(1, min = 50, max = 500))
  geneName <- paste("ran",i, sep="")
  # assign strand randomly
  if(runif(1) >= 0.5){
    # + strand
    geneStrand <- "+"
    geneStart <- floor(runif(1, min = lengthRange, max = lastnt - lengthRange))
    geneEnd <- geneStart + floor(runif(1,min=10,max=lengthRange))
  }
  else{
    geneStrand <- "-"
    geneEnd <- floor(runif(1, min = 1, max = lastnt - lengthRange))
    geneStart <- geneEnd - floor(runif(1,min=10,max=lengthRange))
  }
  gff_random <- gff_random %>%
    add_row(type = "CDS",
            start = geneStart,
            end = geneEnd,
            strand = geneStrand,
            locus_name = geneName)
}

gff_random <- gff_random %>%
  filter(locus_name != "test") %>%
  arrange(start) %>%
  filter(!duplicated(start)) %>%
  filter(!duplicated(end))


# +
gff_plus_random <- gff_random %>%
  filter(strand == "+") %>%
  mutate(
    prev_start = lag(start, order_by=start),
    prev_end = lag(end, order_by=start),
    next_start = lead(start, order_by=start),
    next_end = lead(end, order_by=start),
    overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
    overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(start, .by_group = TRUE)

# Grab defined range of nt upstream
gff_plus_random  <- gff_plus_random  %>%
  mutate(seq_upstream = toString(subseq(fasta,start-nchar(RBS_seq)-RBS_range_max+1,start-RBS_range_min)))
gff_plus_random $seq_upstream_replace <- str_replace_all(gff_plus_random $seq_upstream,replacements)
#gff_plus <- gff_plus %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))

# Grab defined range of nt upstream of next
gff_plus_random <- gff_plus_random %>% mutate(seq_next_upstream = toString(subseq(fasta,next_start-nchar(RBS_seq)-RBS_range_max+1,next_start-RBS_range_min)))
#gff_plus_random$seq_next_upstream <- str_replace_all(gff_plus_random$seq_next_upstream,replacements)
#gff_plus_random <- gff_plus_random %>% mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))

# Get start codon
gff_plus_random <- gff_plus_random %>%
  mutate(start_codon = toString(subseq(fasta,start,start+2)))


# -
gff_minus_random <- gff_random %>%
  filter(strand == "-") %>%
  mutate(
    prev_start = lead(start, order_by=start),
    prev_end = lead(end, order_by=start),
    next_start = lag(start, order_by=start),
    next_end = lag(end, order_by=start),
    overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
    overlap_next = ifelse(start <= next_end & start > next_start,TRUE,FALSE )) %>%
  group_by(start) %>%
  arrange(desc(start))

# Grab defined range of nt upstream
gff_minus_random <- gff_minus_random %>% 
  filter(end + 30 < fasta@ranges@width) %>%
  mutate(seq_upstream = toString(reverseComplement(subseq(fasta,end+RBS_range_min,end+nchar(RBS_seq)+RBS_range_max-1))))
gff_minus_random$seq_upstream_replace <- str_replace_all(gff_minus_random$seq_upstream,replacements)
#gff_minus_random <- gff_minus_random %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))


# Grab defined range of nt upstream of next
gff_minus_random <- gff_minus_random %>%
  dplyr::mutate(seq_next_upstream = ifelse(!is.na(next_start),toString(reverseComplement(subseq(fasta,next_end+RBS_range_min,next_end+nchar(RBS_seq)+RBS_range_max-1))),NA))
#gff_minus_random$seq_next_upstream <- str_replace_all(gff_minus_random$seq_next_upstream,replacements)
#gff_minus_random <- gff_minus_random %>% dplyr::mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))

# Get start codon via mismatched antiSD
gff_minus_random <- gff_minus_random %>%
  mutate(start_codon = toString(reverseComplement(subseq(fasta,end-2,end))))

gff_annot_random <- rbind(gff_plus_random,gff_minus_random)

gff_annot_random <- gff_annot_random %>%
  mutate(testRBS = toString(matchPattern(SD_seq,seq_upstream,max.mismatch = SD_mismatch, with.indels = T)))


# get the free2bind for aSD vs upstream seq
gff_annot_random$deltaG = 0
for(i in 1:nrow(gff_annot_random)){
  msg <- paste("perl", "/Users/polgpeter7/free2bind/free_align.pl", aSD_3first, gff_annot_random$seq_upstream[i])
  kk <- system(msg, intern = T)
  n <- str_match(kk, "(?<=Delta-G).*")
  n <- n[!is.na(n)]
  nn <- as.numeric(str_match(n, "[-0-9.]+$"))
  # get start of overlap
  l <- str_match(kk, "(?<=seq2 binding).*")
  l <- l[!is.na(l)]
  ll <- as.numeric(str_match(l, "[0-9]+$"))
  # get length of overlap
  j <- str_match(kk, "(?<=Length of bound).*")
  j <- j[!is.na(j)]
  jj <- as.numeric(str_match(j, "[0-9]+$"))
  # calc RBS end as: 0-RBS_range_max-RBS-range_min+start of overlap + length of overlap
  free2bind_dist <- RBS_range_max + RBS_range_min - ll - jj
  gff_annot_random$deltaG[i] <- nn
  gff_annot_random$RBS_dist[i] <- free2bind_dist
}

# Check number of purines and number of consecutive purines
gff_annot_random <- gff_annot_random %>%
  mutate(n_R = str_count(seq_upstream_replace, "R"))

consec_R_random <- str_extract_all(gff_annot_random$seq_upstream_replace, "(R)*") %>%
  map_chr(~.x[which.max(nchar(.x))])

gff_annot_random$n_consec_R <- nchar(consec_R_random)
gff_annot_random <- gff_annot_random %>%
  mutate(boolRBS = ifelse(testRBS == "", F,T),
         boolNonC = ifelse(type == "CDS", F, T))

# testplot
ggplot(gff_annot_random, aes(y=n_consec_R, x=deltaG, colour=boolRBS, shape=boolNonC)) +
  geom_point(size = 3) +
  theme_mycopore() +
  scale_x_reverse() +
  scale_colour_manual(values = c(colours$lightlime2, colours$pink))

# testplot
ggplot(gff_annot_random, aes(y=n_R, x=deltaG, colour=boolRBS, shape=boolNonC)) +
  geom_point(size = 3) +
  theme_mycopore() +
  scale_x_reverse() +
  scale_colour_manual(values = c(colours$lightlime, colours$purple))

# testplot - dist vs deltaG - not great
ggplot(gff_annot_random, aes(group=RBS_dist, x=RBS_dist, y=deltaG)) +
  geom_boxplot() +
  theme_mycopore()+
  scale_x_reverse() +
  scale_y_reverse()

#### Sum gff/random values and plot ####

sum_deltaG <- gff_annot %>%
  ungroup() %>%
  group_by(deltaG) %>%
  summarise(n = length(deltaG))

sum_deltaG_random <- gff_annot_random %>%
  ungroup() %>%
  group_by(deltaG) %>%
  summarise(n = length(deltaG))

res1 <- gff_annot %>%
  ungroup() %>%
  select(20:25) %>%
  mutate(set = "real")

res2 <- gff_annot_random %>%
  ungroup() %>%
  select(16:21) %>%
  mutate(set = "fake")

res_full <- rbind(res1,res2)

ggplot(res_full, aes(deltaG, n_consec_R)) + 
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = set)) +
  theme_bw() +
  scale_x_reverse()

ggplot() +
  theme_mycopore() +
  geom_ribbon(data=sum_deltaG, aes(x = deltaG, y = n, ymin=0, ymax=predict(loess(n ~ deltaG, span = 0.4))),
              colour=colours$azure, fill = colours$azure, alpha = 0.6) +
  geom_ribbon(data=sum_deltaG_random, aes(x = deltaG, y = n, ymin=0, ymax=predict(loess(n ~ deltaG, span = 0.4))),
              colour=colours$pink, fill = colours$pink, alpha = 0.6) +
  scale_x_reverse(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-1,125))

# ggplot(res_full, aes(x = deltaG, y = set)) +
#   geom_density_ridges2() +
#   theme_mycopore() +
#   scale_x_reverse()

# expand aSD by 1 each side
# multispecies
# distance?
