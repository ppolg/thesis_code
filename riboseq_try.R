
# PhD Final Stretch: testing riboseq coverage per SD binding calcs
#   • NOTE: not included in final thesis, here for posterity's sake.
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

search_upper = -21 # search cov to the start of gene
search_lower = 20 # nt to search upstream of start for cov

riboseq_threshold = 2
riboret_threshold = 2

# Species folder path (this is output folder of NUGA_multispecies.R)
species_path <- paste0(here(),"/NUGA_Out/species/")

########
# Func #
########

read_full <- function(species, do.GC = T){
  
  if(!dir.exists(file.path(species_path,species))){
    print("you gave a wrong argument!")
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


########
# Init #
########

# get riboseq coverage per nt
riboseq_table <- read.csv("/Users/polgpeter7/NUGA_Data/RiboSEQ_Theresa.txt", sep = "\t") %>%
  dplyr::rename(minus = 1, plus = 2) %>%
  rownames_to_column(var = "pos")

Mtb_full <- read_full("mtb")
mtb_table <- Mtb_full[["annot"]][["gff_annot"]] %>%
  filter(type == "CDS")

########
# MAIN #
########
 
# mtb_table <- mtb_table %>%
#   mutate(ribo_cov = ifelse(strand == "+",
#                            sum(riboseq_table$plus[start-search_lower:start-search_upper]),
#                            sum(riboseq_table$minus[3:4])))


#### PLUS ####

# # plus_bg <- mean(riboseq_table$plus)
# 
# #subset table
# mtb_plus <- mtb_table %>%
#   filter(strand == "+") %>%
#   filter(start-search_lower >= 1)
# 
# plusrow <- c()
# for(i in 1:nrow(mtb_plus)){
#   j <- as.numeric(mtb_plus$start[i])
#   l <- (j-search_lower)
#   z <- (j-search_upper)
#   # print(paste0(l," AND ", z))
#   k <- max(riboseq_table$plus[l:j])
#   plusrow[i] <- k
# }
# mtb_plus$riboseq_upstream <- plusrow
# 
# #### MINUS ####
# 
# #subset table
# mtb_minus <- mtb_table %>%
#   filter(strand == "-") %>%
#   filter(end+search_lower < nrow(riboseq_table))
# 
# minusrow <- c()
# for(i in 1:nrow(mtb_minus)){
#   j <- as.numeric(mtb_minus$end[i])
#   l <- (j+search_lower)
#   z <- (j+search_upper)
#   # print(paste0(l," AND ", z))
#   k <- max(riboseq_table$minus[l:j])
#   minusrow[i] <- k
# }
# mtb_minus$riboseq_upstream <- minusrow
# 
# mtb_combined <- rbind(mtb_plus,mtb_minus) %>%
#   arrange(start) %>%
#   mutate(start_codon = ifelse(start_codon %in% c("ATG","TTG","CTG","GTG"),start_codon,"other"))
# 
# ggplot(mtb_combined,aes(x=deltaG,y=riboseq_upstream, colour = start_codon)) +
#   geom_point() +
#   theme_mycopore() +
#   scale_y_log10() +
#   ggtitle("peaks vs deltaG")

# try 6 - with Sawyer et al 2019's counts
mmc2 <- read.csv("/Users/polgpeter7/NUGA_Data/mmc2.csv") %>%
  mutate(across(everything(), ~replace(., . ==  "#DIV/0!" , 0)),
         ribocc1 = as.numeric(ribocc1),
         ribocc2 = as.numeric(ribocc2),
         ribocc3 = as.numeric(ribocc3)) %>%
  rowwise() %>%
  mutate(
    n = mean(c_across(starts_with("n")), na.rm = TRUE),
    rpkm = mean(c_across(starts_with("rpkm")), na.rm = TRUE),
    cds = mean(c_across(starts_with("cds")), na.rm = TRUE),
    translatome = mean(c_across(starts_with("translatome")), na.rm = TRUE),
    trans_rpkm = mean(c_across(starts_with("trans_rpkm")), na.rm = TRUE),
    tr_rpkm_cds = mean(c_across(starts_with("tr_RPKM_CDS")), na.rm = TRUE),
    ribocc = mean(c_across(starts_with("ribocc")), na.rm = TRUE)) %>%
  ungroup()

h <- mtb_combined$deltaG[match(x = mmc2$gene, table=mtb_combined$Rv_name)]
l <- mtb_combined$RBS_dist[match(x = mmc2$gene, table=mtb_combined$Rv_name)]
z <- mtb_combined$start_codon[match(x = mmc2$gene, table=mtb_combined$Rv_name)]

mmc2 <- mmc2 %>%
  add_column(deltaG = h,
             RBS_dist = l,
             start_codon = z) %>%
  filter(!is.na(h)) %>%
  mutate(start_codon = ifelse(start_codon %in% c("ATG","CTG","TTG","GTG"),start_codon,"other")) %>%
  select(1,23:32)

# ggplot(mmc2,aes(x=deltaG,y=rpkm)) +
#   geom_point() +
#   theme_mycopore() +
#   scale_y_log10() +
#   ggtitle("rpkm vs deltaG")
# 
# ggplot(mmc2,aes(x=deltaG,y=n)) +
#   geom_point() +
#   theme_mycopore() +
#   scale_y_log10() +
#   ggtitle("n vs deltaG")
# 
# ggplot(mmc2,aes(x=deltaG,y=cds)) +
#   geom_point() +
#   theme_mycopore() +
#   scale_y_log10() +
#   ggtitle("cds vs deltaG")

ggplot(mmc2,aes(x=deltaG,y=translatome, colour=start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("translatome vs deltaG")

ggplot(mmc2,aes(x=deltaG,y=trans_rpkm, colour=start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("trans_rpkm vs deltaG")

ggplot(mmc2,aes(x=deltaG,y=tr_rpkm_cds, colour=start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("tr_rpkm_CDS vs deltaG")

ggplot(mmc2,aes(x=deltaG,y=ribocc, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("ribocc vs deltaG")

ggplot(mmc2,aes(x=RBS_dist,y=trans_rpkm, group=RBS_dist)) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("trans_rpkm vs RBS distance")

# what if counting peaks over threshold instead?
limited_riboseq <- riboseq_table %>%
  mutate(plus = ifelse(plus >= riboseq_threshold,1,0),
         minus = ifelse(minus >= riboseq_threshold,1,0))
# calc number of nts above limit

#subset table
mtb_plus_limit <- mtb_table %>%
  filter(strand == "+") %>%
  filter(start-search_lower >= 1)

plusrow_limit <- c()
for(i in 1:nrow(mtb_plus_limit)){
  j <- as.numeric(mtb_plus_limit$start[i])
  l <- (j-search_lower)
  z <- (j-search_upper)
  # print(paste0(l," AND ", z))
  k <- sum(limited_riboseq$plus[l:j])
  plusrow_limit[i] <- k
}
mtb_plus_limit$riboseq_upstream <- plusrow_limit

#### MINUS ####

#subset table
mtb_minus_limit <- mtb_table %>%
  filter(strand == "-") %>%
  filter(end+search_lower < nrow(limited_riboseq))

minusrow_limit <- c()
for(i in 1:nrow(mtb_minus_limit)){
  j <- as.numeric(mtb_minus_limit$end[i])
  l <- (j+search_lower)
  z <- (j+search_upper)
  # print(paste0(l," AND ", z))
  k <- sum(limited_riboseq$minus[l:j])
  minusrow_limit[i] <- k
}
mtb_minus$riboseq_upstream <- minusrow_limit

mtb_combined_limit <- rbind(mtb_plus_limit,mtb_minus_limit) %>%
  arrange(start) %>%
  mutate(riboseq_upstream = ifelse(!is.na(riboseq_upstream),riboseq_upstream,0))

ggplot(mtb_combined_limit,aes(x=deltaG,y=riboseq_upstream)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("count of peaks vs deltaG")


# with wade's bigwigs?

riboseq_minus_1 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep1-fp_minusst.bw")$score,
             width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep1-fp_minusst.bw")))

riboseq_plus_1 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep1-fp_plusst.bw")$score,
                      width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep1-fp_plusst.bw")))

riboseq_minus_2 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep2-fp_minusst.bw")$score,
                       width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep2-fp_minusst.bw")))

riboseq_plus_2 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep2-fp_plusst.bw")$score,
                      width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/MC6020-rep2-fp_plusst.bw")))


#riboseq
riboseq_full <- as.tibble(cbind(riboseq_plus_1,riboseq_minus_1,riboseq_plus_2,riboseq_minus_2)) %>%
  dplyr::rename("plus1" = 1,
         "plus2" = 3,
         "minus1" = 2,
         "minus2" = 4)

k <- (riboseq_full$plus1 + riboseq_full$plus2)/2
l <- (riboseq_full$minus1 + riboseq_full$minus2)/2

riboseq_full$plus <- k
riboseq_full$minus <- l

riboseq_full <- riboseq_full %>%
  select(plus,minus) %>%
  mutate(plus = ifelse(plus >= riboseq_threshold,plus,0),
         minus = ifelse(minus >= riboseq_threshold,minus,0))

#subset table
mtb_plus <- mtb_table %>%
  filter(strand == "+") %>%
  filter(start-search_lower >= 1)

plusrow <- c()
for(i in 1:nrow(mtb_plus)){
  j <- as.numeric(mtb_plus$start[i])
  l <- (j-search_lower)
  z <- (j-search_upper)
  # print(paste0(l," AND ", z))
  k <- mean(riboseq_full$plus[l:j])
  plusrow[i] <- k
}
mtb_plus$riboseq_upstream <- plusrow

#### MINUS ####

#subset table
mtb_minus <- mtb_table %>%
  filter(strand == "-") %>%
  filter(end+search_lower < nrow(riboseq_table))

minusrow <- c()
for(i in 1:nrow(mtb_minus)){
  j <- as.numeric(mtb_minus$end[i])
  l <- (j+search_lower)
  z <- (j+search_upper)
  # print(paste0(l," AND ", z))
  k <- mean(riboseq_full$minus[l:j])
  minusrow[i] <- k
}
mtb_minus$riboseq_upstream <- minusrow

mtb_combined <- rbind(mtb_plus,mtb_minus) %>%
  arrange(start) %>%
  mutate(start_codon = ifelse(start_codon %in% c("ATG","TTG","CTG","GTG"),start_codon,"other"))

ggplot(mtb_combined,aes(x=deltaG,y=riboseq_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("peaks vs deltaG")

ggplot(mtb_combined,aes(x=deltaG,y=riboseq_upstream, group = round(deltaG))) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("peaks vs deltaG")

#### riboret ####

riboret_minus_1 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_minus.bw")$score,
                       width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_minus.bw")))
riboret_minus_1_loc = rtracklayer::start(import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_minus.bw"))

riboret_plus_1 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_plus.bw")$score,
                      width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_plus.bw")))
riboret_plus_1_loc = rtracklayer::start(import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_rep1_retap_3_ends_disp_plus.bw"))

riboret_minus_2 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_minus.bw")$score,
                       width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_minus.bw")))
riboret_minus_2_loc = rtracklayer::start(import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_minus.bw"))

riboret_plus_2 <- rep(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_plus.bw")$score,
                      width(rtracklayer::import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_plus.bw")))
riboret_plus_2_loc = rtracklayer::start(import.bw("/Users/polgpeter7/Ribo_data/mc2_7000_all_3_ends_rep2_disp_plus.bw"))

nts <- tibble(pos = 1:4411532)

min1 <- as.tibble(cbind(riboret_minus_1,riboret_minus_1_loc))
min2 <- as.tibble(cbind(riboret_minus_2,riboret_minus_2_loc))
plus1 <- as.tibble(cbind(riboret_plus_1,riboret_plus_1_loc))
plus2 <- as.tibble(cbind(riboret_plus_2,riboret_plus_2_loc))

riboret_table <- nts %>%
  mutate(minus1 = min1$riboret_minus_1[match(x = nts$pos, table=min1$riboret_minus_1_loc)],
         minus2 = min2$riboret_minus_2[match(x = nts$pos, table=min2$riboret_minus_2_loc)],
         plus1 = plus1$riboret_plus_1[match(x = nts$pos, table=plus1$riboret_plus_1_loc)],
         plus2 = plus2$riboret_plus_2[match(x = nts$pos, table=plus2$riboret_plus_2_loc)])

riboret_table[is.na(riboret_table)] <- 0
  

k <- (riboret_table$plus1 + riboret_table$plus2)/2
l <- (riboret_table$minus1 + riboret_table$minus2)/2

riboret_table$plus <- k
riboret_table$minus <- l

riboret_full <- riboret_table %>%
  select(plus,minus) %>%
  mutate(plus = ifelse(plus >= riboret_threshold,plus,0),
         minus = ifelse(minus >= riboret_threshold,minus,0))

#subset table - plusz
mtb_plus2 <- mtb_table %>%
  filter(strand == "+") %>%
  filter(start-search_lower >= 1)

plusrow <- c()
for(i in 1:nrow(mtb_plus2)){
  j <- as.numeric(mtb_plus2$start[i])
  l <- (j-search_lower)
  z <- (j-search_upper)
  # print(paste0(l," AND ", z))
  k <- mean(riboret_full$plus[l:j])
  plusrow[i] <- k
}
mtb_plus2$riboret_upstream <- plusrow

#subset table - minus
mtb_minus2 <- mtb_table %>%
  filter(strand == "-") %>%
  filter(end+search_lower < nrow(riboret_table))

minusrow <- c()
for(i in 1:nrow(mtb_minus2)){
  j <- as.numeric(mtb_minus2$end[i])
  l <- (j+search_lower)
  z <- (j+search_upper)
  # print(paste0(l," AND ", z))
  k <- mean(riboret_full$minus[l:j])
  minusrow[i] <- k
}
mtb_minus2$riboret_upstream <- minusrow

mtb_combined2 <- rbind(mtb_plus2,mtb_minus2) %>%
  arrange(start) %>%
  mutate(start_codon = ifelse(start_codon %in% c("ATG","TTG","CTG","GTG"),start_codon,"other"))

ggplot(mtb_combined2,aes(x=deltaG,y=riboret_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("peaks vs deltaG")

ggplot(mtb_combined2,aes(x=deltaG,y=riboret_upstream, group = round(deltaG))) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("peaks vs deltaG")

#### correlations, testing ####
#20 on either side riboret
# find riboret internally
# calc 
# see logo around?
# read joe paper to see if they done it before though!

mtb_combined$riboret_upstream <- mtb_combined2$riboret_upstream
riboseq_per_start <- aggregate(riboseq_upstream ~  start_codon, mtb_combined, mean)
riboret_per_start <- aggregate(riboret_upstream ~  start_codon, mtb_combined, mean)



ggplot(mtb_combined,aes(x=riboseq_upstream,y=riboret_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("riboret-riboseq corr")

ggplot(mtb_combined,aes(x=start_codon,y=riboseq_upstream, group = start_codon)) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("riboseq per start")

ggplot(mtb_combined,aes(x=start_codon,y=riboret_upstream, group = start_codon)) +
  geom_boxplot() +
#  geom_text(data = riboret_per_start, aes(label = riboret_upstream, y = riboret_upstream + 0.1)) +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("riboret per start") 

combined_full <- mtb_combined %>%
  ungroup() %>%
  mutate(dataset = "full") %>%
  select(locus_name,start_codon,deltaG,riboseq_upstream,riboret_upstream,dataset)

combined_NUGA <- mtb_combined %>%
  ungroup() %>%
  filter(locus_name %in% Mtb_full[["subset"]][["4nt_5"]][["locus_name"]]) %>%
  mutate(dataset = "NUGA") %>%
  select(locus_name,start_codon,deltaG,riboseq_upstream,riboret_upstream,dataset)

combined_URRUG <- mtb_combined %>%
  ungroup() %>%
  filter(locus_name %in% Mtb_full[["subset"]][["1nt_5"]][["locus_name"]])  %>%
  mutate(dataset = "URRUG") %>%
  select(locus_name,start_codon,deltaG,riboseq_upstream,riboret_upstream,dataset)

combined_noRBS <- mtb_combined %>%
  ungroup() %>%
  filter(locus_name %in% Mtb_full[["subset"]][["4nt_norbs_5"]][["locus_name"]])  %>%
  mutate(dataset = "noRBS") %>%
  select(locus_name,start_codon,deltaG,riboseq_upstream,riboret_upstream,dataset)

combined_yesRBS <- mtb_combined %>%
  ungroup() %>%
  filter(locus_name %in% Mtb_full[["subset"]][["4nt_yesrbs_5"]][["locus_name"]]) %>%
  mutate(dataset = "yesRBS") %>%
  select(locus_name,start_codon,deltaG,riboseq_upstream,riboret_upstream,dataset)

combined_pergroup <- rbind(combined_full,combined_NUGA,combined_URRUG,combined_noRBS,combined_yesRBS)

ggplot(combined_pergroup,aes(x=dataset,y=riboseq_upstream, group = dataset)) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
  ggtitle("riboseq per group")

# combined_forplot <- combined_pergroup %>%
#   filter(riboret_upstream >= 0)
ggplot(combined_pergroup,aes(x=dataset,y=riboret_upstream, group = dataset)) +
  geom_boxplot() +
  theme_mycopore() +
  scale_y_log10() +
#  scale_y_continuous(limits = c(0,2))
  ggtitle("riboret per group")

riboseq_per_group <- aggregate(riboseq_upstream ~  dataset, combined_pergroup, mean)
riboret_per_group <- aggregate(riboret_upstream ~  dataset, combined_pergroup, mean)

# corrs per group
ggplot(combined_full,aes(x=riboseq_upstream,y=riboret_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("full riboret-riboseq corr")

ggplot(combined_noRBS,aes(x=riboseq_upstream,y=riboret_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("norbs riboret-riboseq corr")

ggplot(combined_yesRBS,aes(x=riboseq_upstream,y=riboret_upstream, colour = start_codon)) +
  geom_point() +
  theme_mycopore() +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("yesrbs riboret-riboseq corr")

