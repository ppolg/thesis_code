
# PhD Final Stretch: Analysis for multiple organisms
# Updated with Sd finding and others
# Does not use git after all - too much work, not enough time.

#   • Reads fasta/gff/hive table
#   • Performs analysis 
#   • churns out way too many tables and figures to include
###########################################################################


#### Libs ####


invisible(library(here))

# Contains: themes, get_input, colours hash, nin, packages to load
invisible(source(here("seq/R/mycopore_redux/mycopore_init.R")))


#### Vars ####


# Global folder for gff/fasta/codon table
 read_path <- paste0(here(),"/NUGA_Data/species/")
#read_path <- paste0(here(),"/seq/R/Data/")
 save_path <- paste0(here(),"/NUGA_Out/species/")

# For classifying nucleotides as pyrimidine/purine
replacements <- c("A"="R","G"="R","T"="Y","U"="Y","C"="Y" )

# For checking RBS, the strand of purines & lengths
RBS_seq <- "RRRRR" #largely decripated
RBS_range_min = 5 
RBS_range_max = 15
search_range = 10 # n of nt at 3' to check for aSD
search_mismatch = 1 # mismatches for finding SD
SD_mismatch = 1 #mismatches for pairing SD
aSD_consensus = "ACCTCCT" #CCUCC as per e coli mutagenesis, could change it to ACCTCC (Huber)
SD_consensus = "AGGAGGT"
deltaG_limit = -3.5 # as per Huber et al
# complement switching because string changes are weird :/
complements <- c("A"="1","G"="2","T"="3","C"="4",
                 "1"="T","2"="C","3"="A","4"="G" )


# For checking codons:
n_codons = 11
n_codons_past = 11 #for the "last5" part - this checks after NUGA

aa_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("aa",i)
  aa_positions=c(aa_positions,k)
}

aa_positions_past = c()
for(i in 2:n_codons_past){    # -1 to leave stop codon off
  k = paste0("aa",i)
  aa_positions_past=c(aa_positions_past,k)
}

# For the CGA-ending logo
aa_positions_nolast = c()
for(i in 3:n_codons-2){    # -2 to leave stop codon and CGA ending off
  k = paste0("aa",i)
  aa_positions_nolast=c(aa_positions_nolast,k)
}

rare_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("is_rare",i)
  rare_positions=c(rare_positions,k)
}

rare_positions_past = c()
for(i in 2:n_codons){    # -1 to leave stop codon off
  k = paste0("is_rare",i)
  rare_positions_past=c(rare_positions_past,k)
}

acid_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("is_acid",i)
  acid_positions=c(acid_positions,k)
}

base_positions = c()
for(i in 2:n_codons-1){    # -1 to leave stop codon off
  k = paste0("is_base",i)
  base_positions=c(base_positions,k)
}

# "Name" of last codon for plotting later
last_aa <- paste0("aa",n_codons-1)
last_rare <- paste0("is_rare",n_codons-1)
last_acid <- paste0("is_acid",n_codons-1)
last_base <- paste0("is_base",n_codons-1)

last_aa_past <- paste0("aa",n_codons_past)
last_rare_past <- paste0("is_rare",n_codons_past)
last_acid_past <- paste0("is_acid",n_codons_past)
last_base_past <- paste0("is_base",n_codons_past)

# For checking proline doublets
proline_seq <- "PP"

# For minimum threshold for RER
RER_min <- 10

# For checking nt around NUGA 
nt_range <- 30 # with 0 being the N in NUGA


nt_positions = c()
for(i in -nt_range:nt_range){
  k = paste0("nt",i)
  nt_positions=c(nt_positions,k)
}

# COG GO analysis:
# Get full COG tables
COG_full <- read.csv(paste0(here(),"/COG/cog-20.cog.csv"), header = F) %>%
  dplyr::rename("Gene_ID" = 1,
                "NCBI_ID" = 2,
                "Protein_ID" = 3,
                "length" = 4,
                "footprint_coord" = 5,
                "footprint_length" = 6,
                "COG_ID" = 7,
                "reserved" = 8,
                "COG_class" = 9,
                "PSI_BLAST_bit" = 10,
                "PSI_BLAST_evalue" = 11,
                "COG_profile_length" = 12,
                "protein_footprint_coord" = 13)

# Definition and func cat per COG
COG_def <- read.csv(paste0(here(),"/COG/cog-20.def.tab.txt"), header = F, sep = "\t") %>%
  dplyr::rename("COG_ID" = 1,
                "func_cat" = 2,
                "COG_name" = 3,
                "gene_ass" = 4,
                "pathway_ass" = 5,
                "PubMed_ID" = 6,
                "PDB_ID" = 7)

# COG function group definitions
COG_functions <- read.csv(paste0(here(),"/COG/fun-20.tab.txt"), header = F, sep = "\t") %>%
  dplyr::rename("func_letter" = 1,
                "RGB_hex" = 2,
                "desc" = 3)


# Attempt to fix mismatches in check_multi_P:
placeholder <- tibble(triplet = c("CCA","CCC","CCG","CCT"),
                      codon1 = 0,
                      codon2 = 0,
                      codon3 = 0,
                      codon4 = 0,
                      codon5 = 0,
                      codon6 = 0,
                      codon7 = 0,
                      codon8 = 0,
                      codon9 = 0,
                      codon10 = 0,
                      total = 0,
                      ratio = 0)


#### Functions ####

# Function to change raw codon table into readable form
fix_codon_table <- function(table,out_path,from_file = T, main_directory = T, save_file = F,...){
  
  # Read the codon table (that is just copy-paste from DNA HIVE)
  
  if(from_file ==T){
    codon_path <- ifelse(main_directory ==T, paste0(read_path,table,"_codons.tsv",sep=""),table) 
     codon_table <- read.table(codon_path, sep = "\t")
   }
  
  print("Codon table found, fixing it to a usable format")
  # Merge columns
  codon_table_fixed <- rbind(codon_table[1:3], 
                             setNames(codon_table[4:6], names(codon_table)[1:3]),
                             setNames(codon_table[7:9], names(codon_table)[1:3]),
                             setNames(codon_table[10:12], names(codon_table)[1:3]))
  
  # Rename, filter empty rows, remove parentheses in col3
  codon_table_fixed <- codon_table_fixed %>%
    dplyr::rename("triplet" = 1, "per_thousand" = 2, "all" = 3) %>%
    filter(triplet %nin% c("","\t","\n"," "),!is.na(triplet)) %>%
    mutate(all = as.numeric(str_extract(all, "[0-9]+")))
  
  # Add columns
  codon_table_fixed$amino_acid <- as.character(translate(DNAStringSet(codon_table_fixed$triplet)))
  
  # Sum per aa
  codon_sums <- codon_table_fixed %>%
    group_by(amino_acid) %>%
    summarise(sum = sum(all), n = n())
  
  # Add matching sum to codon table, get fraction
  h <- codon_sums$sum[match(codon_table_fixed$amino_acid, table=codon_sums$amino_acid)]
  codon_table_fixed$sum <- h
  
  codon_table_fixed <- codon_table_fixed %>%
    mutate(fraction = round(all/sum, digits = 2)) %>%
    select(triplet,amino_acid,fraction,per_thousand)
  
  print("Returning fixed table")
   return(codon_table_fixed)
  # return(codon_sums)
  
}


# Read gff and optionally filter for type
get_gff_input <- function(gff_file,filter_type="none",...) {
  input_gff <- paste(read_path, gff_file, ".gff", sep = "") 
  gff <- read.gff(input_gff)
  
  # Filter if argument passed, to get only rRNA, tRNA, CDS etc...
  print("Got GFF")
  if(filter_type != "none"){ 
    gff %>% dplyr::filter(type == filter_type)
  }
  else{invisible(gff)}
}


# Function to read files and return as large list
read_files <- function(species="Mtb",keep_gff_isoforms=F, gff_type = 1, processed_codons = T,multiple_chromosomes = F,...) {
  # Read all files needed - they will be in read_path
  # : FASTA
  fasta_path <- paste0(read_path,species,".fasta", sep="")
  fasta <- readDNAStringSet(fasta_path)
  print("Read FASTA")
  # : GFF
  gff <- get_gff_input(species)
  
  if(multiple_chromosomes == T){
    chromosome_add <- gff %>%
      filter(type == "region")
    
        chromosome_add$add = 0
    
    for(row in 2:nrow(chromosome_add)){
      chromosome_add$add[row] = chromosome_add$end[row-1]
    }
    
    chromosome_add$end = cumsum(chromosome_add$add)
    
     chromosome_add <- chromosome_add %>%
       mutate (add = end)
  }
  
  # Annotate the gff, although would be nice if people started actually adhering to standards :(
  if(gff_type == 1){
    
    #Mtb Mycobrowser
    gff <- gff %>%
      dplyr::mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
      mutate(Rv_name = str_split_fixed(str_split_fixed(attributes, ";Name=", 2)[,1], "Locus=", 2)[,2],) %>%
      mutate(product_name = str_split_fixed(str_split_fixed(attributes, ";Comments=", 2)[,1], "Product=", 2)[,2],) %>%
      mutate(function_group = str_split_fixed(str_split_fixed(attributes, ";Protein Data Bank", 2)[,1], "Functional_Category=", 2)[,2],) %>%
      dplyr::select(type,start,end,strand,attributes,locus_name, Rv_name, product_name, function_group)
  }
  
  else if(gff_type == 2){
    # NCBI refseq, no function groups. Old
    gff <- gff %>%
      filter(type == "gene") %>%
      dplyr::mutate(old_locus_name = str_extract(attributes, "(?<=old_locus_tag=).*"))  %>%
      mutate(gene_name = str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "Name=", 2)[,2],) %>%
      mutate(gene_biotype = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene_biotype=", 2)[,2],) %>%
      mutate(type = ifelse(gene_biotype == "protein_coding", "CDS", gene_biotype)) %>%
      mutate(locus_name = ifelse(nchar(gene_name) < 6, gene_name, old_locus_name)) %>%
      dplyr::select(seqid,type,start,end,strand,attributes,locus_name,old_locus_name, gene_name, gene_biotype)
    
  }
  
  else if(gff_type == 3){
    # NCBI refseq, no function groups :(
    gff <- gff %>%
      filter(type == "gene") %>%
      mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "Name=", 2)[,2],) %>%
      mutate(gene_biotype = str_split_fixed(str_split_fixed(attributes, ";gene_synonym", 2)[,1], "gene_biotype=", 2)[,2],) %>%
      mutate(type = ifelse(gene_biotype == "protein_coding", "CDS", gene_biotype)) %>%
      dplyr::select(seqid,type,start,end,strand,attributes,locus_name,gene_biotype)
    
  }
  
  else if(gff_type == 4){
    # NCBI refseq, no function groups, locus_tag
    gff <- gff %>%
      filter(type == "gene") %>%
      mutate(locus_name = str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "Name=", 2)[,2],) %>%
      mutate(gene_biotype = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene_biotype=", 2)[,2],) %>%
      mutate(type = ifelse(gene_biotype == "protein_coding", "CDS", gene_biotype)) %>%
      dplyr::select(seqid,type,start,end,strand,attributes,locus_name,gene_biotype)
    
  }
  
  if(keep_gff_isoforms==F){
    gff <- gff %>% 
      dplyr::distinct(start, .keep_all = T) %>%
      dplyr::distinct(end, .keep_all=T)
  }

  # Fix start/end with chromosome
  if(multiple_chromosomes == T){
    h <- chromosome_add$add[match(x = gff$seqid, table=chromosome_add$seqid)]
    gff$h <- h
    gff <- gff %>%
      mutate(start = start+h,
             end = end + h) %>%
      select(!(h))
  }
  
  print("GFF modified/fixed")
  # : codon table
  if(processed_codons == T){
    codon_path <- paste0(read_path,species,"_codons.csv",sep="")
    codons <- read.csv(codon_path)
  }
  else{
    codons <- fix_codon_table(table = species, from_file = T, main_directory = T)
  }
  
  
  codons <- codons %>%
    dplyr::mutate(is_rare = ifelse(per_thousand < 5,TRUE,FALSE),
                  is_start = ifelse(triplet %in% c("AUG","GUG","UUG"),TRUE,FALSE),
                  is_acid = ifelse(amino_acid %in% c("D","E"),TRUE,FALSE),
                  is_base = ifelse(amino_acid %in% c("R","K","H"),TRUE,FALSE))
  
    # Replace U with T to allow for direct reading from the coding strand DNA
  codons$triplet <- str_replace_all(codons$triplet,"U","T")
  
  print("Codon table read and fixed")
  
  # Species name
  s <- deparse(substitute(species))
  s <- gsub("\"", "", s, fixed = TRUE)
  
  # List up to return
  files_list = list("fasta"=fasta,
                    "gff"=gff,
                    "codons"=codons,
                    "species"= s
                      )
  return(files_list)
}

# Function for annotating each ORF as overlapped
# includes improved SD finding with free2bind deltaG
annotate_overlaps <- function(input_list, strand="both", from_list = T,  manual_gff, manual_fasta, manual_SD = F,...){
  # Append the gff file with overlaps of plus or minus strand
  # Return the gff file
  print("Annotating overlaps")
  # set gff and fasta based on input
  if(from_list == T){
    gff <- input_list$gff
    fasta <- input_list$fasta
  }
  
  else{
    gff <- manual_gff
    fasta <- manual_fasta
  }
  
  #### Plus strand: ####
  if(strand %in% c("+","both")){
    
    # Get plus, starts, overlaps
    gff_plus <- gff %>%
      filter(strand == "+") %>%
      mutate(
        prev_start = lag(start, order_by=start,default = F),
        prev_end = lag(end, order_by=start,default = F),
        next_start = lead(start, order_by=start,default = F),
        next_end = lead(end, order_by=start,default = F),
        overlap_prev = ifelse(start <= prev_end & end > prev_end,TRUE,FALSE),
        overlap_next = ifelse(end >= next_start & end < next_end,TRUE,FALSE )) %>%
      group_by(start) %>%
      arrange(start, .by_group = TRUE)
    
    # Grab defined range of nt upstream
    gff_plus <- gff_plus %>%
      mutate(seq_upstream = toString(subseq(fasta,start-nchar(RBS_seq)-RBS_range_max+1,start-RBS_range_min)))
    #gff_plus$seq_upstream <- str_replace_all(gff_plus$seq_upstream,replacements)
    #gff_plus <- gff_plus %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))
    
    # Grab defined range of nt upstream of next
    gff_plus <- gff_plus %>% mutate(seq_next_upstream = toString(subseq(fasta,next_start-nchar(RBS_seq)-RBS_range_max+1,next_start-RBS_range_min)))
    #gff_plus$seq_next_upstream <- str_replace_all(gff_plus$seq_next_upstream,replacements)
    #gff_plus <- gff_plus %>% mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))
    
    # Get start codon
    gff_plus <- gff_plus %>%
      mutate(start_codon = toString(subseq(fasta,start,start+2)))
    print("+ strand done")
  }
  

  #### Minus strand: ####
  if(strand %in% c("-","both")){
    
    # Get minus, starts, overlaps
    gff_minus <- gff %>%
      filter(strand == "-") %>%
      mutate(
        prev_start = lead(start, order_by=start,default = F),
        prev_end = lead(end, order_by=start,default = F),
        next_start = lag(start, order_by=start,default = F),
        next_end = lag(end, order_by=start,default = F),
        overlap_prev = ifelse(end >= prev_start & end < prev_end,TRUE,FALSE),
        overlap_next = ifelse(start <= next_end & start > next_start,TRUE,FALSE )) %>%
      group_by(start) %>%
      arrange(desc(start))
    
    # Grab defined range of nt upstream
    gff_minus <- gff_minus %>% 
      filter(end + 30 < fasta@ranges@width) %>%
      mutate(seq_upstream = toString(reverseComplement(subseq(fasta,end+RBS_range_min,end+nchar(RBS_seq)+RBS_range_max-1))))
    #gff_minus$seq_upstream <- str_replace_all(gff_minus$seq_upstream,replacements)
    #gff_minus <- gff_minus %>% mutate(has_RBS = ifelse(str_detect(seq_upstream,RBS_seq),TRUE,FALSE))
    
    
    # Grab defined range of nt upstream of next
    gff_minus <- gff_minus %>%
      dplyr::mutate(seq_next_upstream = ifelse(!is.na(next_start),toString(reverseComplement(subseq(fasta,next_end+RBS_range_min,next_end+nchar(RBS_seq)+RBS_range_max-1))),NA))
    #gff_minus$seq_next_upstream <- str_replace_all(gff_minus$seq_next_upstream,replacements)
    #gff_minus <- gff_minus %>%
    #  dplyr::mutate(next_has_RBS = ifelse(str_detect(seq_next_upstream,RBS_seq),TRUE,FALSE))
    
    # Get start codon
    gff_minus <- gff_minus %>%
      mutate(start_codon = toString(reverseComplement(subseq(fasta,end-2,end))))
    print("- strand done")
  }

  # Combine the two:
  if(strand=="both"){
    gff_annot <- rbind(gff_plus,gff_minus)
    print("strands combined")
  }
  
  
  #### 16S ####
  
  print("getting 16S rRNA for SD determination")
  if(manual_SD == F){
    gff_rrn <- gff %>% 
      filter(type == "rRNA") %>%
      mutate(length = end-start) %>%
      filter(length < 2000 & length > 1000) %>%
      group_by(start) %>% # THIS IS NEEDED FOR SUBSEQ TO WORK FOR SOME BIZARRE REASON
      mutate(end_seq = ifelse(strand == "+",
                              toString(subseq(fasta,end-search_range,end)),
                              toString(reverseComplement(subseq(fasta,start,start+search_range))))) %>%
      mutate(test = toString(matchPattern(aSD_consensus,end_seq,max.mismatch = search_mismatch)))
    
    
    SD_seq <- as.character(reverseComplement(DNAString(get_mode(gff_rrn$test))))
  } else {
    SD_seq <- ""
  }
  
  if(SD_seq == ""){
    SD_seq <- SD_consensus
    print("Cannot find SD, will use predefined consensus")
  }
  
  print(paste0("Extended Consesnus SD is ",SD_seq))
  
  aSD_3first <- str_replace_all(SD_seq,complements)
  
  print(paste0("Finding SD based on deltaG to ",aSD_3first))
  
  gff_annot$deltaG = 0
  print(paste0("Long for loop: ",nrow(gff_annot)," iterations..."))
  for(i in 1:nrow(gff_annot)){
    if(i%%500 == 0){
      print(as.character(i))
    }
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
  
  print("SD sequences analysis complete!!")

  # Determine of SD found based on deltaG
  
  gff_annot <- gff_annot %>%
    mutate(has_RBS = ifelse(deltaG < deltaG_limit,T,F))
  
  work_plus <- gff_annot %>%
    filter(strand == "+") %>%
    ungroup() %>%
    mutate(next_has_RBS = lead(has_RBS,order_by = start,default = F)) %>%
    mutate(next_deltaG = lead(deltaG,order_by = start,default = F)) %>%
    mutate(next_RBS_dist = lead(RBS_dist,order_by = start,default = F)) %>%
    group_by(start) %>%
    arrange(start, .by_group = TRUE)
  
  work_minus <- gff_annot %>%
    filter(strand == "-") %>%
    ungroup() %>%
    mutate(next_has_RBS = lag(has_RBS,order_by = start,default = F)) %>%
    mutate(next_deltaG = lag(deltaG,order_by = start,default = F)) %>%
    mutate(next_RBS_dist = lag(RBS_dist,order_by = start,default = F)) %>%
    group_by(start) %>%
    arrange(start, .by_group = TRUE)
    
  gff_annot <- rbind(work_plus,work_minus) %>%
    ungroup() %>%
    group_by(start) %>%
    arrange(start, .by_group = TRUE) 
  
  #### QC: distance vs deltaG regression ####
  
  fig_dist_deltaG <- gff_annot %>%
    ggplot(aes(y=deltaG, x=RBS_dist, group=RBS_dist)) +
    geom_boxplot() +
    theme_mycopore() +
    xlab("distance to RBS (nt)") +
    ylab("deltaG (kcal/mol)") +
    scale_x_reverse() +
    scale_y_reverse()
  
  print(fig_dist_deltaG)
  
  
  s <- c(strand)
  
  # Conditionally export/return - if list, append, else return
  
  if(from_list == T){
    switch(s,
           "+"={
             input_list$gff_annot <- gff_plus
             input_list$fig_dist_deltaG <- fig_dist_deltaG
             return(input_list)
           },
           "-"={
             input_list$gff_annot <- gff_minus
             input_list$fig_dist_deltaG <- fig_dist_deltaG
             return(input_list)
           },
           "both"={
             input_list$gff_annot <- gff_annot
             input_list$fig_dist_deltaG <- fig_dist_deltaG
             return(input_list)
           },
           stop("Please have '+', '-' or 'both' as strand")
    )
  }
  else{
    switch(s,
           "+"={
             return(gff_plus)
           },
           "-"={
             return(gff_minus)
           },
           "both"={
             return(gff_annot)
           },
           stop("Please have '+', '-' or 'both' as strand")
    )
  }
}

# Function to make subsets from annot_overlaps(overlapping,4nt,4nt_norbs, 1nt)
subset_overlaps <- function(input_list, from_list = T, manual_gff_annot, manual_fasta, filter_inframe = T,...){
  print("Subsetting overlaps")
  # set gff and fasta based on input
  if(from_list == T){
    gff <- input_list$gff_annot
    fasta <- input_list$fasta
  }
  
  else{
    gff <- manual_gff_annot
    fasta <- manual_fasta
  }
  
  input_species <- ifelse(from_list == T,input_list$species,deparse(substitute(fasta)))
  
  gff_overlaps <- gff %>%
    filter(overlap_next == TRUE | overlap_prev == TRUE )
  
  #### Plus strand ####
  gff_OL_plus <- gff_overlaps %>%
    filter(strand=="+") %>%
    dplyr::mutate(
      n_overlap_prev = ifelse(overlap_prev == TRUE,prev_end - start +1,NA),
      n_overlap_next = ifelse(overlap_next == TRUE,end - next_start +1,NA),
      seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(subseq(fasta,start,prev_end)), NA),
      seq_overlap_next = ifelse(overlap_next == TRUE, toString(subseq(fasta,next_start,end)), NA)
    )
  
  gff_OL_plus$group = NA
  gff_OL_plus$group[1] = as.numeric(1)
  for (row in 2:nrow(gff_OL_plus)) {
    if (gff_OL_plus$overlap_prev[row] == TRUE) {
      gff_OL_plus$group[row] <- gff_OL_plus$group[row-1]
    }
    else {
      gff_OL_plus$group[row] <- gff_OL_plus$group[row-1] + 1
    }
  }
  
  print("+ strand subset")
  
  #### Minus strand ####
  gff_OL_minus <- gff_overlaps %>%
    filter(strand=="-") %>%
    dplyr::mutate(
      n_overlap_prev = ifelse(overlap_prev == TRUE,end - prev_start +1,NA),
      n_overlap_next = ifelse(overlap_next == TRUE,next_end - start +1,NA),
      seq_overlap_prev = ifelse(overlap_prev == TRUE, toString(reverseComplement(subseq(fasta,prev_start,end))), NA),
      seq_overlap_next = ifelse(overlap_next == TRUE, toString(reverseComplement(subseq(fasta,start,next_end))), NA)
    )
  
  gff_OL_minus$group = NA
  gff_OL_minus$group[1] = as.numeric(gff_OL_plus$group[nrow(gff_OL_plus)])
  for (row in 2:nrow(gff_OL_minus)) {
    if (gff_OL_minus$overlap_prev[row] == TRUE) {
      gff_OL_minus$group[row] <- gff_OL_minus$group[row-1]
    }
    else {
      gff_OL_minus$group[row] <- gff_OL_minus$group[row-1] + 1
    }
  }
  
  print("- strand subset")
  
  #### Combine ####
  overlaps_all <- rbind(gff_OL_plus,gff_OL_minus)
  
  # overlaps_all <- overlaps_all %>%
  #   mutate(frame = ifelse(is.na(n_overlap_next),0,n_overlap_next%%3))
  
  # if(filter_inframe == T){
  #   overlaps_all <- overlaps_all %>%
  #     filter(frame != 0)
  # }
  
  #### Subsets ####
  

  
  #3' overlaps
  overlaps_3 <- overlaps_all %>%
    filter(overlap_next == T) %>%
    mutate(frame = n_overlap_next%%3)
  
  if(filter_inframe == T){
    overlaps_3 <- overlaps_3 %>%
      filter(frame != 0)
  }
  
  #5' overlaps
  overlaps_5 <- overlaps_all %>%
    filter(overlap_prev == T) %>%
    mutate(frame = n_overlap_prev%%3)
  
  if(filter_inframe == T){
    overlaps_5 <- overlaps_5 %>%
      filter(frame != 0)
  }
  
  #25nt for plotting
  overlaps_25nt <- overlaps_3 %>%
    filter(n_overlap_next <26 & !is.na(n_overlap_next))
  
  # overlaps of 4 nt
  overlaps_4nt_5 <- overlaps_all %>%
    dplyr::filter(n_overlap_prev==4) %>%
    mutate(frame = n_overlap_prev%%3)
  
  overlaps_4nt_3 <-overlaps_all %>%
    dplyr::filter(n_overlap_next==4) %>%
    mutate(frame = n_overlap_next%%3)
  
  overlaps_4nt <- rbind(overlaps_4nt_5,overlaps_4nt_3)
  
  overlaps_4nt_norbs_5 <- overlaps_4nt_5 %>%
    filter((has_RBS == FALSE))
  
  overlaps_4nt_norbs_3 <- overlaps_4nt_3 %>%
    filter((next_has_RBS == FALSE))
  
  overlaps_4nt_norbs <- rbind(overlaps_4nt_norbs_5,overlaps_4nt_norbs_3)
  
  overlaps_4nt_yesrbs_5 <- overlaps_4nt_5 %>%
    filter((has_RBS == T))
  
  overlaps_4nt_yesrbs_3 <- overlaps_4nt_3 %>%
    filter((next_has_RBS == T))
  
  overlaps_4nt_yesrbs <- rbind(overlaps_4nt_yesrbs_5,overlaps_4nt_yesrbs_3)
  
  # overlaps of 1nt (URRUG)
  overlaps_1nt_5 <- overlaps_all %>%
    dplyr::filter(n_overlap_prev==1) %>%
    mutate(frame = n_overlap_prev%%3)
  
  overlaps_1nt_3 <-overlaps_all %>%
    dplyr::filter(n_overlap_next==1) %>%
    mutate(frame = n_overlap_next%%3)
  
  
  # "25nt","5","3","4nt_5","4nt_3","4nt_norbs_5","4nt_norbs_3"
  return_list <- list(
    "species" = as.character(input_species),
    "full" = input_list$gff_annot,
    "test" = gff_overlaps,
    "25nt" = overlaps_25nt,
    "overlaps" = overlaps_all,
    "5" = overlaps_5,
    "3" = overlaps_3,
    "4nt" = overlaps_4nt,
    "4nt_5" = overlaps_4nt_5,
    "4nt_3" = overlaps_4nt_3,
    "4nt_norbs" = overlaps_4nt_norbs,
    "4nt_norbs_5" = overlaps_4nt_norbs_5,
    "4nt_norbs_3" = overlaps_4nt_norbs_3,
    "4nt_yesrbs" = overlaps_4nt_yesrbs,
    "4nt_yesrbs_5" = overlaps_4nt_yesrbs_5,
    "4nt_yesrbs_3" = overlaps_4nt_yesrbs_3,
    "1nt_5" = overlaps_1nt_5,
    "1nt_3" = overlaps_1nt_3
  )
  
  print("All subsets complete, returning list")
  # Return list
  return(return_list)
}

# Function to draw graphs from subsets
# includes QC figs for SD strength and distance per type
draw_graphs <- function(subset_list, has_functions = F,...){
  print("creating subset graphs")
  list25 <- subset_list[["25nt"]]
  start_codon_list <- subset_list[["4nt_5"]] %>% 
    mutate(start_codon = ifelse(start_codon %in% c("ATG","GTG","TTG"),start_codon,"other"))
  input_species <- noquote(subset_list$species)
  pics <- c()
  
  #### Pre-work ####
  
  # for pie chart - 4nt
  overlaps_4_count <- subset_list[["4nt_5"]] %>%
    dplyr::select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
    filter(seq_overlap_prev %in% c("ATGA","GTGA","TTGA")) %>%
    ungroup() %>%
    arrange(seq_overlap_prev) %>%
    group_by(seq_overlap_prev) %>%
    summarise(n = length(seq_overlap_prev)) %>%
    arrange(n) %>%
    mutate(prop = n / sum(n) *100) %>%
    mutate(ypos = cumsum(n)- 0.5*n )
  
  # for pie chart - noRBS
  
  if(nrow(subset_list[["4nt_norbs_5"]]) != 0){
    overlaps_4_count_noRBS <- subset_list[["4nt_norbs_5"]] %>%
      dplyr::select(type,start,end,strand,locus_name,seq_overlap_prev) %>%
      ungroup() %>%
      arrange(seq_overlap_prev) %>%
      group_by(seq_overlap_prev) %>%
      summarise(n = length(seq_overlap_prev)) %>%
      arrange(n) %>%
      mutate(prop = n / sum(n) *100) %>%
      mutate(ypos = cumsum(n)- 0.5*n )
  }
  else{
    print("WARNING, noRBS dataset does not work!!!")
    overlaps_4_count_noRBS <- overlaps_4_count
  }
  
  
  # for pie chart - start codons
  
  gff_weirds <- subset_list[["full"]] %>%
    filter(start_codon %nin% c("ATG","GTG","TTG","CTG")) 
  
  start_codon_count <- subset_list[["full"]] %>%
    filter(start_codon %in% c("ATG","GTG","TTG","CTG")) %>%
    dplyr::select(type,start,end,strand,locus_name,start_codon) %>%
    ungroup() %>%
    filter(!is.na(start_codon)) %>%
    arrange(start_codon) %>%
    group_by(start_codon) %>%
    summarise(n = length(start_codon)) %>%
    add_row(start_codon = "other", n=nrow(gff_weirds)) %>%
    mutate(prop = n / sum(n) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
    mutate(start_codon = factor(start_codon, levels=c("ATG","GTG","TTG","other")))%>%
    filter(!is.na(start_codon))
  
  # For function groups
  if(has_functions == T){
    total_func <- subset_list[["full"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "all")
    
    overlaps_func <- subset_list[["3"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "overlap")
    
    overlaps_func_4 <- subset_list[["4nt_3"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "4nt")
    
    overlaps_func_4_noRBS <- subset_list[["4nt_norbs_3"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "noRBS")
    
    overlaps_func_4_yesRBS <- subset_list[["4nt_yesrbs_3"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "yesRBS")
    
    overlaps_func_1 <- subset_list[["1nt_3"]] %>%
      dplyr::filter(function_group != "") %>%
      select(type,start,end,strand,locus_name, function_group) %>%
      ungroup() %>%
      arrange(function_group) %>%
      group_by(function_group) %>%
      summarise(n = length(function_group)) %>%
      mutate(n_norm = n/sum(n), dataset= "1nt")
    
    func_all <- rbind(total_func, overlaps_func, overlaps_func_4, 
                      overlaps_func_4_noRBS, overlaps_func_4_yesRBS, overlaps_func_1) %>%
      mutate(dataset = factor(dataset, levels=c("all","overlap","4nt","noRBS","yesRBS","1nt")))
  }

  print("Pre-work complete!")
  
  #### stats ####
  
  list_ATG <- start_codon_list %>%
    filter(start_codon == "ATG")
  
  list_GTG <- start_codon_list %>%
    filter(start_codon == "GTG")
  
  list_TTG <- start_codon_list %>%
    filter(start_codon == "TTG")
  if(nrow(list_TTG) > 5 &
     nrow(list_GTG) > 5 &
     nrow(list_ATG) > 5){
      
      deltaG_ATG_GTG <- t.test(list_ATG$deltaG,list_GTG$deltaG)
      deltaG_ATG_TTG <- t.test(list_ATG$deltaG,list_TTG$deltaG)
      deltaG_GTG_TTG <- t.test(list_GTG$deltaG,list_TTG$deltaG)
      dist_ATG_GTG <- t.test(list_ATG$RBS_dist,list_GTG$RBS_dist)
      dist_ATG_TTG <- t.test(list_ATG$RBS_dist,list_TTG$RBS_dist)
      dist_GTG_TTG <- t.test(list_GTG$RBS_dist,list_TTG$RBS_dist)
      
      pvals <- list(
        p_dG_A_G = deltaG_ATG_GTG$p.value,
        p_dG_A_T = deltaG_ATG_TTG$p.value,
        p_dG_G_T = deltaG_GTG_TTG$p.value,
        p_dist_A_T = dist_ATG_GTG$p.value,
        p_dist_A_T = dist_ATG_TTG$p.value,
        p_dist_A_T = dist_GTG_TTG$p.value
      )
      
      pics$pvals <- pvals
    }
  print("Stats complete!")
  
  #### DRAW THE FIGS ####
  # Draws all the graphs
  
  # Graph - length  of overlaps
  fig_l_dist <- ggplot(list25, aes(n_overlap_next)) +
    geom_histogram(data=subset(list25,frame==1),binwidth = 1, color = "black", fill = colours$alex_R3, size=1) +
    geom_histogram(data=subset(list25,frame==2),binwidth = 1, color = "black", fill = colours$alex_R1, size=1) +
    geom_histogram(data=subset(list25,frame%nin%c(1,2)),binwidth = 1, color = "black", fill = colours$alex_R2, size=1) +
    theme_mycopore() +
    xlab("Overlap length (nt)") +
    ylab("Count") +
    ggtitle(paste0("Overlaps of ", input_species)) +
    scale_linetype_manual(values = c(4,1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits = c(0,25), expand = c(0,0), breaks = c(1,4,7,10,13,16,19,22,25))
  
  print(fig_l_dist)
  pics$fig_l_dist <- fig_l_dist
  
  # NTGA distribution of overlaps
  fig_NTGA_dist <- ggplot(overlaps_4_count,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
    geom_bar(width = 1,stat = "identity") +
#    geom_text(aes(y = ypos, label = n), color = "black", size=5) +
    coord_polar("y", start = 0) +
    theme_mycopore() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    scale_colour_manual(values = c("grey20","grey20","grey20", "grey20")) +
    scale_fill_manual(values = c(colours$brown,colours$red,colours$orange, colours$yellow)) +
    scale_linetype_manual(values = c(4,1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0("NUGA sequences of ", input_species)) +
    scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
    theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))
  
  print(fig_NTGA_dist)
  pics$fig_NTGA_dist <- fig_NTGA_dist
  
  # NTGA distribution of overlaps - no RBS
  fig_NTGA_dist_noRBS <- ggplot(overlaps_4_count_noRBS,aes(x="",y=n, colour=seq_overlap_prev, fill=seq_overlap_prev))+
    geom_bar(width = 1,stat = "identity") +
#    geom_text(aes(y = ypos, label = n), color = "black", size=5) +
    coord_polar("y", start = 0) +
    theme_mycopore() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    scale_colour_manual(values = c("grey20","grey20","grey20", "grey20")) +
    scale_fill_manual(values = c(colours$brown,colours$red,colours$orange, colours$yellow)) +
    scale_linetype_manual(values = c(4,1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0("NUGA sequences of ", input_species, " - no RBS")) +
    scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
    theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))
  
  print(fig_NTGA_dist_noRBS)
  pics$fig_NTGA_dist_noRBS <- fig_NTGA_dist_noRBS
  
  
  print(fig_NTGA_dist_noRBS)
  pics$fig_NTGA_dist_noRBS <- fig_NTGA_dist_noRBS
  
  # Start codon distribution
  fig_start_codon <- ggplot(start_codon_count,aes(x="",y=n, colour=start_codon, fill=start_codon))+
    geom_bar(width = 1,stat = "identity") +
#    geom_text(aes(y = ypos, label = n), color = "black", size=5) +
    coord_polar("y", start = 0) +
    theme_mycopore() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    scale_colour_manual(values = c("grey20","grey20","grey20","grey20")) +
    scale_fill_manual(values = rev(c(colours$yellow,colours$orange,colours$red,colours$brown))) +
    scale_linetype_manual(values = c(4,1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0("Start codons of ", input_species)) +
    scale_y_discrete(expand = c(0,0,0, 0), breaks = c(0)) +
    theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.8)))
  
  print(fig_start_codon)
  pics$fig_start_codon <- fig_start_codon
  
  # Functional groups
  if(has_functions == T) {
    fig_func_dist <- ggplot(func_all, aes(x=function_group, y=n_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      theme_mycopore() +
      coord_flip() +
      xlab("") +
      ylab("") +
      ggtitle(paste0("ORF functions of ", input_species)) +
      theme(legend.title = element_blank()) +
      scale_x_discrete(limits = rev) +
      scale_y_continuous(breaks=c(10,20,30), labels = c("10%","20%","30%"), expand=c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,colours$alex_B1_edge,colours$grey))
    
    print(fig_func_dist)
    pics$fig_func_dist <- fig_func_dist
  }
  
  # deltaG per group
  fig_deltaG_per_group <- ggplot(start_codon_list, aes(x = start_codon, y = deltaG)) +
    geom_boxplot() +
    theme_mycopore() +
    scale_y_reverse() +
    xlab("") +
    ylab("deltaG (kcal/mol)") +
    ggtitle(paste0("deltaG per start - ", input_species)) 
  
  print(fig_deltaG_per_group)
  pics$fig_deltaG_per_group <- fig_deltaG_per_group
  
  # distance per group
  fig_dist_per_group <- ggplot(start_codon_list, aes(x = start_codon, y = RBS_dist)) +
    geom_boxplot() +
    theme_mycopore() +
    scale_y_reverse() +
    coord_flip() +
    xlab("") +
    ylab("distance from TTS (nt)") +
    ggtitle(paste0("RBS distance per start - ", input_species)) 
  
  print(fig_dist_per_group)
  pics$fig_dist_per_group <- fig_dist_per_group
  
  
  #### FINAL ####
  
  print("Figures complete")
  
  # Return figs
  return(pics)
  
}

# Function for last10 analysis
do_last10 <- function(from_subset = T,draw_fig = T, input_subset, 
                      input_import, input_tasks=c("acid","base","rare","nP","strP","Arg"),...){
  # takes: gff_annot from from subset_overlaps, or a separate gff_annot
  # takes: fasta and codons input from read_files
  # does all defined tasks in input_tasks
  
  print("last 10 codons analysis")
  
 if(from_subset == T){
    gff_annot <- input_subset[["full"]]
    overlaps_all <- input_subset[["overlaps"]]
    overlaps_4_3 <- input_subset[["4nt_3"]]
    overlaps_1_3 <- input_subset[["1nt_3"]]
  }
  else{
    gff_annot <- input_subset
  }
  
  fasta <- input_import$fasta
  codons <- input_import$codons
  
  # Define the  actions to save lines
  defined_tasks <- c("acid","base","rare","nP","strP","Arg")
  copypaste_tasks <- c("")
  
  for(k in 1:length(input_tasks)){
    if(input_tasks[k] %nin% defined_tasks){
      print("Sorry, I have no idea what that task is:")
      print(input_tasks[k])
      print(k)
      return("")
    }
  }
  
  #### Mandatory pre-work ####
  
  print("Pre-work:")
  
  # Filter out outside fasta to prevent weird break from misannotation
  limit <- nchar(toString(subset(fasta)))
  gff_annot <- gff_annot %>%
    filter(start < limit & end < limit)
  
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
      dplyr::mutate(!!k := toString(subseq(fasta,(end-3*(n_codons-i+1)+1),(end-3*(n_codons-i+1)+3))))
  }
  
  for(i in 1:n_codons){    # - strand
    k <- paste0("codon",i)
    codon_test_minus <- codon_test_minus %>%
      dplyr::mutate(!!k := toString(reverseComplement(subseq(fasta,(start+3*(n_codons-i)),(start+3*(n_codons-i)+2)))))
  }
  
  # Merge
  codon_test <- rbind(codon_test_plus,codon_test_minus)
  
  # Get aa
  for(i in 1:n_codons){ 
    k <- paste0("codon",i)
    m <- paste0("aa",i)
    h <- codons$amino_acid[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
    codon_test <- codon_test %>%
      add_column(!!m := h)
  }

  
  # Filter out stuff not ending in stop, because clearly wrong
  codon_test <- codon_test %>%
    filter(aa11=="*")
  
  # Get aa string
  codon_test <- codon_test %>%
    unite(col = "full_aa",!!aa_positions,sep="",remove = FALSE, na.rm=TRUE)
  
  # Filter out stuff not 10 aa, due to N in FASTA
  codon_test <- codon_test %>%
    filter(nchar(full_aa) == 10)
  
  print("Last 10 aa and their nt grabbed")
  

  #### Conditional pre-work ####
  
  # Check acidity?
  if("acid" %in% input_tasks){
    for(i in 1:n_codons){ 
      k <- paste0("codon",i)
      m <- paste0("is_acid",i)
      h <- codons$is_acid[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
      codon_test <- codon_test %>%
        add_column(!!m := h)
    }
    codon_test <- codon_test %>%
      rowwise() %>%
      mutate(n_acid = sum(c_across(starts_with("is_acid"))))
  }
  
  # Check rarity?
  if("rare" %in% input_tasks){
    for(i in 1:n_codons){ 
      k <- paste0("codon",i)
      m <- paste0("is_rare",i)
      h <- codons$is_rare[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
      codon_test <- codon_test %>%
        add_column(!!m := h)
    }
    
    codon_test <- codon_test %>%
      rowwise() %>%
      mutate(n_rare = sum(c_across(starts_with("is_rare")))-1) # -1 for stop codon
  }
  
  # Check basicity?
  if("base" %in% input_tasks){
    for(i in 1:n_codons){ 
      k <- paste0("codon",i)
      m <- paste0("is_base",i)
      h <- codons$is_base[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
      codon_test <- codon_test %>%
        add_column(!!m := h)
    }
    
    codon_test <- codon_test %>%
      rowwise() %>%
      mutate(n_base = sum(c_across(starts_with("is_base"))))
  }
  
  # Check proline count?
  if("nP" %in% input_tasks){
    codon_test <- codon_test %>%
      rowwise() %>%
      mutate(n_proline = sum(c_across(!!aa_positions) == "P"))
  }
  
  # Check aa string for P stretch?
  if("strP" %in% input_tasks){
    codon_test <- codon_test %>%
      dplyr::mutate(has_P_stretch = ifelse(str_detect(full_aa,proline_seq),TRUE,FALSE))
  }
  
  print("Acid, base, rare, proline grabbed")
  
  #### Main work ####
  
  # Define return list
  return_list <- c(species = input_subset$species)
  return_list$codon_test <- codon_test
  
  ##### Make subsets: #####
  # Select overlap only
  codon_test_overlap <- codon_test %>%
    dplyr::filter(locus_name %in% overlaps_all$locus_name & overlap_next == TRUE)
  
  # Select 4-overlap only
  codon_test_overlap_4 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_4_3$locus_name)
  
  # Select 4-overlap, no-RBS only
  codon_test_overlap_4_noRBS <- codon_test_overlap_4 %>%
    dplyr::filter(next_has_RBS == FALSE)
  
  # Select 4-overlap, has RBS only
  codon_test_overlap_4_yesRBS <- codon_test_overlap_4 %>%
    dplyr::filter(next_has_RBS == TRUE)
  
  # Select 4-overlap, ending in CGA only
  # drop codon10
  codon_test_overlap_4_CGA <- codon_test_overlap_4 %>%
    dplyr::filter(codon10 == "CGA")
  
  if(nrow(codon_test_overlap_4_CGA) == 0){
    print("no CGA found. Will copy NUGA data to avoid breaks.")
    print("please note this down, it is not indicated downstream!!")
    codon_test_overlap_4_CGA <- codon_test_overlap_4[1]
  }
  
  codon_test_overlap_4_CGA <- select(codon_test_overlap_4_CGA, -any_of(c("is_rare10","is_acid10","is_base10","full_aa"))) %>%
    mutate(is_rare10 = F, is_acid10 = F, is_base10 = F) %>%
    unite(col = "full_aa",!!aa_positions_nolast,sep="",remove = FALSE, na.rm=TRUE)
  
  # select 1nt overlaps only
  codon_test_overlap_1 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_1_3$locus_name)
  
  print("Found data subsetted")
  
  ###### RARE ######
  if("rare" %in% input_tasks){
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
    
    rare_codons_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      dplyr::filter(n_rare >= 0) %>%
      ungroup() %>%
      arrange(n_rare) %>%
      group_by(n_rare) %>%
      summarise(count = length(n_rare)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_yesRBS), dataset = "overlap_4_yesRBS")
    
    rare_codons_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      dplyr::filter(n_rare >= 0) %>%
      ungroup() %>%
      arrange(n_rare) %>%
      group_by(n_rare) %>%
      summarise(count = length(n_rare)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
    
    rare_codons_overlap_1 <- codon_test_overlap_1 %>%
      dplyr::filter(n_rare >= 0) %>%
      ungroup() %>%
      arrange(n_rare) %>%
      group_by(n_rare) %>%
      summarise(count = length(n_rare)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
    
    rare_codons_sum = rbind(rare_codons_all,rare_codons_overlap,rare_codons_overlap_4,
                            rare_codons_overlap_4_noRBS, rare_codons_overlap_4_yesRBS,
                            rare_codons_overlap_4_CGA, rare_codons_overlap_1)
    
    print("rare codons counted")
    
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
    
    rare_test_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      select(starts_with("is_rare")) %>%
      gather(position,rare) %>%
      group_by(position,rare) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4_yesRBS))) %>%
      select(rare,!!rare_positions) %>%
      filter(rare == TRUE) %>%
      mutate(category = "overlap_4_yesRBS")
    
    rare_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      select(starts_with("is_rare")) %>%
      gather(position,rare) %>%
      group_by(position,rare) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
      select(rare,!!rare_positions) %>%
      filter(rare == TRUE) %>%
      mutate(category = "overlap_4_CGA")
    
    rare_test_overlap_1 <- codon_test_overlap_1 %>%
      select(starts_with("is_rare")) %>%
      gather(position,rare) %>%
      group_by(position,rare) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_1))) %>%
      select(rare,!!rare_positions) %>%
      filter(rare == TRUE) %>%
      mutate(category = "overlap_1")
    
    rare_test_all <- rbind(rare_test,rare_test_overlap,rare_test_overlap_4,
                           rare_test_overlap_4_noRBS, rare_test_overlap_4_yesRBS,
                           rare_test_overlap_4_CGA,rare_test_overlap_1) %>%
      select(category, !!rare_positions)
    
    # Add to return
    return_list$rare <- rare_codons_sum
    return_list$rare_all <- rare_test_all
  }
  
  print("rare codons per position counted")
  
  ###### ACID ######
  if("acid" %in% input_tasks){
    # Check acid codon counts for plotting
    acid_codons_all <- codon_test %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test), dataset = "all")
    
    acid_codons_overlap <- codon_test_overlap %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")
    
    acid_codons_overlap_4 <- codon_test_overlap_4 %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")
    
    acid_codons_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")
    
    acid_codons_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_yesRBS), dataset = "overlap_4_yesRBS")
    
    acid_codons_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
    
    acid_codons_overlap_1 <- codon_test_overlap_1 %>%
      dplyr::filter(n_acid >= 0) %>%
      ungroup() %>%
      arrange(n_acid) %>%
      group_by(n_acid) %>%
      summarise(count = length(n_acid)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
    
    acid_codons_sum = rbind(acid_codons_all,acid_codons_overlap,acid_codons_overlap_4,
                            acid_codons_overlap_4_noRBS, acid_codons_overlap_4_yesRBS,
                            acid_codons_overlap_4_CGA,acid_codons_overlap_1)
    
    print("acid codons counted")
    
    # For acid codon per position
    # Get acid composition for all
    acid_test <- codon_test %>%
      select(starts_with("is_acid")) %>%
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
    
    acid_test_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      select(starts_with("is_acid")) %>%
      gather(position,acid) %>%
      group_by(position,acid) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap_4_yesRBS))) %>%
      select(acid,!!acid_positions) %>%
      filter(acid == TRUE) %>%
      mutate(category = "overlap_4_yesRBS")
    
    acid_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      select(starts_with("is_acid")) %>%
      gather(position,acid) %>%
      group_by(position,acid) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
      select(acid,!!acid_positions) %>%
      filter(acid == TRUE) %>%
      mutate(category = "overlap_4_CGA")
    
    acid_test_overlap_1 <- codon_test_overlap_1 %>%
      select(starts_with("is_acid")) %>%
      gather(position,acid) %>%
      group_by(position,acid) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-acid),funs(. / nrow(codon_test_overlap_1))) %>%
      select(acid,!!acid_positions) %>%
      filter(acid == TRUE) %>%
      mutate(category = "overlap_1")
    
    acid_test_all <- rbind(acid_test,acid_test_overlap,acid_test_overlap_4,
                           acid_test_overlap_4_noRBS, acid_test_overlap_4_yesRBS,
                           acid_test_overlap_4_CGA,acid_test_overlap_1) %>%
      select(category, !!acid_positions)
    
    print("acid codons per position counted")
    
    # Add to return
    return_list$acid <- acid_codons_sum
    return_list$acid_all <- acid_test_all
  }
  
  ###### BASE ######
  if("base" %in% input_tasks){
    # Check base codon counts for plotting
    base_codons_all <- codon_test %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test), dataset = "all")
    
    base_codons_overlap <- codon_test_overlap %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")
    
    base_codons_overlap_4 <- codon_test_overlap_4 %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")
    
    base_codons_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")
    
    base_codons_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_yesRBS), dataset = "overlap_4_yesRBS")
    
    base_codons_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
    
    base_codons_overlap_1 <- codon_test_overlap_1 %>%
      dplyr::filter(n_base >= 0) %>%
      ungroup() %>%
      arrange(n_base) %>%
      group_by(n_base) %>%
      summarise(count = length(n_base)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
    
    base_codons_sum = rbind(base_codons_all,base_codons_overlap,base_codons_overlap_4,
                            base_codons_overlap_4_noRBS, base_codons_overlap_4_yesRBS,
                            base_codons_overlap_4_CGA,base_codons_overlap_1)
    
    print("base codons counted")
    
    # For base codon per position
    # Get base composition for all
    base_test <- codon_test %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "all")
    
    base_test_overlap <- codon_test_overlap %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap")
    
    base_test_overlap_4 <- codon_test_overlap_4 %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap_4))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap_4")
    
    base_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap_4_noRBS")
    
    base_test_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap_4_yesRBS))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap_4_yesRBS")
    
    base_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap_4_CGA")
    
    base_test_overlap_1 <- codon_test_overlap_1 %>%
      select(starts_with("is_base")) %>%
      gather(position,base) %>%
      group_by(position,base) %>%
      summarise(no = n()) %>%
      spread(position,no) %>%
      mutate_at(vars(-base),funs(. / nrow(codon_test_overlap_1))) %>%
      select(base,!!base_positions) %>%
      filter(base == TRUE) %>%
      mutate(category = "overlap_1")
    
    #DEBUG!!!
    base_test_all <- rbind(base_test,base_test_overlap,base_test_overlap_4,
                           base_test_overlap_4_noRBS, base_test_overlap_4_yesRBS,
                           base_test_overlap_4_CGA, base_test_overlap_1) %>%
      select(category, !!base_positions)
    
    # Add to return
    return_list$base <- base_codons_sum
    return_list$base_all <- base_test_all
  }
  
  print("base codons per position counted")
  
  ###### PROLINE ######
  
  if("nP" %in% input_tasks){
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
    
    proline_codons_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      ungroup() %>%
      arrange(n_proline) %>%
      group_by(n_proline) %>%
      summarise(count = length(n_proline)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_yesRBS), dataset = "overlap_4_yesRBS")
    
    proline_codons_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      ungroup() %>%
      arrange(n_proline) %>%
      group_by(n_proline) %>%
      summarise(count = length(n_proline)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_CGA")
    
    proline_codons_overlap_1 <- codon_test_overlap_1 %>%
      ungroup() %>%
      arrange(n_proline) %>%
      group_by(n_proline) %>%
      summarise(count = length(n_proline)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
    
    proline_codons_sum = rbind(proline_codons_all,proline_codons_overlap,proline_codons_overlap_4,
                               proline_codons_overlap_4_noRBS, proline_codons_overlap_4_yesRBS,
                               proline_codons_overlap_4_CGA, proline_codons_overlap_1)
    
    print("prolines counted")
    
    return_list$proline <- proline_codons_sum
  }
  ##### Other tasks #####
  
  # Everything else
  
  ###### Proline doublets ######
  
  if("strP" %in% input_tasks){
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
      mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")
    
    proline_stretch_overlap_4 <- codon_test_overlap_4 %>%
      ungroup() %>%
      arrange(has_P_stretch) %>%
      group_by(has_P_stretch) %>%
      summarise(count = length(has_P_stretch)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")
    
    proline_stretch_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
      ungroup() %>%
      arrange(has_P_stretch) %>%
      group_by(has_P_stretch) %>%
      summarise(count = length(has_P_stretch)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")
    
    proline_stretch_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
      ungroup() %>%
      arrange(has_P_stretch) %>%
      group_by(has_P_stretch) %>%
      summarise(count = length(has_P_stretch)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_yesRBS), dataset = "overlap_4_yesRBS")
    
    proline_stretch_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
      ungroup() %>%
      arrange(has_P_stretch) %>%
      group_by(has_P_stretch) %>%
      summarise(count = length(has_P_stretch)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
    
    proline_stretch_overlap_1 <- codon_test_overlap_1 %>%
      ungroup() %>%
      arrange(has_P_stretch) %>%
      group_by(has_P_stretch) %>%
      summarise(count = length(has_P_stretch)) %>%
      mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
    
    proline_stretch_sum <- rbind(proline_stretch_all,
                                 proline_stretch_overlap,
                                 proline_stretch_overlap_4,
                                 proline_stretch_overlap_4_noRBS,
                                 proline_stretch_overlap_4_yesRBS,
                                 proline_stretch_overlap_4_CGA,
                                 proline_stretch_overlap_1) %>%
      dplyr::filter(has_P_stretch == T)
    
    return_list$proline_stretch <- proline_stretch_sum
    
    print("proline doublets counted")
  }
  
  ###### AA analysis ######
  # Get aa composition at all
  aa_test_all <- codon_test %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "all")
  
  # Get aa composition at overlap
  aa_test_overlap <- codon_test_overlap %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "overlap")
  
  # Get aa composition at 4nt-overlap
  aa_test_overlap_4 <- codon_test_overlap_4 %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "overlap_4")
  
  # Get aa composition at 4nt-overlap, no RBS
  aa_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "overlap_4_noRBS")
  
  # Get aa composition at 4nt-overlap, yes RBS
  aa_test_overlap_4_yesRBS <- codon_test_overlap_4_yesRBS %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_yesRBS))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "overlap_4_yesRBS")
  
  # Get aa composition at 4nt-overlap, last CGA
  aa_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
    select(amino,!!aa_positions_nolast) %>%
    mutate(dataset = "overlap_4_CGA")
  
  # Get aa composition at 1nt-overlap
  aa_test_overlap_1 <- codon_test_overlap_1 %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_1))) %>%
    select(amino,!!aa_positions) %>%
    mutate(dataset = "overlap_1")
  
  aa_test_sum <- rbind( aa_test_all,
                        aa_test_overlap,
                        aa_test_overlap_4,
                        aa_test_overlap_4_noRBS,
                        aa_test_overlap_4_yesRBS,
                        aa_test_overlap_1)

  return_list$aa_test <- aa_test_sum
  
  print("aa analysis complete")
  
  ###### ARGININE ######
  
  if("Arg" %in% input_tasks){
    
    #Get list of all Arg codons
    
    arg_codons <- codons %>%
      filter(amino_acid == "R")
    
    arg_codons <- arg_codons$triplet
    
    # Run Arg test - all final
    arg_test <- codon_test %>%
      filter(aa10 == "R")
    
    arg_sum <- arg_test %>%
      group_by(codon10) %>%
      summarise(n = length(codon10))
    
    for(codon in 1:length(arg_codons)){
      if(arg_codons[codon] %nin% arg_sum$codon10){
        arg_sum <- arg_sum %>%
          add_row(codon10 = arg_codons[codon], n=0)
      }
    }
    
    arg_sum <- arg_sum %>% 
      mutate(k = n/nrow(arg_test))
    
    h <- codons$fraction[match(x = arg_sum$codon10, table=codons$triplet)]
    arg_sum_all <- arg_sum %>%
      add_column(h = h) %>%
      select(codon10,k,h) %>%
      dplyr::rename("Last codon" = k,
                    "All codons" = h)
    
    # Run Arg test - NUGA
    arg_test_4 <- codon_test_overlap_4 %>%
      filter(aa10 == "R")
    
    arg_sum_4 <- arg_test_4 %>%
      group_by(codon10) %>%
      summarise(n = length(codon10))
    
    for(codon in 1:length(arg_codons)){
      if(arg_codons[codon] %nin% arg_sum_4$codon10){
        arg_sum_4 <- arg_sum_4 %>%
          add_row(codon10 = arg_codons[codon], n=0)
      }
    }
    
    arg_sum_4 <- arg_sum_4 %>% 
      mutate(k = n/nrow(arg_test_4))
    
    arg_sum_all <- arg_sum_all %>%
      add_column("NUGA" = arg_sum_4$k)
    
    # Run Arg test - no RBS
    arg_test_noRBS <- codon_test_overlap_4_noRBS %>%
      filter(aa10 == "R")
    
    arg_sum_noRBS <- arg_test_noRBS %>%
      group_by(codon10) %>%
      summarise(n = length(codon10))
    
    for(codon in 1:length(arg_codons)){
      if(arg_codons[codon] %nin% arg_sum_noRBS$codon10){
        arg_sum_noRBS <- arg_sum_noRBS %>%
          add_row(codon10 = arg_codons[codon], n=0)
      }
    }
    
    arg_sum_noRBS <- arg_sum_noRBS %>% 
      mutate(k = n/nrow(arg_test_noRBS))
    
    
    arg_sum_all <- arg_sum_all %>%
      add_column("noRBS" = arg_sum_noRBS$k)
  
    return_list$arg_sum <- arg_sum_all
    
  }
  
  # Run Arg test - yes RBS
  arg_test_yesRBS <- codon_test_overlap_4_yesRBS %>%
    filter(aa10 == "R")
  
  arg_sum_yesRBS <- arg_test_yesRBS %>%
    group_by(codon10) %>%
    summarise(n = length(codon10))
  
  for(codon in 1:length(arg_codons)){
    if(arg_codons[codon] %nin% arg_sum_yesRBS$codon10){
      arg_sum_yesRBS <- arg_sum_yesRBS %>%
        add_row(codon10 = arg_codons[codon], n=0)
    }
  }
  
  arg_sum_yesRBS <- arg_sum_yesRBS %>% 
    mutate(k = n/nrow(arg_test_yesRBS))
  
  
  arg_sum_all <- arg_sum_all %>%
    add_column("yesRBS" = arg_sum_yesRBS$k)
  
  # Run Arg test - URRUG
  arg_test_1 <- codon_test_overlap_1 %>%
    filter(aa10 == "R")
  
  arg_sum_1 <- arg_test_1 %>%
    group_by(codon10) %>%
    summarise(n = length(codon10))
  
  for(codon in 1:length(arg_codons)){
    if(arg_codons[codon] %nin% arg_sum_1$codon10){
      arg_sum_1 <- arg_sum_1 %>%
        add_row(codon10 = arg_codons[codon], n=0)
    }
  }
  
  arg_sum_1 <- arg_sum_1 %>% 
    mutate(k = n/nrow(arg_test_1))
  
  arg_sum_all <- arg_sum_all %>%
    add_column("URRUG" = arg_sum_1$k)
  
  return_list$arg_sum <- arg_sum_all
  
  print("argninines counted")
  
  ##### Plotting #####
  
  #Break function if no need to draw figs
  if(draw_fig == F){
    return(return_list)
  }
  print("Plotting...")
  
  #Define list of figs?
  pics <- c()
  
  ###### FIG - RARE ######
  
  if("rare" %in% input_tasks){
    
    fig_rare_sum <- ggplot(rare_codons_sum, aes(x=n_rare, y=count_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      xlab("Number of rare codons") +
      ylab("Relative abundance") +
      ggtitle(paste0("Rare codons in last 10 - ",input_subset$species)) +
      theme_mycopore() +
      theme(legend.title = element_blank()) +
      scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,20,40), labels = c("0%","20%","40%"), expand=c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_rare_sum)
    pics$fig_rare_sum <- fig_rare_sum
    
    fig_rare_test <- rare_test_all %>%
      gather(pos,value,is_rare1:!!last_rare) %>%
      ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
      geom_bar(stat="identity", position = "dodge") +
      theme_mycopore() +
      xlab("Codon position") +
      ylab("Percentage of rare codons") +
      scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,5,10,15,20,25,30), labels = c("0%","5%","10%","15%","20%","25%","30%"), expand=c(0,0)) +
      ggtitle(paste0("Rare codons in last 10 - ",input_subset$species)) +
      theme(legend.title = element_blank()) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_rare_test)
    pics$fig_rare_test <- fig_rare_test
    
  }
 
  ###### FIG - ACID ######
  
  if("acid" %in% input_tasks){
    
    fig_acid_sum <- ggplot(acid_codons_sum, aes(x=n_acid, y=count_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      xlab("Number of acidic codons") +
      ylab("Relative abundance") +
      ggtitle(paste0("Acidic codons in last 10 - ",input_subset$species)) +
      theme_mycopore() +
      theme(legend.title = element_blank()) +
      scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,20,40), labels = c("0%","20%","40%"), expand=c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,
                                   colours$blue,colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_acid_sum)
    pics$fig_acid_sum <- fig_acid_sum
    
    fig_acid_test <- acid_test_all %>%
      gather(pos,value,is_acid1:!!last_acid) %>%
      ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
      geom_bar(stat="identity", position = "dodge") +
      theme_mycopore() +
      xlab("Codon position") +
      ylab("Percentage of acidic codons") +
      scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,5,10,15,20,25,30), labels = c("0%","5%","10%","15%","20%","25%","30%"), expand=c(0,0)) +
      ggtitle(paste0("Acidic codons in last 10 - ",input_subset$species)) +
      theme(legend.title = element_blank()) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_acid_test)
    pics$fig_acid_test <- fig_acid_test
    
  }
  
  ###### FIG - BASE ######
  
  if("base" %in% input_tasks){
    
    fig_base_sum <- ggplot(base_codons_sum, aes(x=n_base, y=count_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      xlab("Number of basic codons") +
      ylab("Relative abundance") +
      ggtitle(paste0("Basic codons in last 10 - ",input_subset$species)) +
      theme_mycopore() +
      theme(legend.title = element_blank()) +
      scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,20,40), labels = c("0%","20%","40%"), expand=c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_base_sum)
    pics$fig_base_sum <- fig_base_sum
    
    fig_base_test <- base_test_all %>%
      gather(pos,value,is_base1:!!last_base) %>%
      ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
      geom_bar(stat="identity", position = "dodge") +
      theme_mycopore() +
      xlab("Codon position") +
      ylab("Percentage of basic codons") +
      scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
      scale_y_continuous(breaks = c(0,5,10,15,20,25,30), labels = c("0%","5%","10%","15%","20%","25%","30%"), expand=c(0,0)) +
      ggtitle(paste0("Basic codons in last 10 - ",input_subset$species)) +
      theme(legend.title = element_blank()) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_base_test)
    pics$fig_base_test <- fig_base_test
    
  }
  
  ###### FIG - PROLINE ######
  if("nP" %in% input_tasks){ 
    fig_proline_n <- ggplot(proline_codons_sum, aes(x=n_proline, y=count_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      theme_mycopore() +
      xlab("Number of prolines") +
      ylab("Relative abundance") +
      ggtitle(paste0("Prolines in last 10 - ",input_subset$species)) +
      theme(legend.title = element_blank()) +
      scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9)) +
      scale_y_continuous(breaks = c(0,10,20,30,40,50), labels = c("0%","10%","20%","30%","40%","50%"), expand = c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1,colours$grey,"black"))
    
    print(fig_proline_n)
    pics$fig_proline_n <- fig_proline_n
  }
  
  if("strP" %in% input_tasks){ 

    fig_proline_str <- ggplot(proline_stretch_sum, aes(x=has_P_stretch, y=count_norm*100, fill=dataset)) +
      geom_bar(position="dodge",stat="identity") +
      coord_flip() +
      theme_mycopore() +
      xlab("") +
      ylab("") +
      ggtitle(paste0("Multi-prolines in last 10 - ",input_subset$species)) +
      scale_y_continuous(breaks=c(1,2,3,4,5), labels = c("1%","2%","3%","4%","5%")) +
      theme(legend.title = element_blank(),
            axis.text.y = element_blank()) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,
                                   colours$alex_B1_edge,colours$grey,"black"))
    
    print(fig_proline_str)
    pics$fig_proline_str <- fig_proline_str
  }
  
  ###### FIG - ARGININE ######
  if("Arg" %in% input_tasks){
    arg_melt <- arg_sum_all %>%
      reshape2::melt()
    
    return_list$arg_melt <- arg_melt
    
    fig_arg_codons <- ggplot(arg_melt, aes(x=codon10, y=value*100, fill=variable)) +
      geom_bar(position="dodge",stat="identity") +
      theme_mycopore() +
      xlab("") +
      ylab("Relative abundance (%)") +
      ggtitle(paste0("Arginine codons - ",input_subset$species)) +
      theme(legend.title = element_blank()) +
      scale_y_continuous(breaks = c(0,10,20,30,40,50,60), labels = c("0%","10%","20%","30%","40%","50%","60%"), expand = c(0,0)) +
      scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,
                                   colours$blue,colours$grey,"black"))
    
    print(fig_arg_codons)
    pics$fig_arg_codons <- fig_arg_codons
  }
  
  ###### FIG - AA logo? ######
  
  # All
  fig_logo_all <- ggplot() +
    geom_logo(codon_test$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of all ORFs - ",input_subset$species)) 
  
  print(fig_logo_all)
  pics$fig_logo_all <- fig_logo_all
  
  # Overlap
  fig_logo_overlap <- ggplot() +
    geom_logo(codon_test_overlap$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap)
  pics$fig_logo_overlap <- fig_logo_overlap
  
  # NUGA
  fig_logo_overlap_4 <- ggplot() +
    geom_logo(codon_test_overlap_4$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of NUGA overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4)
  pics$fig_logo_overlap_4 <- fig_logo_overlap_4
  
  # noRBS
  fig_logo_overlap_4_noRBS <- ggplot() +
    geom_logo(codon_test_overlap_4_noRBS$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of noRBS overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4_noRBS)
  pics$fig_logo_overlap_4_noRBS <- fig_logo_overlap_4_noRBS
  
  # yesRBS
  fig_logo_overlap_4_yesRBS <- ggplot() +
    geom_logo(codon_test_overlap_4_yesRBS$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of YES RBS overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4_yesRBS)
  pics$fig_logo_overlap_4_yesRBS <- fig_logo_overlap_4_yesRBS
  
  # CGA
  fig_logo_overlap_4_CGA <- ggplot() +
    geom_logo(codon_test_overlap_4_CGA$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of NUGA ending in CGA - ",input_subset$species)) 
  
  print(fig_logo_overlap_4_CGA)
  pics$fig_logo_overlap_4_CGA <- fig_logo_overlap_4_CGA
  
  # URRUG
  fig_logo_overlap_1 <- ggplot() +
    geom_logo(codon_test_overlap_1$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("aa of URRUG overlaps - ",input_subset$species))
  
  print(fig_logo_overlap_1)
  pics$fig_logo_overlap_1 <- fig_logo_overlap_1
  
  ##### Finalise, return #####
  return_list$figures <- pics
  
  print("last 10 analysis complete!")
  
  return(return_list)
}


# Test if aa after the overlap show patterns
do_first5 <- function(input_subset,input_import){
  # define from subset
  print("performing first aa analysis")
  
  gff_annot <- input_subset[["full"]]
  overlaps_all <- input_subset[["5"]]
  overlaps_4_5 <- input_subset[["4nt_5"]]
  overlaps_1_5 <- input_subset[["1nt_5"]]
  fasta <- input_import$fasta
  codons <- input_import$codons
  
  pics <- list()
  
  # limit for preventing breaks
  limit <- nchar(toString(subset(fasta)))
  gff_annot <- gff_annot %>%
    filter(start < limit & end < limit)
  
  # Grab + and -
  codon_test_plus <- gff_annot %>%
    dplyr::filter(strand == "+" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  codon_test_minus <- gff_annot %>%
    dplyr::filter(strand == "-" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  #### PRE ####
  # Get first n codons of all
  
  for(i in 1:n_codons_past){    # + strand
    k <- paste0("codon",i)
    codon_test_plus <- codon_test_plus %>%
      dplyr::mutate(!!k := toString(subseq(fasta,(start+3*(i-1)),(2+start+3*(i-1)))))
  }
  
  codon_test_plus <- codon_test_plus %>%
    filter(start > 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(subseq(fasta,start-2,start)) == "CGA",T,F))

  for(i in 1:n_codons_past){    # - strand
    k <- paste0("codon",i)
    codon_test_minus <- codon_test_minus %>%
      dplyr::mutate(!!k := toString(reverseComplement(subseq(fasta,(end-3*(i-1)-2),(end-3*(i-1))))))
  }
  
  codon_test_minus <- codon_test_minus %>%
    filter(end < limit - 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(reverseComplement(subseq(fasta,end,end+2))) == "CGA",T,F))
  
  # Merge
  codon_test <- rbind(codon_test_plus,codon_test_minus) %>%
    ungroup()
  
  # Get aa
  for(i in 1:n_codons_past){ 
    k <- paste0("codon",i)
    m <- paste0("aa",i)
    h <- codons$amino_acid[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
    codon_test <- codon_test %>%
      add_column(!!m := h)
  }
  
  print("aa found")
  
  # Filter out stuff not starting in NTG, because clearly wrong
  # codon_test <- codon_test %>%
  #   filter(aa1 %in% c("ATG","TTG","CTG","TTG"))
  
  # Get aa string
  codon_test <- codon_test %>%
    unite(col = "full_aa",!!aa_positions_past,sep="",remove = FALSE, na.rm=TRUE)
  
  # Filter out stuff not 10 aa, due to N in FASTA
  codon_test <- codon_test %>%
  filter(nchar(full_aa) == n_codons_past-1)
  
  # find P stretches
   codon_test <- codon_test %>%
     dplyr::mutate(has_P_stretch = ifelse(str_detect(full_aa,proline_seq),TRUE,FALSE))
   
  # find rare 
  for(i in 1:n_codons){ 
    k <- paste0("codon",i)
    m <- paste0("is_rare",i)
    h <- codons$is_rare[match(x = codon_test[[rlang::as_name(k)]], table=codons$triplet)]
    codon_test <- codon_test %>%
      add_column(!!m := h)
  }
   
  codon_test <- codon_test %>%
    rowwise() %>%
    mutate(n_rare = sum(c_across(starts_with("is_rare")))) # foo
  
  print("rare, P, string complete")
  
  #### SUBSET ####
  # Select overlap only
  codon_test_overlap <- codon_test %>%
    dplyr::filter(locus_name %in% overlaps_all$locus_name & overlap_prev == TRUE)
   
  # Select 4-overlap only
  codon_test_overlap_4 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_4_5$locus_name)
   
  # Select 4-overlap, no-RBS only
  codon_test_overlap_4_noRBS <- codon_test_overlap_4 %>%
    dplyr::filter(has_RBS == FALSE)
   
  # Select 4-overlap, ending in CGA only
  codon_test_overlap_4_CGA <- codon_test_overlap_4 %>%
    dplyr::filter(test_CGA == TRUE)
  
  # Select 1-overlap only
  codon_test_overlap_1 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_1_5$locus_name)

  print("subsetting done")
  
  #### AA analysis ####
  # Get aa composition at all
  aa_test_all <- codon_test %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "all")
    

  # Get aa composition at overlap
  aa_test_overlap <- codon_test_overlap %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "overlap")
  
  # Get aa composition at 4nt-overlap
  aa_test_overlap_4 <- codon_test_overlap_4 %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "overlap_4")
  
  # Get aa composition at 4nt-overlap, no RBS
  aa_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "overlap_4_noRBS")
  
  # Get aa composition at 4nt-overlap, past CGA
  aa_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "overlap_4_CGA")
  
  aa_test_overlap_1 <- codon_test_overlap_1 %>%
    select(starts_with("aa")) %>%
    gather(position,amino) %>%
    group_by(position,amino) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-amino),funs(. / nrow(codon_test_overlap_1))) %>%
    select(amino,!!aa_positions_past) %>%
    mutate(dataset = "overlap_1")
  
  aa_test_sum <- rbind( aa_test_all,
                        aa_test_overlap,
                        aa_test_overlap_4,
                        aa_test_overlap_4_noRBS,
                        aa_test_overlap_4_CGA,
                        aa_test_overlap_1)
  
  "aa analysis complete"
  
  #### AA logo ####
  
  print("Creating figures...")
  
  # All
  fig_logo_all <- ggplot() +
    geom_logo(codon_test$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of all ORFs - ",input_subset$species)) 
  
  print(fig_logo_all)
  pics$fig_logo_all <- fig_logo_all
  
  # Overlap
  fig_logo_overlap <- ggplot() +
    geom_logo(codon_test_overlap$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap)
  pics$fig_logo_overlap <- fig_logo_overlap
  
  # NUGA
  fig_logo_overlap_4 <- ggplot() +
    geom_logo(codon_test_overlap_4$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of NUGA overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4)
  pics$fig_logo_overlap_4 <- fig_logo_overlap_4
  
  # noRBS
  fig_logo_overlap_4_noRBS <- ggplot() +
    geom_logo(codon_test_overlap_4_noRBS$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of noRBS overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4_noRBS)
  pics$fig_logo_overlap_4_noRBS <- fig_logo_overlap_4_noRBS
  
  # CGA
  fig_logo_overlap_4_CGA <- ggplot() +
    geom_logo(codon_test_overlap_4_CGA$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of CGA overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_4_CGA)
  pics$fig_logo_overlap_4_CGA <- fig_logo_overlap_4_CGA
  
  # URRUG
  fig_logo_overlap_1 <- ggplot() +
    geom_logo(codon_test_overlap_1$full_aa) +
    theme_mycopore() +
    theme(legend.title = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank()) +
    xlab("Amino acid position") +
    ylab("Bits") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)) +
    ggtitle(paste0("first aa of URRUG overlaps - ",input_subset$species)) 
  
  print(fig_logo_overlap_1)
  pics$fig_logo_overlap_1 <- fig_logo_overlap_1
  
  print("Done!")
  
  #### PP stretch ####
  
  print("PP analysis and figures...")
  
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
    mutate(count_norm = count/nrow(codon_test_overlap), dataset = "overlap")
  
  proline_stretch_overlap_4 <- codon_test_overlap_4 %>%
    ungroup() %>%
    arrange(has_P_stretch) %>%
    group_by(has_P_stretch) %>%
    summarise(count = length(has_P_stretch)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_4), dataset = "overlap_4")
  
  proline_stretch_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
    ungroup() %>%
    arrange(has_P_stretch) %>%
    group_by(has_P_stretch) %>%
    summarise(count = length(has_P_stretch)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_4_noRBS), dataset = "overlap_4_noRBS")
  
  proline_stretch_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    ungroup() %>%
    arrange(has_P_stretch) %>%
    group_by(has_P_stretch) %>%
    summarise(count = length(has_P_stretch)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
  
  proline_stretch_overlap_1 <- codon_test_overlap_1 %>%
    ungroup() %>%
    arrange(has_P_stretch) %>%
    group_by(has_P_stretch) %>%
    summarise(count = length(has_P_stretch)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
  
  proline_stretch_sum <- rbind(proline_stretch_all,
                               proline_stretch_overlap,
                               proline_stretch_overlap_4,
                               proline_stretch_overlap_4_noRBS,
                               proline_stretch_overlap_4_CGA,
                               proline_stretch_overlap_1) %>%
    dplyr::filter(has_P_stretch == T)
  
  fig_proline_str <- ggplot(proline_stretch_sum, aes(x=has_P_stretch, y=count_norm*100, fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    coord_flip() +
    theme_mycopore() +
    xlab("") +
    ylab("") +
    ggtitle(paste0("Multi-prolines in first 10 - ",input_subset$species)) +
    scale_y_continuous(breaks=c(1,2,3,4,5), labels = c("1%","2%","3%","4%","5%")) +
    theme(legend.title = element_blank(),
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,colours$alex_B1_edge,colours$grey))
  
  print(fig_proline_str)
  pics$fig_proline_str <- fig_proline_str

  print("Done!")
  #### RARE? ####
  
  print("rare analysis and figures...")
  
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
  
  rare_codons_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    dplyr::filter(n_rare >= 0) %>%
    ungroup() %>%
    arrange(n_rare) %>%
    group_by(n_rare) %>%
    summarise(count = length(n_rare)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_4_CGA), dataset = "overlap_4_CGA")
  
  rare_codons_overlap_1 <- codon_test_overlap_1 %>%
    dplyr::filter(n_rare >= 0) %>%
    ungroup() %>%
    arrange(n_rare) %>%
    group_by(n_rare) %>%
    summarise(count = length(n_rare)) %>%
    mutate(count_norm = count/nrow(codon_test_overlap_1), dataset = "overlap_1")
  
  rare_codons_sum = rbind(rare_codons_all,rare_codons_overlap,rare_codons_overlap_4,rare_codons_overlap_4_noRBS,rare_codons_overlap_4_CGA,rare_codons_overlap_1)
  
  # For rare codon per position
  # Get rare composition for all
  rare_test <- codon_test %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "all")
  
  rare_test_overlap <- codon_test_overlap %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "overlap")
  
  rare_test_overlap_4 <- codon_test_overlap_4 %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "overlap_4")
  
  rare_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "overlap_4_noRBS")
  
  rare_test_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "overlap_4_CGA")
  
  rare_test_overlap_1 <- codon_test_overlap_1 %>%
    select(starts_with("is_rare")) %>%
    gather(position,rare) %>%
    group_by(position,rare) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-rare),funs(. / nrow(codon_test_overlap_1))) %>%
    select(rare,!!rare_positions_past) %>%
    filter(rare == TRUE) %>%
    mutate(category = "overlap_1")
  
  rare_test_all <- rbind(rare_test,rare_test_overlap,rare_test_overlap_4,rare_test_overlap_4_noRBS,rare_test_overlap_4_CGA,rare_test_overlap_1) %>%
    select(category, !!rare_positions_past)

  # figs
  
  fig_rare_sum <- ggplot(rare_codons_sum, aes(x=n_rare, y=count_norm*100, fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    xlab("Number of rare codons") +
    ylab("Relative abundance") +
    ggtitle(paste0("Rare codons in first 10 - ",input_subset$species)) +
    theme_mycopore() +
    theme(legend.title = element_blank()) +
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
    scale_y_continuous(breaks = c(0,20,40), labels = c("0%","20%","40%"), expand=c(0,0)) +
    scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,colours$alex_B1_edge,colours$grey))
  
  print(fig_rare_sum)
  pics$fig_rare_sum <- fig_rare_sum
  
  fig_rare_test <- rare_test_all %>%
    gather(pos,value,is_rare2:!!last_rare_past) %>%
    ggplot(aes(x = fct_inorder(pos), y=value*100, fill=category)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_mycopore() +
    xlab("Codon position") +
    ylab("Percentage of rare codons") +
    scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10)) +
    scale_y_continuous(breaks = c(0,5,10,15,20,25,30), labels = c("0%","5%","10%","15%","20%","25%","30%"), expand=c(0,0)) +
    ggtitle(paste0("Rare codons in first 10 - ",input_subset$species)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = c(colours$orange,colours$red,colours$purple,colours$blue,colours$alex_B1_edge,colours$grey))
  
  print(fig_rare_test)
  pics$fig_rare_test <- fig_rare_test
  
  print("Done!")
  
  
  #### RETURN ####
  
  # Define return list
   return_list <- c(species = input_subset$species)
   return_list$codon_test <- codon_test
   return_list$codon_test_overlap <- codon_test_overlap
   return_list$codon_test_overlap_4 <- codon_test_overlap_4
   return_list$codon_test_overlap_4_noRBS <- codon_test_overlap_4_noRBS
   return_list$codon_test_overlap_4_CGA <- codon_test_overlap_4_CGA
   return_list$codon_test_overlap_1 <- codon_test_overlap_1
   return_list$aa_test <- aa_test_sum
   return_list$proline_stretch <- proline_stretch_sum
   return_list$pics <- pics
   return_list$rare <- rare_codons_sum
   return_list$rare_all <- rare_test_all
   # return_list$codon_test_overlap_4_CGA <- codon_test_overlap_4_CGA
   
   print("Analysis done, returning list")
   
  return(return_list)
}


# relative enrichment ratio (RER)
# "is a certain last codon enriched in 3' NUGA overlapped ORFs?"

calc_NUGA_RER <- function(input_last10,input_subset){
  print("calculating normalised  relative enrichment ratio (RER)")
  # Things to take into account in normalising:
  #   • log2 is good to account for lower ones
  #   • similarly, filter out very very rare codons?
  #   • no CUGA (no CTG start), so only calc frequency of last codons that aren't ending in C
  codon_test <- input_last10$codon_test 
  codon_test_NUGA <- codon_test %>%
    filter(locus_name %in% input_subset$"4nt_3"[["locus_name"]])
  codon_test_URRUG <- codon_test %>%
    filter(locus_name %in% input_subset$"1nt_3"[["locus_name"]])
  
  # nonC for RER of NUGA - N is never C (no CTG start)
  nonC <- codon_test %>%
    filter(substring(codon10,3,3) != "C")
  
  nonC_summed <- nonC %>%
    group_by(codon10) %>%
    summarise(n = length(codon10)) %>%
    mutate(k = (n/nrow(nonC))*1000)
  
  
  all_summed <- codon_test %>%
    group_by(codon10) %>%
    summarise(n2 = length(codon10)) %>%
    mutate(k2 = (n2/nrow(codon_test)) * 1000)
  
  
  # NUGA sum
  
  NUGA_summed <- codon_test_NUGA %>%
    group_by(codon10) %>%
    summarise(n = length(codon10)) %>%
    mutate(k = (n/nrow(codon_test_NUGA))*1000)
  
  # all_table <- nonC_summed %>%
  #   dplyr::rename("C_freq" = "k") %>%
  #   mutate(NUGA_freq = ifelse(codon10 %in% NUGA_summed$codon10,
  #                             "a",
  #                             0))
  
  h <-  NUGA_summed$k[match(x=nonC_summed$codon10,table = NUGA_summed$codon10)]
  h2 <-  NUGA_summed$n[match(x=nonC_summed$codon10,table = NUGA_summed$codon10)]
  
  print("NUGA done!")
  
  # URRUG sum
  URRUG_summed <- codon_test_URRUG%>%
    group_by(codon10) %>%
    summarise(n = length(codon10)) %>%
    mutate(k = (n/nrow(codon_test_URRUG))*1000)
  
  l <-  URRUG_summed$k[match(x=all_summed$codon10,table = URRUG_summed$codon10)]
  l2 <-  URRUG_summed$n[match(x=all_summed$codon10,table = URRUG_summed$codon10)]
  
  
  all_summed_NUGA <- nonC_summed %>%
    add_column(h = h) %>%
    add_column(h2 = h2) %>%
    mutate(h = ifelse(is.na(h),0,h)) %>%
    mutate(h2 = ifelse(is.na(h2),0,h2)) %>%
    dplyr::rename("triplet"= codon10,
                  "n_C" = n,
                  "C_freq" = k,
                  "NUGA_freq" = h,
                  "n_NUGA" = h2) %>%
    mutate(RER = NUGA_freq/C_freq,
           RERlog = log2(RER))
  
  
  all_summed_URRUG <- all_summed  %>%
    add_column(l = l) %>%
    add_column(l2 = l2) %>%
    mutate(n2 = ifelse(is.na(n2),0,n2)) %>%
    mutate(k2 = ifelse(is.na(k2),0,k2)) %>%
    mutate(l = ifelse(is.na(l),0,l)) %>%
    mutate(l2 = ifelse(is.na(l2),0,l2)) %>%
    dplyr::rename("triplet" = codon10,
                  "n_R" = n2,
                  "R_freq" = k2,
                  "URRUG_freq" = l,
                  "n_URRUG" = l2) %>%
    mutate(RER_URRUG = URRUG_freq/R_freq,
           RERlog_URRUG = log2(RER_URRUG))
  
  print("URRUG done!")
  
  # Hypergeo?
  print("performing stats analysis...")
  
  sum_C <- sum(all_summed_NUGA$n_C)
  sum_R <- sum(all_summed_URRUG$n_R)
  sum_NUGA <- sum(all_summed_NUGA$n_NUGA)
  sum_URRUG <- sum(all_summed_URRUG$n_URRUG)
  
  all_summed_NUGA <- all_summed_NUGA %>%
    mutate(pvalue = 1-(phyper(n_NUGA,n_C,sum_C - sum_NUGA,sum_NUGA))) %>%
    mutate(padj = p.adjust(pvalue, method="BH"))
  
  all_summed_URRUG <- all_summed_URRUG %>%
    mutate(pvalue = 1-(phyper(n_URRUG,n_R,sum_R - sum_URRUG,sum_URRUG))) %>%
    mutate(padj = p.adjust(pvalue, method="BH"))
  
  print("done!")
  
  # Plot - occurrences per codon per thousand
  print("plotting...")
  
  fig_occ_NUGA <- all_summed_NUGA %>%
    select(triplet,C_freq,NUGA_freq) %>%
    reshape2::melt() %>%
    ggplot(aes(y=triplet, x=value, fill=variable)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("codons per thousand") +
    ylab("") +
    ggtitle(paste0("Last codons in ",input_subset$species)) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_fill_manual(values = c(colours$red,colours$blue)) +
    theme(panel.grid.major.y = element_blank())
  
  print(fig_occ_NUGA)
  
  fig_occ_URRUG <- all_summed_URRUG %>%
    select(triplet,R_freq,URRUG_freq) %>%
    reshape2::melt() %>%
    ggplot(aes(y=triplet, x=value, fill=variable)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("codons per thousand") +
    ylab("") +
    ggtitle(paste0("Last codons in ",input_subset$species)) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_fill_manual(values = c(colours$red,colours$blue)) +
    theme(panel.grid.major.y = element_blank())
  
  print(fig_occ_URRUG)
  
  # Plot - RER
  
  fig_RER <- ggplot(all_summed_NUGA, aes(y=triplet, x=RER)) +
    geom_bar(position="dodge",stat="identity",fill = colours$azure) +
    theme_mycopore() +
    xlab("Enrichment in NUGA") +
    ylab("") +
    ggtitle(paste0("Enrichment ratios in ",input_subset$species)) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand=c(0,0)) + 
    theme(panel.grid.major.y = element_blank())
  
  print(fig_RER)
  
  fig_RER_URRUG <- all_summed_URRUG %>%
    filter(n_R != "0") %>%
    ggplot(aes(y=triplet, x=RER_URRUG)) +
    geom_bar(position="dodge",stat="identity",fill = colours$azure) +
    theme_mycopore() +
    xlab("Enrichment in URRUG") +
    ylab("") +
    ggtitle(paste0("Enrichment ratios in ",input_subset$species)) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand=c(0,0)) + 
    theme(panel.grid.major.y = element_blank())
  
  print(fig_RER_URRUG)
  
  # Plot - volcano??
  
  fig_volcano <- ggplot(all_summed_NUGA) +
    geom_point(aes(x=RERlog,
                   y=-log10(padj),
                   color=ifelse(padj<0.05,"grey60",colours$alex_R1),
                   fill=ifelse(padj<0.05,"grey60",colours$alex_R1))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("Enrichment volcano plot -  ",input_subset$species))
  
  print(fig_volcano)
  
  fig_volcano_URRUG <- ggplot(all_summed_URRUG) +
    geom_point(aes(x=RERlog_URRUG,
                   y=-log10(padj),
                   color=ifelse(padj<0.05,"grey60",colours$alex_R1),
                   fill=ifelse(padj<0.05,"grey60",colours$alex_R1))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("URRUG volcano plot -  ",input_subset$species))
  
  print(fig_volcano_URRUG)
  
  print("done!")
  
  print("RER analysis done, returning list")
  
  return(list("NUGA" = NUGA_summed, "C" = nonC_summed, "R" = all_summed, "URRUG" = URRUG_summed, 
              "final_NUGA" = all_summed_NUGA, "final_URRUG" = all_summed_URRUG,
              "fig_occ_NUGA" = fig_occ_NUGA, "fig_occ_URRUG" = fig_occ_URRUG,
              "fig_RER" = fig_RER, "fig_volcano" = fig_volcano,
              "fig_RER_URRUG" = fig_RER_URRUG, "fig_volcano_URRUG" = fig_volcano_URRUG))
}

# Alan's stop codon dist function
# this creates a table with just codon 11
show_stop <- function(codon_table){
  print("checking stop codons...")
  codon_test <- codon_table$codon_test %>%
    group_by(codon11)%>%
    summarise(n=length(codon11))
  
  #making the graphs
  
  fig_bar <- ggplot(codon_test, aes(x=codon11, y=n, fill=codon11)) +
    theme_mycopore()+
    geom_bar(stat = "identity")+
    xlab("Stop codons")+
    ylab("Number of codons")+
    ggtitle(paste0("Distribution of stop codons in ",codon_table$species)) +
    # scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values=c(colours$green, colours$azure, colours$alex_R3)) +
    theme(legend.position = "none")
  
  print(fig_bar)
  
  fig_pie <- ggplot(codon_test, aes(x="", y=n, fill=codon11))+
    
    geom_bar(stat="identity", width=1)+
    coord_polar("y", start=0)+
    theme_mycopore() +
    ggtitle(paste0("Distribution of stop codons in ",codon_table$species))+
    ylab("Number of codons")+
    xlab("") +
    scale_fill_manual(values=c(colours$green, colours$azure, colours$alex_R3)) +
    labs(fill = "Stop codon") +
    scale_y_continuous(labels = c("","1000","2000","3000","4000"), breaks=c(0,1000,2000,3000,4000)) +
    theme(panel.grid.major = element_blank())
  
  
  print(fig_pie)
  
  print("done!")
  
  #This is just for the output, so that we can export multiple things
  return_list <- list("table" = codon_test, "figure_bar" = fig_bar, "figure_pie" = fig_pie)
  return(return_list)
  
}


# Test QC of SD finding on subsets
do_SD_QC <- function(input_subset,input_last10){
  
  input_species <- noquote(input_subset$species)
  NUGA_list <- input_subset[["4nt_5"]]
  norbs_list <- input_subset[["4nt_norbs_5"]]
  yesrbs_list <- input_subset[["4nt_yesrbs_5"]]
  URRUG_list <- input_subset[["1nt_5"]]
  
  CGA_test <- input_last10[["codon_test"]] %>%
    filter(codon10 == "CGA")
  
  CGA_list <- input_subset[["4nt_3"]] %>%
    filter(locus_name %in% CGA_test$locus_name)
  
  final_list <- list()
  
  # for getting the deltaG of noRBS and yesRBS
  
  deltaG_4 <- tibble(deltaG = NUGA_list[["deltaG"]]) %>%
    mutate(dataset = "NUGA")
  
  deltaG_4_noRBS <- tibble(deltaG = norbs_list[["deltaG"]]) %>%
    mutate(dataset = "noRBS")
  
  deltaG_4_yesRBS <- tibble(deltaG = yesrbs_list[["deltaG"]]) %>%
    mutate(dataset = "yesRBS")
  
  deltaG_CGA <- tibble(deltaG = CGA_list[["next_deltaG"]]) %>%
    mutate(dataset = "CGA")
  
  deltaG_URRUG <- tibble(deltaG = URRUG_list[["deltaG"]]) %>%
    mutate(dataset = "URRUG")
  
  deltaG_table <- rbind(deltaG_4,deltaG_4_noRBS,deltaG_4_yesRBS,deltaG_CGA,deltaG_URRUG)
  
  # for getting the average distance of noRBS and yesRBS
  dist_4 <- tibble(dist = NUGA_list[["RBS_dist"]]) %>%
    mutate(dataset = "NUGA")
  
  dist_4_noRBS <- tibble(dist = norbs_list[["RBS_dist"]]) %>%
    mutate(dataset = "noRBS")
  
  dist_4_yesRBS <- tibble(dist = yesrbs_list[["RBS_dist"]]) %>%
    mutate(dataset = "yesRBS")
  
  dist_CGA <- tibble(dist = CGA_list[["next_RBS_dist"]]) %>%
    mutate(dataset = "CGA")
  
  dist_URRUG <- tibble(dist = URRUG_list[["RBS_dist"]]) %>%
    mutate(dataset = "URRUG")
  
  dist_table <- rbind(dist_4,dist_4_noRBS,dist_4_yesRBS,dist_CGA,dist_URRUG)
  
  ##### Plot #####
  
  # deltaG per group
  fig_deltaG_per_group <- ggplot(deltaG_table, aes(x = dataset, y = deltaG)) +
    geom_boxplot() +
    theme_mycopore() +
    scale_y_reverse() +
    xlab("") +
    ylab("deltaG (kcal/mol)") +
    ggtitle(paste0("deltaG comparison - ", input_species)) 
  
  print(fig_deltaG_per_group)
  final_list$fig_deltaG_per_group <- fig_deltaG_per_group
  
  # distance per group
  fig_dist_per_group <- ggplot(dist_table, aes(x = dataset, y = dist)) +
    geom_boxplot() +
    theme_mycopore() +
    scale_y_reverse() +
    coord_flip() +
    xlab("") +
    ylab("distance from TTS (nt)") +
    ggtitle(paste0("RBS distance comparison - ", input_species)) 
  
  print(fig_dist_per_group)
  final_list$fig_dist_per_group <- fig_dist_per_group
  
  final_list$NUGA <- NUGA_list
  final_list$noRBS <- norbs_list
  final_list$yesRBS <- yesrbs_list
  final_list$CGA <- CGA_list
  final_list$URRUG <- URRUG_list
  
  #return
  return(final_list)
  
}

# nt context of TeRe overlaps (and all start and end)
analyse_nt <- function(input_subset,input_import){
  print("analysing nt content around overlaps...")
  # check the nt context around NUGA?
  
  gff_annot <- input_subset[["full"]]
  overlaps_all <- input_subset[["5"]]
  overlaps_4_5 <- input_subset[["4nt_5"]]
  overlaps_1_5 <- input_subset[["1nt_5"]]
  fasta <- input_import$fasta
  
  pics <- list()
  
  #### calc GC ####
  num_g <- str_count(toString(input_import$fasta), "G")
  num_c <- str_count(toString(input_import$fasta), "C")
  num_a <- str_count(toString(input_import$fasta), "A")
  num_t <- str_count(toString(input_import$fasta), "T")
  
  full <- str_length(toString(input_import$fasta))
  
  ratio_g <- num_g/full
  ratio_c <- num_c/full
  ratio_a <- num_a/full
  ratio_t <- num_t/full
  
  gc <- ((num_g + num_c)/full)
  gc_percent <- gc*100
  
  print(paste0("GC% calculated as ",gc_percent,"%"))
  
  #### get nt per pos ####
  
  print("Acquiring nt/position...")
  
  # filter to avoid errors
  limit <- nchar(toString(subset(fasta)))
  gff_annot <- gff_annot %>%
    filter(start > nt_range & end < limit - nt_range)
  
  # Grab + and -
  codon_test_plus <- gff_annot %>%
    dplyr::filter(strand == "+" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  codon_test_minus <- gff_annot %>%
    dplyr::filter(strand == "-" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  print(paste0("+ strand loop: ",2*nt_range+1," iterations..."))
  for(i in -nt_range:nt_range){    # + strand
    if(i%%10 == 0){
      print(as.character(i))
    }
    k <- paste0("nt",i)
    codon_test_plus <- codon_test_plus %>%
      dplyr::mutate(!!k := toString(subseq(fasta,start+i,start+i)))
  }
  
  codon_test_plus <- codon_test_plus %>%
    filter(start > 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(subseq(fasta,start-2,start)) == "CGA",T,F))
  
  print(paste0("- strand loop: ",2*nt_range+1," iterations..."))
  for(i in -nt_range:nt_range){    # - strand
    if(i%%10 == 0){
      print(as.character(i))
    }
    k <- paste0("nt",i)
    codon_test_minus <- codon_test_minus %>%
      dplyr::mutate(!!k := toString(reverseComplement(subseq(fasta,end-i,end-i))))
  }
  
  codon_test_minus <- codon_test_minus %>%
    filter(end < limit - 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(reverseComplement(subseq(fasta,end,end+2))) == "CGA",T,F))
  
  # Merge
  codon_test <- rbind(codon_test_plus,codon_test_minus) %>%
    ungroup()
  
  # nt string
  codon_test <- codon_test %>%
    unite(col = "full_nt",!!nt_positions,sep="",remove = FALSE, na.rm=TRUE)
  
  print("merged and acquired full string!")
  
  # At stop codons to see 4th nt effect?
  print("checking nt at stop codons...")
  # Grab + and -
  stop_codon_test_plus <- gff_annot %>%
    dplyr::filter(strand == "+" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  stop_codon_test_minus <- gff_annot %>%
    dplyr::filter(strand == "-" & type == "CDS") %>%
    select(locus_name,strand,start,end, overlap_prev, has_RBS)
  
  print("+ for loop")
  for(i in -nt_range:nt_range){    # + strand
    k <- paste0("nt",i)
    stop_codon_test_plus <- stop_codon_test_plus %>%
      dplyr::mutate(!!k := toString(subseq(fasta,end+i,end+i)))
  }
  
  stop_codon_test_plus <- stop_codon_test_plus %>%
    filter(start > 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(subseq(fasta,end-2,end)) == "CGA",T,F))
  
  print("- for loop")
  for(i in -nt_range:nt_range){    # - strand
    k <- paste0("nt",i)
    stop_codon_test_minus <- stop_codon_test_minus %>%
      dplyr::mutate(!!k := toString(reverseComplement(subseq(fasta,start-i,start-i))))
  }
  
  stop_codon_test_minus <- stop_codon_test_minus %>%
    filter(end < limit - 2) %>%
    dplyr::mutate(test_CGA = ifelse(toString(reverseComplement(subseq(fasta,start,start+2))) == "CGA",T,F))
  
  # Merge
  stop_codon_test <- rbind(stop_codon_test_plus,stop_codon_test_minus) %>%
    ungroup()
  
  # nt string
  stop_codon_test <- stop_codon_test %>%
    unite(col = "full_nt",!!nt_positions,sep="",remove = FALSE, na.rm=TRUE)
  
  print("done!")
  # Subset:
    # Select overlap only
  codon_test_overlap <- codon_test %>%
    dplyr::filter(locus_name %in% overlaps_all$locus_name & overlap_prev == TRUE)
  
   # Select 4-overlap only
  codon_test_overlap_4 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_4_5$locus_name)
  
   # Select 4-overlap, no-RBS only
  codon_test_overlap_4_noRBS <- codon_test_overlap_4 %>%
    dplyr::filter(has_RBS == FALSE)
  
   # Select 4-overlap, ending in CGA only
  codon_test_overlap_4_CGA <- codon_test_overlap_4 %>%
    dplyr::filter(test_CGA == TRUE)
  
  # Select 1-overlap only
  codon_test_overlap_1 <- codon_test_overlap %>%
    dplyr::filter(locus_name %in% overlaps_1_5$locus_name)
  
  print("subsetting complete!")
  
  # sum to plot - all
  
  print("summarising and creating figures...")
  count_nuc <- codon_test %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc <- count_nuc %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc$pos <- sum_nuc$variable %>% str_replace("nt","")
  sum_nuc$pos <- as.numeric(sum_nuc$pos)
  
  fig_nt_all <- ggplot(sum_nuc,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in all - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T, colours$grey))+
    theme(legend.title = element_blank())
  
  print(fig_nt_all)
  pics$nt_all <- fig_nt_all
  
  # sum to plot - overlap
  count_nuc_overlap <- codon_test_overlap %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test_overlap))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_overlap <- count_nuc_overlap %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_overlap$pos <- sum_nuc_overlap$variable %>% str_replace("nt","")
  sum_nuc_overlap$pos <- as.numeric(sum_nuc_overlap$pos)
  
  fig_nt_overlap <- ggplot(sum_nuc_overlap,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in overlap - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_overlap)
  pics$nt_overlap <- fig_nt_overlap
  
  # sum to plot - NUGA
  count_nuc_overlap_4 <- codon_test_overlap_4 %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test_overlap_4))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_overlap_4 <- count_nuc_overlap_4 %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_overlap_4$pos <- sum_nuc_overlap_4$variable %>% str_replace("nt","")
  sum_nuc_overlap_4$pos <- as.numeric(sum_nuc_overlap_4$pos)
  
  fig_nt_overlap_4 <- ggplot(sum_nuc_overlap_4,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in NUGA - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_overlap_4)
  pics$nt_overlap_4 <- fig_nt_overlap_4
  
  # sum to plot - noRBS
  count_nuc_overlap_4_norbs <- codon_test_overlap_4_noRBS %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test_overlap_4_noRBS))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_overlap_4_norbs <- count_nuc_overlap_4_norbs %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_overlap_4_norbs$pos <- sum_nuc_overlap_4_norbs$variable %>% str_replace("nt","")
  sum_nuc_overlap_4_norbs$pos <- as.numeric(sum_nuc_overlap_4_norbs$pos)
  
  fig_nt_overlap_4_norbs <- ggplot(sum_nuc_overlap_4_norbs,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in no RBS - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_overlap_4_norbs)
  pics$nt_overlap_4_norbs <- fig_nt_overlap_4_norbs
  
  # sum to plot - CGA
  count_nuc_overlap_4_CGA <- codon_test_overlap_4_CGA %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test_overlap_4_CGA))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_overlap_4_CGA <- count_nuc_overlap_4_CGA %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_overlap_4_CGA$pos <- sum_nuc_overlap_4_CGA$variable %>% str_replace("nt","")
  sum_nuc_overlap_4_CGA$pos <- as.numeric(sum_nuc_overlap_4_CGA$pos)
  
  fig_nt_overlap_4_CGA <- ggplot(sum_nuc_overlap_4_CGA,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in CGA - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_overlap_4_CGA)
  pics$nt_overlap_4_CGA <- fig_nt_overlap_4_CGA
  
  # sum to plot - URRUG
  count_nuc_overlap_1 <- codon_test_overlap_1 %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test_overlap_1))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_overlap_1 <- count_nuc_overlap_1 %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_overlap_1$pos <- sum_nuc_overlap_1$variable %>% str_replace("nt","")
  sum_nuc_overlap_1$pos <- as.numeric(sum_nuc_overlap_1$pos)
  
  fig_nt_overlap_1 <- ggplot(sum_nuc_overlap_1,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution in URRUG - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_overlap_1)
  pics$nt_overlap_1 <- fig_nt_overlap_1
  
  # sum to plot - stop codons all
  count_nuc_stop <- stop_codon_test %>%
    select(starts_with("nt")) %>%
    gather(position,nuc) %>%
    group_by(position,nuc) %>%
    summarise(no = n()) %>%
    spread(position,no) %>%
    mutate_at(vars(-nuc),funs(. / nrow(codon_test))) %>%
    select(nuc,!!nt_positions)
  
  sum_nuc_stop <- count_nuc_stop %>%
    reshape2::melt() %>%
    mutate(value = ifelse(is.na(value),0,value))
  sum_nuc_stop$pos <- sum_nuc$variable %>% str_replace("nt","")
  sum_nuc_stop$pos <- as.numeric(sum_nuc$pos)
  
  fig_nt_stop <- ggplot(sum_nuc_stop,aes(x=pos,y=value*100,colour=nuc)) +
    geom_point() +
    geom_line() +
    theme_mycopore() +
    ggtitle(paste0("nt distribution at stops - ",input_subset$species)) +
    xlab("nt position") +
    ylab("% of total nt") +
    scale_x_continuous(breaks = seq(-30,30,5)) +
    scale_colour_manual(values = c(colours$nt_A,colours$nt_C,colours$nt_G,colours$nt_T))+
    theme(legend.title = element_blank())
  
  print(fig_nt_stop)
  pics$nt_stop <- fig_nt_stop
  
  print("done!")
  # return
  return_list <- c(species = input_subset$species)
  return_list$codon_test <- codon_test
  return_list$stop_codon_test <- stop_codon_test
  return_list$nt <- count_nuc
  return_list$fig <- pics
  
  print("nt analysis complete, returning list")
  
  return(return_list)
}

# Check proline distribution in 1/4nt overlap last10?
check_multi_P_seq <- function(input_subset,input_last10){
  # get inputs
  print("analysing proline distribution...")
  
  P_list <- c("CCT","CCC","CCG","CCA")
  
  codon_test <- input_last10$codon_test
  codon_test_NUGA <- codon_test %>%
    filter(locus_name %in% input_subset$"4nt_3"[["locus_name"]])
  codon_test_URRUG <- codon_test %>%
    filter(locus_name %in% input_subset$"1nt_3"[["locus_name"]])
  
  codon_test_str <- codon_test %>%
    filter(has_P_stretch == T)
  codon_test_NUGA_str <- codon_test_NUGA %>%
    filter(has_P_stretch == T)
  codon_test_URRUG_str <- codon_test_URRUG %>%
    filter(has_P_stretch == T)
  
  
  # calc number of all codons
  if(nrow(codon_test) > 5){
    P_test_full <- codon_test %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_full <- placeholder
  }
  
  if(nrow(P_test_full < 4)){
    P_test_full <- placeholder
  }
  
  if(nrow(codon_test_NUGA) > 5){
    P_test_NUGA <- codon_test_NUGA %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_NUGA <- placeholder
  }
  
  if(nrow(P_test_NUGA < 4)){
    P_test_NUGA <- placeholder
  }
  
  if(nrow(codon_test_URRUG) > 5){
    P_test_URRUG <- codon_test_URRUG %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_URRUG <- placeholder
  }
  
  
  if(nrow(P_test_URRUG < 4)){
    P_test_URRUG <- placeholder
  }
  
  
  if(nrow(codon_test_str) > 5){
    P_test_full_str <- codon_test_str %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_full_str <- placeholder
  }

  if(nrow(P_test_full_str < 4)){
    P_test_full_str <- placeholder
  }
  
  if(nrow(codon_test_NUGA_str) > 5){
    P_test_NUGA_str <- codon_test_NUGA_str %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_NUGA_str <- placeholder
  }
  
  if(nrow(P_test_NUGA_str < 4)){
    P_test_NUGA_str <- placeholder
  }
  
  if(nrow(codon_test_URRUG_str) > 5){
    P_test_URRUG_str <- codon_test_URRUG_str %>%
      select(contains("codon")) %>%
      gather(pos,triplet) %>%
      dplyr::group_by(pos,triplet) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      spread(pos,n) %>%
      select(-codon11) %>%
      filter(triplet %in% P_list) %>%
      replace(is.na(.),0) %>%
      rowwise %>%
      mutate(total = sum(c_across(2:11))) %>%
      ungroup() %>%
      mutate(ratio = total/sum(total))
  } else {
    P_test_URRUG_str <- placeholder
  }
  
  if(nrow(P_test_URRUG_str < 4)){
    P_test_URRUG_str <- placeholder
  }
  
  print("prolines counted!")
  
  
  # plot prep#
  sum_P <- tibble(triplet = P_test_full$triplet)
  sum_P$full_rat <- P_test_full$ratio
  sum_P$NUGA_rat <- P_test_NUGA$ratio
  sum_P$URRUG_rat <- P_test_URRUG$ratio
  
  melt_P <- sum_P %>%
    reshape2::melt()
  
  fig_P_dist <- ggplot(melt_P, aes(x=triplet,y=value,colour=variable,fill=variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_mycopore() +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(paste0("Proline ratios last 10 - ",input_subset$species)) +
    xlab("triplet") +
    ylab("relative abundance") +
    scale_fill_manual(values=c(colours$azure,colours$purple,colours$red)) +
    scale_colour_manual(values=c(colours$azure,colours$purple,colours$red)) +
    theme(legend.title = element_blank())
  
  print(fig_P_dist)
  
  print("figures drawn!")
  
  return_list <- list()
  return_list$full <- P_test_full
  return_list$NUGA <- P_test_NUGA
  return_list$URRUG <- P_test_URRUG
  return_list$full_str <- P_test_full_str
  return_list$NUGA_str <- P_test_NUGA_str
  return_list$URRUG_str <- P_test_URRUG_str
  return_list$fig_P_dist <- fig_P_dist
  
  print("analysis complete, returning list")
  
  return(return_list)
}

# Analyse NUGA/URRUG GO enrichment using COG
perform_GOA <- function(input_subset, myco_fix = F){
  print("performing GO analysis based on COG...")
  codon_test <- input_subset$full
  
  return_list <- list()
  
  # Get COG
  if(myco_fix == T){
    codon_test$COG <- COG_full$COG_ID[match(x = codon_test$Rv_name, table=COG_full$Gene_ID)]
  }
  else{
    codon_test <- codon_test %>%
      #dplyr::mutate(test_name = str_split_fixed(str_split_fixed(attributes, ";Function=", 2)[,1], "Name=", 2)[,2],) %>%
      dplyr::mutate(test_name = str_extract(attributes, "[A-Za-z0-9_]+$"))
    codon_test$COG <- COG_full$COG_ID[match(x = codon_test$test_name, table=COG_full$Gene_ID)]
  }
  codon_test$COG_str <- COG_def$func_cat[match(x = codon_test$COG, table=COG_def$COG_ID)]
  codon_test$COG_main <- substring(codon_test$COG_str,1,1)
  codon_test$COG_main_desc <- COG_functions$desc[match(x = codon_test$COG_main, table=COG_functions$func_letter)]
  
  codon_test <- codon_test %>%
    mutate(COG_main_desc = ifelse(is.na(COG_main_desc),"Not in COG database",COG_main_desc))
  
  print("COG tables acquired")
  
  # subset after annotating from GO
  codon_test_NUGA <- codon_test %>%
    filter(locus_name %in% input_subset$"4nt_3"[["locus_name"]])
  
  codon_test_noRBS <- codon_test %>%
    filter(locus_name %in% input_subset$"4nt_norbs_3"[["locus_name"]])
  
  codon_test_yesRBS <- codon_test %>%
    filter(locus_name %in% input_subset$"4nt_yesrbs_3"[["locus_name"]])
  
  codon_test_URRUG <- codon_test %>%
    filter(locus_name %in% input_subset$"1nt_3"[["locus_name"]])
  
  # sum
  full_sum <- codon_test %>%
    group_by(COG_main_desc) %>%
    summarise(n = length(COG_main_desc))
  
  NUGA_sum <- codon_test_NUGA %>%
    group_by(COG_main_desc) %>%
    summarise(n_NUGA = length(COG_main_desc))
  
  final_sum <- merge(full_sum, NUGA_sum)
  
  noRBS_sum <- codon_test_noRBS %>%
    group_by(COG_main_desc) %>%
    summarise(n_noRBS = length(COG_main_desc))
  
  final_sum <- merge(final_sum, noRBS_sum)
  
  yesRBS_sum <- codon_test_yesRBS %>%
    group_by(COG_main_desc) %>%
    summarise(n_yesRBS = length(COG_main_desc))
  
  final_sum <- merge(final_sum, yesRBS_sum)
  
  URRUG_sum <- codon_test_URRUG %>%
    group_by(COG_main_desc) %>%
    summarise(n_URRUG = length(COG_main_desc))
  
  final_sum <- merge(final_sum, URRUG_sum)
  
  final_sum <- final_sum %>%
    mutate(k = n/sum(n),
           k_NUGA = n_NUGA/sum(n_NUGA),
           k_noRBS = n_noRBS/sum(n_noRBS),
           k_yesRBS = n_yesRBS/sum(n_yesRBS),
           k_URRUG = n_URRUG/sum(n_URRUG))
  
  print("data subset and summarised")
  
  #pval,l2fc
  print("performing stats and calculating l2fc...")
  
  final_sum <- final_sum %>%
      mutate(l2fc_NUGA = log2(k_NUGA/k),
             l2fc_noRBS = log2(k_noRBS/k),
             l2fc_yesRBS = log2(k_yesRBS/k),
             l2fc_URRUG = log2(k_URRUG/k)) %>%
    mutate(p_NUGA = phyper(n_NUGA,n,sum(n) - sum(n_NUGA),sum(n_NUGA), lower.tail=F),
           p_noRBS = phyper(n_noRBS,n,sum(n) - sum(n_noRBS),sum(n_noRBS), lower.tail=F),
           p_yesRBS = phyper(n_yesRBS,n,sum(n) - sum(n_yesRBS),sum(n_yesRBS), lower.tail=F),
           p_URRUG = phyper(n_URRUG,n,sum(n) - sum(n_URRUG),sum(n_URRUG), lower.tail=F)) %>%
    mutate(padj_NUGA = p.adjust(p_NUGA, method="BH"),
           padj_noRBS = p.adjust(p_noRBS, method="BH"),
           padj_yesRBS = p.adjust(p_yesRBS, method="BH"),
           padj_URRUG = p.adjust(p_URRUG, method="BH"))
  
  print("done!")
  
  
  print("plotting...")
  # plots
  fig_COG_ratio <- final_sum %>%
    select(1,8:11) %>%
    reshape2::melt() %>%
    ggplot(aes(x=COG_main_desc,y=value, colour=variable,fill=variable)) +
    geom_bar(stat= "identity", position = "dodge") +
    theme_mycopore() +
    xlab("") +
    ylab("% of all ORFs") +
    ggtitle("Percentage of ORFs in cat") +
    scale_fill_manual(values = c(colours$red,colours$purple,colours$azure,colours$blue,colours$grey)) +
    scale_colour_manual(values = c(colours$red,colours$purple,colours$azure,colours$blue,colours$grey)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_flip() +
    theme(legend.title = element_blank())
  
  print(fig_COG_ratio)
  
  fig_l2fc <- final_sum %>%
    select(1,12:15) %>%
    reshape2::melt() %>%
    mutate(value = ifelse(value == -Inf,0,value)) %>%
    ggplot(aes(x=COG_main_desc,y=value, colour=variable,fill=variable)) +
    geom_bar(stat= "identity", position = "dodge") +
    theme_mycopore() +
    xlab("") +
    ylab("log2-fold enrichment") +
    ggtitle("COG enrichment log2-fold") +
    scale_fill_manual(values = c(colours$purple, colours$azure,colours$blue,colours$grey)) +
    scale_colour_manual(values = c(colours$purple, colours$azure,colours$blue,colours$grey)) +
    scale_y_continuous(expand=c(0,0)) +
    coord_flip() +
    theme(legend.title = element_blank())
  
  print(fig_l2fc)
  
  fig_volcano_NUGA <- ggplot(final_sum, 
                           aes(x = l2fc_NUGA, y = -log10(padj_NUGA),
                               color=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1),
                               fill=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("volcano - enrichment for NUGA"))
  
  
  print(fig_volcano_NUGA)
  
  fig_volcano_noRBS <- ggplot(final_sum, 
                             aes(x = l2fc_noRBS, y = -log10(padj_noRBS),
                                 color=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1),
                                 fill=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("volcano - enrichment for no RBS"))
  
  
  print(fig_volcano_noRBS)
  
  fig_volcano_yesRBS <- ggplot(final_sum, 
                             aes(x = l2fc_yesRBS, y = -log10(padj_yesRBS),
                                 color=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1),
                                 fill=ifelse(padj_NUGA<0.05,"grey60",colours$alex_R1))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("volcano - enrichment for yes RBS"))
  
  
  print(fig_volcano_yesRBS)
  
  fig_volcano_URRUG <- ggplot(final_sum, 
                             aes(x = l2fc_URRUG, y = -log10(padj_URRUG),
                                 color=ifelse(padj_URRUG<0.05,"grey60",colours$alex_R1),
                                 fill=ifelse(padj_URRUG<0.05,"grey60",colours$alex_R1))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_manual(values = c("grey60",colours$alex_R3)) +
    theme(legend.position = "none") +
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle(paste0("volcano - enrichment for URRUG"))
  
  
  print(fig_volcano_URRUG)
  
  print("plots finished!")
  
  # return
  return_list$full <-codon_test
  return_list$NUGA <- codon_test_NUGA
  return_list$noRBS <- codon_test_noRBS
  return_list$yesRBS <- codon_test_yesRBS
  return_list$URRUG <- codon_test_URRUG
  return_list$sum <- final_sum
  return_list$fig_ratio <- fig_COG_ratio
  return_list$fig_l2fc <- fig_l2fc
  return_list$fig_volcano_NUGA <- fig_volcano_NUGA
  return_list$fig_volcano_noRBS <- fig_volcano_noRBS
  return_list$fig_volcano_yesRBS <- fig_volcano_yesRBS
  return_list$fig_volcano_URRUG <- fig_volcano_URRUG
  
  print("COG analysis complete, returning list")
  
  return(return_list)
  
}



# Save all function - save last10 and draw_figs output into save_path

save_all <- function(annotated,subset,figs,last10,first5,RER,stop,SD_QC,nts,multiP,GOA,...){

  print("saving all data...")
  
  #Get species from last10
  species <- last10$species
    
  # Make dir if not exists
  if(!dir.exists(file.path(save_path,species))){
    dir.create(file.path(save_path,species))
    print(paste0("creating directory for ", species))
  }
  
  
  # Save - annotated
  fwrite(annotated$gff_annot,file=paste0(save_path,"/",species,"/",species,"_","annot.csv"))
  print("annot saved")
  
  # Sace - deltaG fig
  ggsave(paste0(species,"_","deltaG_dist.png"),annotated$fig_dist_deltaG,path=paste0(save_path,species,"/"))
  print("deltaG QC fig saved")
  
  # Save - subsets: 3' all, 4nt, 4nt RBS, 1nt
  fwrite(subset$"3",file=paste0(save_path,"/",species,"/",species,"_","overlaps.csv"))
  fwrite(subset$"4nt_3",file=paste0(save_path,"/",species,"/",species,"_","NUGA_overlaps.csv"))
  fwrite(subset$"4nt_norbs_3",file=paste0(save_path,"/",species,"/",species,"_","noRBS_overlaps.csv"))
  fwrite(subset$"1nt_3",file=paste0(save_path,"/",species,"/",species,"_","URRUG_overlaps.csv"))
  print("subsets saved")
  
  # Save - figs from draw_graphs
  
  ggsave(paste0(species,"_","length_dist.png"),figs$fig_l_dist,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","NUGA_pie.png"),figs$fig_NTGA_dist,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","NUGA_noRBS_pie.png"),figs$fig_NTGA_dist_noRBS,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","start_codons_pie.png"),figs$fig_start_codon,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","deltaG_per_start.png"),figs$fig_deltaG_per_group,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","dist_per_start.png"),figs$fig_dist_per_group,path=paste0(save_path,species,"/"))
  print("subset figures saved!")
  
  # Save - tables from last10
  
  fwrite(last10$codon_test,file=paste0(save_path,"/",species,"/",species,"_","codon_test.csv"))
  fwrite(last10$aa_test,file=paste0(save_path,"/",species,"/",species,"_","aa_test.csv"))
  print("last10 data saved")
  
  # Save - figures from last10 (for loop)
  for(fig in 1:length(last10[["figures"]])){
    #iterate through them, print
    figure <- last10[["figures"]][fig]
    name <- names(last10[["figures"]][fig])
    fig_name <- gsub("fig",species,name)
    print(figure)
    ggsave(filename = paste0(fig_name,".png"),path=paste0(save_path,species,"/"))
    
  }
  print("last10 figures saved")
  
  # Save - tables from first5
  fwrite(first5$codon_test,file=paste0(save_path,"/",species,"/",species,"_","first5_codon_test.csv"))
  fwrite(first5$aa_test,file=paste0(save_path,"/",species,"/",species,"_","first5_aa_test.csv"))
  print("first5 data saved")
  
  # Save - figures from first5 (for loop)
  for(fig in 1:length(first5[["pics"]])){
    #iterate through them, print
    figure <- first5[["pics"]][fig]
    name <- names(first5[["pics"]][fig])
    fig_name <- gsub("fig",species,name)
    print(figure)
    ggsave(filename = paste0("first5_",fig_name,".png"),path=paste0(save_path,species,"/"))
    
  }
  print("first5 figures saved")
  
  # Save SD_QC figs
  ggsave(paste0(species,"_","deltaG_comp.png"),SD_QC$fig_deltaG_per_group,path=paste0(save_path,species,"/"))
  ggsave(paste0(species,"_","dist_comp.png"),SD_QC$fig_dist_per_group,path=paste0(save_path,species,"/"))
  
  
  # Save - RER
  fwrite(RER$final_NUGA,file=paste0(save_path,"/",species,"/",species,"_","RER_NUGA.csv"))
  fwrite(RER$final_URRUG,file=paste0(save_path,"/",species,"/",species,"_","RER_URRUG.csv"))
  ggsave(paste0(species,"_","last_codons_NUGA.png"),RER$fig_occ_NUGA,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","last_codons_URRUG.png"),RER$fig_occ_URRUG,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","RER_NUGA.png"),RER$fig_RER,path=paste0(save_path,species,"/"), height = 10)
  ggsave(paste0(species,"_","RER_volcano_NUGA.png"),RER$fig_volcano,path=paste0(save_path,species,"/"), height = 10)
  ggsave(paste0(species,"_","RER_URRUG.png"),RER$fig_RER_URRUG,path=paste0(save_path,species,"/"), height = 10)
  ggsave(paste0(species,"_","RER_volcano_URRUG.png"),RER$fig_volcano_URRUG,path=paste0(save_path,species,"/"), height = 10)
  
  print("RER data/figures saved")
  
  # Save - stop codons
  fwrite(stop$table,file=paste0(save_path,"/",species,"/",species,"_","stop_codons.csv"))
  ggsave(paste0(species,"_","stop_codons_bar.png"),stop$figure_bar,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","stop_codons_pie.png"),stop$figure_pie,path=paste0(save_path,species,"/"),height = 10)
  
  print("stop codons saved")
  
  # Save - nt analysis tables
  fwrite(nts$codon_test,file=paste0(save_path,"/",species,"/",species,"_","codon_test_nts.csv"))
  fwrite(nts$stop_codon_test,file=paste0(save_path,"/",species,"/",species,"_","stop_codon_test_nts.csv"))
  fwrite(nts$nt,file=paste0(save_path,"/",species,"/",species,"_","nts_per_position.csv"))
  
  print("nt data saved")
  
  # Save - figures from nt (for loop)
  for(fig in 1:length(nts[["fig"]])){
    #iterate through them, print
    figure <- nts[["fig"]][fig]
    name <- names(nts[["fig"]][fig])
    fig_name <- paste0(gsub("nt",species,name),"_nt_fig")
    print(figure)
    ggsave(filename = paste0(fig_name,".png"),path=paste0(save_path,species,"/"))
    
  }
  print("nt figures saved")
  
  # Save - multiP analysis
  fwrite(multiP$full,file=paste0(save_path,"/",species,"/",species,"_","multiP_full.csv"))
  fwrite(multiP$NUGA,file=paste0(save_path,"/",species,"/",species,"_","multiP_NUGA.csv"))
  fwrite(multiP$URRUG,file=paste0(save_path,"/",species,"/",species,"_","multiP_URRUG.csv"))
  fwrite(multiP$full_str,file=paste0(save_path,"/",species,"/",species,"_","multiP_full_str.csv"))
  fwrite(multiP$NUGA_str,file=paste0(save_path,"/",species,"/",species,"_","multiP_NUGA_str.csv"))
  fwrite(multiP$URRUG_str,file=paste0(save_path,"/",species,"/",species,"_","multiP_URRUG_str.csv"))
  ggsave(paste0(species,"_","P_ratios.png"),multiP$fig_P_dist,path=paste0(save_path,species,"/"),height = 10)
  
  print("multiP data/figures saved")
  
  # Save - COG/GO analysis
  fwrite(GOA$full,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_full.csv"))
  fwrite(GOA$NUGA,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_NUGA.csv"))
  fwrite(GOA$noRBS,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_noRBS.csv"))
  fwrite(GOA$yesRBS,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_yesRBS.csv"))
  fwrite(GOA$URRUG,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_URRUG.csv"))
  fwrite(GOA$sum,file=paste0(save_path,"/",species,"/",species,"_","COG_GOA_sum.csv"))
  
  ggsave(paste0(species,"_","COG_percents.png"),GOA$fig_ratio,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","COG_enrich.png"),GOA$fig_l2fc,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","COG_NUGA_volcano.png"),GOA$fig_volcano_NUGA,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","COG_noRBS_volcano.png"),GOA$fig_volcano_noRBS,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","COG_yesRBS_volcano.png"),GOA$fig_volcano_yesRBS,path=paste0(save_path,species,"/"),height = 10)
  ggsave(paste0(species,"_","COG_URRUG_volcano.png"),GOA$fig_volcano_URRUG,path=paste0(save_path,species,"/"),height = 10)
  
  print("COG/GOA data and figures saved")
  
  # finish up
  print("all data is saved! Yay!")
  print("find your data at:")
  print(paste0(save_path,species))
}

# Full analysis function
analyse_overlaps <- function(species = "Mtb", processed_codons = F, gff_type = 3,
                             tasks=c("acid","base","rare","nP","strP","Arg"),
                             multi_chromosome=F, save_files = T,
                             do_GOA = T, fix_myco = F, skip_SD = F, ...){
  # Run everything in pipe
  
  file_list <- read_files(species = species, 
                          processed_codons = processed_codons,
                          multiple_chromosomes = multi_chromosome,
                          gff_type = gff_type)
  
  # Species name
  s <- deparse(substitute(species))
  s <- gsub("\"", "", s, fixed = TRUE)
  
  file_list$species <- s
  
  annotate_list <- annotate_overlaps(file_list, manual_SD = skip_SD)
  subset_list <- subset_overlaps(annotate_list)
  fig_list <- draw_graphs(subset_list)
  last10_list <- do_last10(input_subset = subset_list, input_import = file_list, input_task = tasks)
  first5_list <- do_first5(input_subset = subset_list, input_import = file_list)
  RER_list <- calc_NUGA_RER(last10_list,subset_list)
  stop_list <- show_stop(last10_list)
  SD_QC_list <- do_SD_QC(input_subset = subset_list, input_last10 = last10_list)
  nt_list <- analyse_nt(input_subset = subset_list, input_import = file_list)
  multiP_list <- check_multi_P_seq(input_subset = subset_list, input_last10 = last10_list)
  GOA_list <- perform_GOA(input_subset = subset_list, myco_fix = fix_myco)
  
  # list
  final_list <- list("files" = file_list,
                  "annot" = annotate_list,
                  "subset" = subset_list,
                  "fig" = fig_list,
                  "last10" = last10_list,
                  "first5" = first5_list,
                  "RER" = RER_list,
                  "stops" = stop_list,
                  "SD_QC_list" = SD_QC_list,
                  "nts" = nt_list,
                  "multiP " = multiP_list,
                  "GOA" = GOA_list)
  
  # Save?
  if(save_files == T){
    save_all(annotate_list,subset_list,fig_list,last10_list,first5_list,
             RER_list,stop_list,SD_QC_list,nt_list,multiP_list,GOA_list)
    
    # also save as RDS
    saveRDS(final_list,paste0(save_path,"/",s,"/",s,"_full.rds"))
    
  }
  
  return(final_list)
}


#### Main ####


##### read files list #####
q <- read_files("Mtb",processed_codons = T)
q2 <- read_files("Msm", gff_type = 2,processed_codons = T)
q3 <- read_files("Ecoli",gff_type=3,processed_codons = F)
# q4 <- read_files("vibrio",gff_type=2,processed_codons = F,multiple_chromosomes = T)
# q5 <- read_files("pseudomonas",gff_type=4,processed_codons = F)
q6 <- read_files("listeria",gff_type=4,processed_codons = F)
# q7 <- read_files("bovis",gff_type=4,processed_codons = F)
# q8 <- read_files("Cdipht",gff_type=4,processed_codons = F)
# q9 <- read_files("Cglut",gff_type=4,processed_codons = F)
# q10 <- read_files("desulfo",gff_type=4,processed_codons = F)
# q11 <- read_files("leprae",gff_type=4,processed_codons = F)
# q12 <- read_files("leptothrix",gff_type=4,processed_codons = F)
# q13 <- read_files("Mabs",gff_type=4,processed_codons = F,multiple_chromosomes = T)
# q14 <- read_files("salmonella",gff_type=4,processed_codons = F)
q15 <- read_files("staph",gff_type=4,processed_codons = F)
# q16 <- read_files("strep",gff_type=4,processed_codons = F)
q17 <- read_files("subtilis",gff_type=4,processed_codons = F)
# q18 <- read_files("Xcamp",gff_type=4,processed_codons = F)
# q19 <- read_files("Cdiff",gff_type=4,processed_codons = F)
# q20 <- read_files("Suso",gff_type=4,processed_codons = F)



##### test sequence #####
l <- annotate_overlaps(q3, manual_SD = F)
m <- subset_overlaps(l)
p <- draw_graphs(m,has_functions = T)
j <- do_last10(input_subset = m,input_import = q3, input_tasks = c("acid","base","rare","nP","strP","Arg"))
v <- do_first5(m,q)
z <- calc_NUGA_RER(j,m) #add NUGA split by RBS
y <- show_stop(j)
f <- do_SD_QC(m,j)
b <- analyse_nt(m,q) # better print annotation - for loops
h <- check_multi_P_seq(m,j)
t <- perform_GOA(m,F)

save_all(l,m,p,j,v,z,y,f,b,h,t)


##### run analysis code per species #####

x <- analyse_overlaps(species = "Mtb",processed_codons = T, gff_type = 1,do_GOA = T,fix_myco = T)
x2 <- analyse_overlaps(species = "Msm",processed_codons = T,gff_type = 2)
x3 <- analyse_overlaps(species = "Ecoli",processed_codons = F,gff_type = 3)
x4 <- analyse_overlaps(species = "vibrio",processed_codons = F,multi_chromosome=T,gff_type = 2)
x5 <- analyse_overlaps(species = "pseudomonas",processed_codons = F,gff_type = 4)
x6 <- analyse_overlaps(species = "listeria",processed_codons = F,gff_type = 4)
#x7 <- analyse_overlaps(species = "bovis",processed_codons = F,gff_type = 4)
x8 <- analyse_overlaps(species = "Cdipht",processed_codons = F,multi_chromosome=T,gff_type = 4)
x9 <- analyse_overlaps(species = "Cglut",processed_codons = F,gff_type = 4)
x10 <- analyse_overlaps(species = "desulfo",processed_codons = F,multi_chromosome=T,gff_type = 4)
x11 <- analyse_overlaps(species = "leprae",processed_codons = F,gff_type = 4)
x12 <- analyse_overlaps(species = "leptothrix",processed_codons = F,gff_type = 4)
x13 <- analyse_overlaps(species = "Mabs",processed_codons = F,multi_chromosome=T,gff_type = 4)
#x14 <- analyse_overlaps(species = "salmonella",processed_codons = F,multi_chromosome=T,gff_type = 4)
x15 <- analyse_overlaps(species = "staph",processed_codons = F,gff_type = 4)
x16 <- analyse_overlaps(species = "strep",processed_codons = F,gff_type = 4)
x17 <- analyse_overlaps(species = "subtilis",processed_codons = F,gff_type = 4,skip_SD = T)
x18 <- analyse_overlaps(species = "Xcamp",processed_codons = F,gff_type = 4)
#x19 <- analyse_overlaps(species = "Cdiff",processed_codons = F,multi_chromosome=T,gff_type = 4)
x20 <- analyse_overlaps(species = "Suso",processed_codons = F,gff_type = 4)


# Todo: marinum, avium, janaschii, acidocaldarius, campylobacter, Helicobacter pylori,
# ...chlamydia, yersinia, rickettsia, trepomena pallidum

### NEED TO REMOVE GENE_SYNONYM WITH THE BIOTYPE STR_SEARCH!!
### MERGE FASTA IN CODE NOT IN EDITOR?
### Breaks: bovis, salmonella, cdiff
### Get warnings muted 

# quick getting fig for how much more RBS detected?
# asda <- data.table("subset" = c("full","overlap","NUGA","URRUG"), "n" = 0)
# zaza <- x[["subset"]]
# asda$n[1] <- nrow(zaza[["full"]] %>% filter(has_RBS == TRUE)) /
#   (nrow(zaza[["full"]] %>% filter(has_RBS == F)) + 
#      nrow(zaza[["full"]] %>% filter(has_RBS == TRUE)))
# asda$n[2] <- nrow(zaza[["overlaps"]] %>% filter(has_RBS == TRUE)) /
#   (nrow(zaza[["overlaps"]] %>% filter(has_RBS == F)) + 
#      nrow(zaza[["overlaps"]] %>% filter(has_RBS == TRUE)))
# asda$n[3] <- nrow(zaza[["4nt_3"]] %>% filter(has_RBS == TRUE)) /
#   (nrow(zaza[["4nt_3"]] %>% filter(has_RBS == F)) + 
#      nrow(zaza[["4nt_3"]] %>% filter(has_RBS == TRUE)))
# asda$n[4] <- nrow(zaza[["1nt_3"]] %>% filter(has_RBS == TRUE)) /
#   (nrow(zaza[["1nt_3"]] %>% filter(has_RBS == F)) + 
#      nrow(zaza[["1nt_3"]] %>% filter(has_RBS == TRUE)))
# 
# asda %>%
#   ggplot(aes(x = n, y = subset, fill=subset)) +
#   geom_bar(position = "dodge", stat="identity") +
#   theme_mycopore() +
#   scale_fill_manual(values = c(colours$blue,colours$purple,colours$red, colours$orange)) +
#   xlab("") +
#   ylab("")

