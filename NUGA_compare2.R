# PhD Final Stretch: Comparisons between organisms

#   • Reads NUGA_multispecies outputs
#   • Performs analysis 
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

# Species folder path (this is output folder of NUGA_multispecies.R)
species_path <- paste0(here(),"/NUGA_Out/species/")

#How many of the most enriched codons to take per species?
n_sample_RER <- 3

# What counts as mycobacteria
myco_list <- c("Mtb","Msm","Mabs","bovis","leprae","avium","marinum")

# Is reading working?
# test_read <- readRDS(paste0(species_path,"Mtb/Mtb_full.rds"))



########
# Func #
########

# Function to read a single full_table from species name
# Calculate GC content for downstream?
read_full <- function(species, do.GC = T){
  
  if(!dir.exists(file.path(species_path,species))){
    print("you gave a wrong species!")
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

jj <- read_full("Mtb")
# Go through list of species, read all, return one gigatable
# also have list of species in return list
read_all <- function(species_list){
  return_list <- list()
  for(item in 1:length(species_list)){
    to_add <- read_full(species_list[item])
    name <- to_add$files[["species"]]
    return_list[[name]] <- to_add
  }
  return(return_list)
}

klmn <- read_all(c("Mtb","Msm","Mabs","vibrio","Ecoli",
                   "subtilis","staph","strep","leprae","Cdipht",
                   "Xcamp", "pseudomonas", "listeria", "CGlut", "desulfo",
                   "leptothrix", "Suso"))


# Select only the tables needed from gigatable
# merge them into long format tables for specific
# Return specific tables asked
process_tables <- function(all_list){
  species_list <- list()
  names <- unlist(attributes(all_list),use.names = F)
  #### RER ####
  RER_table <- tibble()
  RER2_table <- tibble()
  RER_head_table <- tibble()
  RER2_head_table <- tibble()
  for(i in seq_along(all_list)){
    RER <- all_list[[i]][["RER"]][["final_NUGA"]] %>%
      mutate(dataset = names[[i]])
    RER_short <- RER %>%
      arrange(padj) %>%
      head(n_sample_RER)
    RER2 <- all_list[[i]][["RER"]][["final_URRUG"]] %>%
      mutate(dataset = names[[i]])
    RER2_short <- RER2 %>%
      arrange(padj) %>%
      head(n_sample_RER)
    RER_table <- rbind(RER_table,RER)
    RER_head_table <- rbind(RER_head_table,RER_short)
    RER2_table <- rbind(RER2_table,RER2)
    RER2_head_table <- rbind(RER2_head_table,RER2_short)
  }
  
  # Sum for "how many times the codon occurs in these i.e. the "most likely to be in top5 enrich"
  RER_sum <- RER_head_table %>%
    group_by(triplet) %>%
    summarise(n = length(triplet))
  
  RER2_sum <- RER2_head_table %>%
    group_by(triplet) %>%
    summarise(n = length(triplet))
  
  #### stop codons ####
  
  stop_table <- tibble()
  for(i in seq_along(all_list)){
    stop_codon <- all_list[[i]][["stops"]][["table"]] %>%
      mutate(dataset = names[[i]]) %>%
      mutate(k = n/sum(n))
    stop_table <- rbind(stop_table,stop_codon)
  }
  
  #### Rare per slot ####
  rare_table <- tibble()
  rare_table_CGA <- tibble()
  rare_table_NUGA <- tibble()
  for(i in seq_along(all_list)){
    rare <- all_list[[i]][["last10"]][["rare_all"]] %>%
      filter(category == "all") %>%
      mutate(dataset = names[[i]])
    rare_CGA <- all_list[[i]][["last10"]][["rare_all"]] %>%
      filter(category == "overlap_4_CGA") %>%
      mutate(dataset = names[[i]])
    rare_NUGA <- all_list[[i]][["last10"]][["rare_all"]] %>%
      filter(category == "overlap_4") %>%
      mutate(dataset = names[[i]])
    rare_table <- rbind(rare_table,rare)
    rare_table_CGA <- rbind(rare_table_CGA,rare_CGA)
    rare_table_NUGA <- rbind(rare_table_NUGA,rare_NUGA)
  }
  
  rare_table_CGA <- rare_table_CGA %>%
    mutate_all(~replace(., is.na(.), 0))
  
  rare_table_NUGA <- rare_table_NUGA %>%
    mutate_all(~replace(., is.na(.), 0))
  #### Acid per slot ####
  acid_table <- tibble()
  acid_table_CGA <- tibble()
  acid_table_NUGA <- tibble()
  for(i in seq_along(all_list)){
    acid <- all_list[[i]][["last10"]][["acid_all"]] %>%
      filter(category == "all") %>%
      mutate(dataset = names[[i]])
    acid_CGA <- all_list[[i]][["last10"]][["acid_all"]] %>%
      filter(category == "overlap_4_CGA") %>%
      mutate(dataset = names[[i]])
    acid_NUGA <- all_list[[i]][["last10"]][["acid_all"]] %>%
      filter(category == "overlap_4") %>%
      mutate(dataset = names[[i]])
    acid_table <- rbind(acid_table,acid)
    acid_table_CGA <- rbind(acid_table_CGA,acid_CGA)
    acid_table_NUGA <- rbind(acid_table_NUGA,acid_NUGA)
  }
  
  acid_table_CGA <- acid_table_CGA %>%
    mutate_all(~replace(., is.na(.), 0))
  acid_table_NUGA <- acid_table_NUGA %>%
    mutate_all(~replace(., is.na(.), 0))
  
  #### Base per slot ####
  base_table <- tibble()
  base_table_CGA <- tibble()
  base_table_NUGA <- tibble()
  for(i in seq_along(all_list)){
    base <- all_list[[i]][["last10"]][["base_all"]] %>%
      filter(category == "all") %>%
      mutate(dataset = names[[i]])
    base_CGA <- all_list[[i]][["last10"]][["base_all"]] %>%
      filter(category == "overlap_4_CGA") %>%
      mutate(dataset = names[[i]])
    base_NUGA <- all_list[[i]][["last10"]][["base_all"]] %>%
      filter(category == "overlap_4") %>%
      mutate(dataset = names[[i]])
    base_table <- rbind(base_table,base)
    base_table_CGA <- rbind(base_table_CGA,base_CGA)
    base_table_NUGA <- rbind(base_table_NUGA,base_NUGA)
  }
  
  base_table_CGA <- base_table_CGA %>%
    mutate_all(~replace(., is.na(.), 0))
  base_table_NUGA <- base_table_NUGA %>%
    mutate_all(~replace(., is.na(.), 0))
  
  #### SD QCs ####
  QC_table <- tibble()
  for(i in seq_along(all_list)){
    NUGA_list <- all_list[[i]][["SD_QC_list"]][["NUGA"]] %>%
      ungroup %>%
      select("deltaG","RBS_dist") %>%
      mutate(dataset = names[[i]],
             type = "NUGA")
    URRUG_list <- all_list[[i]][["SD_QC_list"]][["URRUG"]] %>%
      ungroup %>%
      select("deltaG","RBS_dist") %>%
      mutate(dataset = names[[i]],
             type = "URRUG")
    noRBS_list <- all_list[[i]][["SD_QC_list"]][["noRBS"]] %>%
      ungroup %>%
      select("deltaG","RBS_dist") %>%
      mutate(dataset = names[[i]],
             type = "noRBS")
    yesRBS_list <- all_list[[i]][["SD_QC_list"]][["yesRBS"]] %>%
      ungroup %>%
      select("deltaG","RBS_dist") %>%
      mutate(dataset = names[[i]],
             type = "yesRBS")
    CGA_list <- all_list[[i]][["SD_QC_list"]][["CGA"]] %>%
      ungroup %>%
      select("deltaG","RBS_dist") %>%
      mutate(dataset = names[[i]],
             type = "CGA")
    QC_table <- rbind(QC_table,NUGA_list,URRUG_list,noRBS_list,yesRBS_list,CGA_list)
  }
  
  #### GOA ####
  GOA_table <- tibble()
  for(i in seq_along(all_list)){
    GOA_list <- all_list[[i]][["GOA"]][["sum"]] %>%
      ungroup %>%
      mutate(dataset = names[[i]])
    GOA_table <- rbind(GOA_table,GOA_list)
  }
  
  
  #### Save GC ####
  GC_table <- tibble("species" = names)
  GC_col <- c()
  for(i in seq_along(all_list)){
    GC <- all_list[[i]][["files"]][["GC"]]
    GC_col <- append(GC_col,GC)
  }
  
  GC_table$n <- GC_col
  
  #### Return all ####
  return_list <- list(
    "species" = names,
    "RER" = RER_table,
    "RER_URRUG" = RER2_table,
    "RER_head" = RER_head_table,
    "RER_head_URRUG" = RER2_head_table,
    "RER_sum" = RER_sum,
    "RER_sum_URRUG" = RER2_sum,
    "stops" = stop_table,
    "rare" = rare_table,
    "rare_CGA" = rare_table_CGA,
    "rare_NUGA" = rare_table_NUGA,
    "acid" = acid_table,
    "acid_CGA" = acid_table_CGA,
    "acid_NUGA" = acid_table_NUGA,
    "base" = base_table,
    "base_CGA" = base_table_CGA,
    "base_NUGA" = base_table_NUGA,
    "GC" = GC_table,
    "QC" = QC_table,
    "GOA" = GOA_table
  )
  
  return(return_list)

}

l <- process_tables(klmn)

# Plot based on processed
plot_merged <- function(processed){
  
  #### Prep ####
  
  # Get tables
  RER_head_table <- processed$RER_head %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  RER_table <- processed$RER %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  RER2_head_table <- processed$RER_head_URRUG %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  RER2_table <- processed$RER_URRUG %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  stop_table <- processed$stops %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  rare_table <- processed$rare %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  rare_table_CGA <- processed$rare_CGA %>%
    select(2:10,12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  rare_table_NUGA <- processed$rare_NUGA %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  acid_table <- processed$acid %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  acid_table_CGA <- processed$acid_CGA %>%
    select(2:10,12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  acid_table_NUGA <- processed$acid_NUGA %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  base_table <- processed$base %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  base_table_CGA <- processed$base_CGA %>%
    select(2:10,12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  base_table_NUGA <- processed$base_NUGA %>%
    select(2:12) %>%
    reshape2::melt(id = "dataset") %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F)) %>%
    group_by(is_myco)
  taa_table <- stop_table %>%
    filter(codon11 == "TAA") %>%
    dplyr::rename("species" = dataset)
  tag_table <- stop_table %>%
    filter(codon11 == "TAG") %>%
    dplyr::rename("species" = dataset)
  tga_table <- stop_table %>%
    filter(codon11 == "TGA") %>%
    dplyr::rename("species" = dataset)
  GC_table <- merge(processed$GC,taa_table,by="species") %>%
    dplyr::rename("GC" = n.x) %>%
    select(species,GC,k,is_myco)
  GC_table2 <- merge(processed$GC,tag_table,by="species") %>%
    dplyr::rename("GC" = n.x) %>%
    select(species,GC,k,is_myco)
  GC_table3 <- merge(processed$GC,tga_table,by="species") %>%
    dplyr::rename("GC" = n.x) %>%
    select(species,GC,k,is_myco)
  QC_table <- processed$QC %>%
    mutate(is_myco = ifelse(dataset %in% myco_list,T,F))
  
  #### Figure ####
  # RER per codon
  fig_RER <- ggplot(RER_head_table, aes(x = triplet, y = RERlog, fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("") +
    ylab("log2 enrichment") +
    ggtitle("Top 5 enriched NUGA last codons") +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.title = element_blank())
  
  
  print(fig_RER)
  
  # RER2 per codon
  fig_RER2 <- ggplot(RER_head_table, aes(x = triplet, y = RERlog, fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("") +
    ylab("log2 enrichment") +
    ggtitle("Top 5 enriched URRUG last codons") +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.title = element_blank())
  
  
  print(fig_RER2)
  
  # RER per codon for padj
  fig_RER_padj <- ggplot(RER_head_table, aes(x = triplet, y = -log10(padj), fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("") +
    ylab("-log10 p(adjust)") +
    ggtitle("Most significant NUGA last codons") +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.title = element_blank())
  
  
  print(fig_RER_padj)
  
  # RER2 per codon for padj
  fig_RER2_padj <- ggplot(RER2_head_table, aes(x = triplet, y = -log10(padj), fill=dataset)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("") +
    ylab("-log10 p(adjust)") +
    ggtitle("Most significant URRUG last codons") +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.title = element_blank())
  
  
  print(fig_RER2_padj)
  
  # RER heatmap?
  fig_RER_heat <- ggplot(RER_table, aes(x=reorder(dataset, -is_myco), y=triplet, fill=-log10(padj))) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$azure, na.value = colours$azure) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("")
  
  print(fig_RER_heat)
  
  # RER URRUG heatmap?
  fig_RER2_heat <- ggplot(RER2_table, aes(x=reorder(dataset, -is_myco), y=triplet, fill=-log10(padj))) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$red, na.value = colours$azure) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("")
  
  print(fig_RER2_heat)
  
  # stop codon dist
  fig_stop_dist <- ggplot(stop_table, aes(x=reorder(dataset, -is_myco), y = k*100, fill=codon11)) +
    geom_bar(position="dodge",stat="identity") +
    theme_mycopore() +
    xlab("") +
    ylab("% of stop codons") +
    ggtitle("Stop codon distribution") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values = c(colours$blue,colours$purple,colours$green)) +
    theme(legend.title = element_blank())
  
  print(fig_stop_dist)
  

  
  # volcano on top5 RER
  fig_RER_volcano <- ggplot(RER_head_table) +
    geom_point(aes(x=RERlog,
                   y=-log10(padj),
                   color=triplet,
                   fill=triplet,
                   shape=is_myco,
                   size=1.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour=colours$alex_R1) +
    theme_mycopore() +
    scale_color_brewer(palette = "Paired") +
    scale_shape_manual(values=c(16,15))+
    xlab("log2-fold enrichment") +
    ylab("-log10 p(adjusted)") +
    ggtitle("Last codon enrichment") +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill="grey95"))
  
  print(fig_RER_volcano)
  
  # heatmap for rare
  fig_rare_heat <- ggplot(rare_table, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$purple, na.value = colours$purple) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Rare per position - all")
  
  print(fig_rare_heat)
  
  fig_rare_CGA_heat <- ggplot(rare_table_CGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$darkgold, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Rare per position - CGA")
  
  print(fig_rare_CGA_heat)
  
  fig_rare_NUGA_heat <- ggplot(rare_table_NUGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$green, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Rare per position - NUGA")
  
  print(fig_rare_NUGA_heat)
  
  # heatmap for acid
  fig_acid_heat <- ggplot(acid_table, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$purple, na.value = colours$purple) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Acidic per position - all")
  
  print(fig_acid_heat)
  
  fig_acid_CGA_heat <- ggplot(acid_table_CGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$darkgold, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Acid per position - CGA")
  
  print(fig_acid_CGA_heat)
  
  fig_acid_NUGA_heat <- ggplot(acid_table_NUGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$green, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Acidic per position - NUGA")
  
  print(fig_acid_NUGA_heat)
  
  # heatmap for base
  fig_base_heat <- ggplot(base_table, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$purple, na.value = colours$purple) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Basic per position - all")
  
  print(fig_base_heat)
  
  fig_base_CGA_heat <- ggplot(base_table_CGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$darkgold, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Basic per position - CGA")
  
  print(fig_base_CGA_heat)
  
  fig_base_NUGA_heat <- ggplot(base_table_NUGA, aes(y=reorder(dataset, -is_myco), x=variable, fill=value)) +
    geom_tile() +
    scale_fill_gradient(low="white",high=colours$green, na.value = "white") +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("Position") +
    scale_x_discrete(labels=(-10:-1)) +
    ggtitle("Basic per position - NUGA")
  
  print(fig_base_NUGA_heat)
  
  # Stop codon TAA regressed against GC%?
  fig_taa_regress <- ggplot(GC_table, aes(x = k*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("TAA% as stop codons") +
    ylab("GC%") +
    ggtitle("GC content vs TAA stop codons")

  print(fig_taa_regress)
  
  # Stop codon TAG regressed against GC%?
  fig_tag_regress <- ggplot(GC_table2, aes(x = k*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("TAG% as stop codons") +
    ylab("GC%") +
    ggtitle("GC content vs TAG stop codons")
  
  print(fig_tag_regress)
  
  # Stop codon TGA regressed against GC%?
  fig_tga_regress <- ggplot(GC_table3, aes(x = k*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("TGA% as stop codons") +
    ylab("GC%") +
    ggtitle("GC content vs TGA stop codons")
  
  print(fig_tga_regress)
  
  # QC heatmap for deltaG
  fig_QC_deltaG <- ggplot(QC_table, aes(y=reorder(dataset, is_myco), x=type, fill=deltaG)) +
    geom_tile() +
    scale_fill_gradient(low=colours$orange,high="white", na.value = colours$azure) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("")
  
  print(fig_QC_deltaG)
  
  fig_QC_dist <- ggplot(QC_table, aes(y=reorder(dataset, is_myco), x=type, fill=RBS_dist)) +
    geom_tile() +
    scale_fill_gradient(low=colours$lightlime2,high="white", na.value = colours$azure) +
    theme_heatmap() +
    theme(
      axis.text.y = element_text(face="bold",size = rel(1.3)), 
      axis.line.y = element_line(colour="black",size=1),
      axis.ticks.y = element_line(colour="black",size = 1),
      legend.position = "right") +
    guides(fill = guide_legend_heatmap()) +
    xlab("")
  
  print(fig_QC_dist)
  
}

plot_merged(l)


# What% of ORFs are overlap? NUGA? NUGA noRBS? 

calc_NUGA_percents <- function(all_list,processed_list){
  species_list <- list()
  names <- unlist(attributes(all_list),use.names = F)
  
  perc_table <- tibble()
  
  
  # Get the numbers
  for(i in seq_along(all_list)){
    all <- nrow(all_list[[i]][["subset"]][["full"]])
    OL <- nrow(all_list[[i]][["subset"]][["3"]])
    NUGA <- nrow(all_list[[i]][["subset"]][["4nt_3"]])
    noRBS <- nrow(all_list[[i]][["subset"]][["4nt_norbs_3"]])
    
    #CGA needs to be grabbed separately
    CGA_table <- all_list[[i]][["last10"]][["codon_test"]] %>%
      filter(locus_name %in% all_list[[i]][["subset"]][["4nt_3"]][["locus_name"]]) %>%
      filter(codon10 == "CGA")
    CGA <- nrow(CGA_table)
    
    
    add <- tibble(n_all = all,
                  n_OL = OL,
                  n_NUGA = NUGA,
                  n_noRBS = noRBS,
                  n_CGA = CGA,
                  species = names[i])
    
    perc_table <- rbind(perc_table,add)
  }
  
  # add ratios
  perc_table <- perc_table %>%
    dplyr::mutate(r_OL = n_OL/n_all,
                  r_NUGA = n_NUGA/n_all,
                  r_noRBS = n_noRBS/n_all,
                  r_CGA_total = n_CGA/n_all,
                  r_CGA_of_NUGA = n_CGA/n_NUGA) %>%
    mutate(is_myco = ifelse(species %in% myco_list,T,F)) %>%
    group_by(is_myco)
  
  # calc the GC content related
  GC_table <-processed_list[["GC"]]
  
  perc_table <- merge(perc_table,GC_table, by="species") %>%
    rename(n = "GC")
  
  
  # melt?
  perc_melt <- perc_table %>%
    select(starts_with("r_"),"species") %>%
    select(-"r_CGA_of_NUGA") %>%
    dplyr::rename("Overlapping" = "r_OL",
                  "NUGA" = "r_NUGA",
                  "NUGA, no RBS" = "r_noRBS",
                  "NUGA with CGA" = "r_CGA_total") %>%
    reshape2::melt(id.vars = c("species")) %>%
    mutate(is_myco = ifelse(species %in% myco_list,T,F)) %>%
    group_by(is_myco)
  
  # PLOTS:
  # all %
  perc_plot <- ggplot(perc_melt, aes(x=reorder(species,-is_myco),y=value*100, colour=variable, fill = variable)) +
    geom_bar(stat= "identity", position = "dodge") +
    theme_mycopore() +
    xlab("Species") +
    ylab("% of all ORFs") +
    ggtitle("Percentage of ORFs") +
    scale_fill_manual(values = c(colours$red,colours$purple,colours$blue,colours$alex_B1_edge)) +
    scale_colour_manual(values = c(colours$red,colours$purple,colours$blue,colours$alex_B1_edge)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.title = element_blank())
  
  print(perc_plot)
  
  # CGA as part of NUGA
  CGA_NUGA_plot <- ggplot(perc_table,aes(x=reorder(species, -r_CGA_of_NUGA),y=100*r_CGA_of_NUGA, fill= is_myco)) +
    geom_bar(stat="identity") +
    theme_mycopore() +
    xlab("Species") +
    ylab("% of NUGA ORFS") +
    ggtitle("NUGA ORFs containing CGA") +
    scale_fill_manual(values = c(colours$grey,colours$red)) +
    scale_colour_manual(values = c(colours$grey,colours$red)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position = "none",
          axis.text.x = element_text(size=rel(1.6)))
  
  print(CGA_NUGA_plot)
  
  # The role of GC% - NUGA
  fig_NUGA_regress <- ggplot(perc_table, aes(x = r_NUGA*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("% ORFS with NUGA overlap") +
    ylab("GC%") +
    ggtitle("GC content vs NUGA ORFs")
  
  print(fig_NUGA_regress)
  
  # The role of GC% - CGA
  fig_CGA_regress <- ggplot(perc_table, aes(x = r_CGA_total*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("% ORFS with CGA(NUGA) overlap") +
    ylab("GC%") +
    ggtitle("GC content vs CGA ORFs")
  
  print(fig_CGA_regress)
  
  # The role of GC% - CGA
  fig_CGA_regress_2 <- ggplot(perc_table, aes(x = r_CGA_of_NUGA*100, y = GC)) +
    geom_point(aes(shape=is_myco)) +
    stat_smooth(method='lm', formula= y~x, fullrange = T, 
                color=colours$alex_B3,fill=colours$alex_B1) +
    theme_mycopore() +
    scale_x_continuous(expand=c(0,1)) +
    scale_y_continuous(expand=c(0,1)) +
    scale_shape_manual(values=c(16,15)) +
    xlab("% NUGA ORFS with CGA overlap") +
    ylab("GC%") +
    ggtitle("GC content vs CGA ORFs per NUGA")
  
  print(fig_CGA_regress_2)

  
  # return
  return(perc_table)
  
}


xyz <- calc_NUGA_percents(klmn,l)
