###########################################################################

# Additional data analysis for D'Halluin et al., 2023
#   • Note - most of this is copy-paste and not solidly defined functions,
#   • due to lack of time (and energy). Apologies in advance if you check it!

#   • p-value corrections for RT-score

###########################################################################

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","reshape2","DESeq2")
invisible(lapply(packages, require, character.only = TRUE))

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

# ggplot2 theme for figures
theme_alex <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2), hjust = 0.5, vjust=3),
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
           strip.text = element_text(face="bold"),
           aspect.ratio = 1/1
   ))
}

theme_heatmap2 <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2), hjust = 0.5),
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
           axis.ticks.length = unit(.20, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.background= element_rect(color = NA, fill = "white"),
           legend.spacing.x = unit(0.2, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="bold.italic", color = "black", size = rel(1.2)),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.8, "cm"),
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
  alex_B1_edge = "#BBD0FA",
  alex_R3_edge = "#520101",
  ashenskin = "#e8dcde"
)

# tsv loader
get_tsv_input <- function(tsv_file, norm_to) {
  input_tsv <- paste(here::here("NUGA_data/"), tsv_file, ".tsv", sep = "")
  column_name <- as.character(tsv_file)
  n <- nreads[[norm_to]]
  tsv <- read.csv(input_tsv, sep='\t') %>%
    dplyr::select(2:3) %>%
    dplyr::rename("pos"=1,!!column_name:=2)
  tsv[[column_name]] = round(((tsv[[column_name]]/n)*1000000),3) # to convert to CPM
  k=round(seq.int(0,4411532, by=1),1)
  k_table <-  tibble("pos" = k, !!column_name := 0 )
  tsv <- rbind(tsv,k_table %>% filter(pos %nin% tsv$pos)) %>%
    arrange(pos)
  return(tsv)
}


#### Main ####

nt_100 <- read.csv(here::here("NUGA_Data/100nt_table.csv")) %>%
  dplyr::select(ID,Position,Strand,Locus,Class,RT.score.0H.3H,RT.score.0H.4.5H,RT.score.0H.6H) %>%
  dplyr::rename("RT0_3" = RT.score.0H.3H,
                "RT0_45" = RT.score.0H.4.5H,
                "RT0_6" = RT.score.0H.6H) %>%
  filter(!is.na(as.numeric(RT0_3)),
         !is.na(as.numeric(RT0_45)),
         !is.na(as.numeric(RT0_6)),
         Class != "A")

rho <- read.csv(here::here("NUGA_Data/rho_table.csv")) %>%
  dplyr::mutate(LTF = ifelse(Leading == "Found", "Leading",
                             ifelse(Trailing == "Found", "Trailing", "Final"))) %>%
  dplyr::select(ID,TTS.position,Strand,Locus,Category,Class,LTF,Rho.Log2.0H.3H.,Rho.Log2.0H.4.5H.,Rho.Log2.0H.6H.) %>%
  dplyr::rename("Rho3" = Rho.Log2.0H.3H.,
                "Rho45" = Rho.Log2.0H.4.5H.,
                "Rho6" = Rho.Log2.0H.6H.) %>%
  filter(!is.na(as.numeric(Rho3)),
         !is.na(as.numeric(Rho45)),
         !is.na(as.numeric(Rho6)),
         Class != "A") 

# Old numbers, as calced by Alex
tts_old <- read.csv (here::here("NUGA_Data/TTS_score_new.csv")) %>%
  dplyr::rename("fdr3" = "X3H.FDR",
                "fdr45" = "X4.5H.FDR",
                "fdr6" = "X6H.FDR",
                "p3"= "p.val.3H",
                "p45"= "p.val.4.5H",
                "p6"= "p.val.6H") %>%
  dplyr::filter(class != "A")

# n-reads for tsv_loader() - Term-Seq
nreads <- hash::hash(
  h0_r1 = 12663175,
  h0_r2 = 12837536,
  h0_r3 = 13049714,
  h3_r1 = 12126378, 
  h3_r2 = 14140123,
  h3_r3 = 12193795,
  h45_r1 = 10700392,
  h45_r2 = 11200870,
  h6_r1 = 10977966, 
  h6_r2 = 11532819,
  rna_h0_r1 = 10436615,
  rna_h0_r2 = 10319815,
  rna_h3_r1 = 9694862,
  rna_h3_r2 = 9728650,
  rna_h45_r1 = 9410193,
  rna_h45_r2 = 10430338,
  rna_h6_r1 = 8428220,
  rna_h6_r2 = 8825242
)

# region limits for calc of RT
#
# • Will calc bu dividing min1-to-max1 by min2-to-max2
#
RT_min1 <- -75
RT_max1 <- -50
RT_min2 <- 1
RT_max2 <- 25
nt_range <- 80 # should be larger both ways than min/max!

norm_limit = 2

norm_range <- hash::hash(
  min1 = nt_range+RT_min1+1,
  max1 = nt_range+RT_max1+1,
  min2 = nt_range+RT_min2+1,
  max2 = nt_range+RT_max2+1
)


##### Term-seq depths #####
# t = 0h
h0_r1_f_tsv <- get_tsv_input("term_h0_r1_f","h0_r1")
h0_r1_r_tsv <- get_tsv_input("term_h0_r1_r","h0_r1")
h0_r2_f_tsv <- get_tsv_input("term_h0_r2_f","h0_r2")
h0_r2_r_tsv <- get_tsv_input("term_h0_r2_r","h0_r2")
h0_r3_f_tsv <- get_tsv_input("term_h0_r3_f","h0_r3")
h0_r3_r_tsv <- get_tsv_input("term_h0_r3_r","h0_r3")

h0_tsv_f <- left_join(left_join(h0_r1_f_tsv, h0_r2_f_tsv), h0_r3_f_tsv)
h0_tsv_f <- h0_tsv_f %>% 
  mutate(n = rowMeans(select(h0_tsv_f,2:3)))

h0_tsv_r <- left_join(left_join(h0_r1_r_tsv, h0_r2_r_tsv), h0_r3_r_tsv)
h0_tsv_r <- h0_tsv_r %>% 
  mutate(n = rowMeans(select(h0_tsv_r,2:3)))

# t = 3h
h3_r1_f_tsv <- get_tsv_input("term_h3_r1_f","h3_r1")
h3_r1_r_tsv <- get_tsv_input("term_h3_r1_r","h3_r1")
h3_r2_f_tsv <- get_tsv_input("term_h3_r2_f","h3_r2")
h3_r2_r_tsv <- get_tsv_input("term_h3_r2_r","h3_r2")
h3_r3_f_tsv <- get_tsv_input("term_h3_r3_f","h3_r3")
h3_r3_r_tsv <- get_tsv_input("term_h3_r3_r","h3_r3")

h3_tsv_f <- left_join(left_join(h3_r1_f_tsv, h3_r2_f_tsv), h3_r3_f_tsv)
h3_tsv_f <- h3_tsv_f %>% 
  mutate(n = rowMeans(select(h3_tsv_f,2:3)))

h3_tsv_r <- left_join(left_join(h3_r1_r_tsv, h3_r2_r_tsv), h3_r3_r_tsv)
h3_tsv_r <- h3_tsv_r %>% 
  mutate(n = rowMeans(select(h3_tsv_r,2:3)))

# t = 4.5h
h45_r1_f_tsv <- get_tsv_input("term_h45_r1_f","h45_r1")
h45_r1_r_tsv <- get_tsv_input("term_h45_r1_r","h45_r1")
h45_r2_f_tsv <- get_tsv_input("term_h45_r2_f","h45_r2")
h45_r2_r_tsv <- get_tsv_input("term_h45_r2_r","h45_r2")

h45_tsv_f <- left_join(h45_r1_f_tsv, h45_r2_f_tsv)
h45_tsv_f <- h45_tsv_f %>% 
  mutate(n = rowMeans(select(h45_tsv_f,2:3)))

h45_tsv_r <- left_join(h45_r1_r_tsv, h45_r2_r_tsv)
h45_tsv_r <- h45_tsv_r %>% 
  mutate(n = rowMeans(select(h45_tsv_r,2:3)))

# t = 6h
h6_r1_f_tsv <- get_tsv_input("term_h6_r1_f","h6_r1")
h6_r1_r_tsv <- get_tsv_input("term_h6_r1_r","h6_r1")
h6_r2_f_tsv <- get_tsv_input("term_h6_r2_f","h6_r2")
h6_r2_r_tsv <- get_tsv_input("term_h6_r2_r","h6_r2")

h6_tsv_f <- left_join(h6_r1_f_tsv, h6_r2_f_tsv)
h6_tsv_f <- h6_tsv_f %>% 
  mutate(n = rowMeans(select(h6_tsv_f,2:3)))

h6_tsv_r <- left_join(h6_r1_r_tsv, h6_r2_r_tsv)
h6_tsv_r <- h6_tsv_r %>% 
  mutate(n = rowMeans(select(h6_tsv_r,2:3)))

##### RNA-seq depths #####

# All .tsv in NUGA_Data, generated by samtools depth

# t = 0h
h0_r1_f_rna <- get_tsv_input("h0_r1_f","rna_h0_r1")
h0_r1_r_rna <- get_tsv_input("h0_r1_r","rna_h0_r1")
h0_r2_f_rna <- get_tsv_input("h0_r2_f","rna_h0_r2")
h0_r2_r_rna <- get_tsv_input("h0_r2_r","rna_h0_r2")

h0_rna_f <- left_join(h0_r1_f_rna, h0_r2_f_rna)
h0_rna_f <- h0_rna_f %>% 
  mutate(n = rowMeans(select(h0_rna_f,2:3)))

h0_rna_r <- left_join(h0_r1_r_rna, h0_r2_r_rna)
h0_rna_r <- h0_rna_r %>% 
  mutate(n = rowMeans(select(h0_rna_r,2:3)))

# t = 3h
h3_r1_f_rna <- get_tsv_input("h3_r1_f","rna_h3_r1")
h3_r1_r_rna <- get_tsv_input("h3_r1_r","rna_h3_r1")
h3_r2_f_rna <- get_tsv_input("h3_r2_f","rna_h3_r2")
h3_r2_r_rna <- get_tsv_input("h3_r2_r","rna_h3_r2")

h3_rna_f <- left_join(h3_r1_f_rna, h3_r2_f_rna)
h3_rna_f <- h3_rna_f %>% 
  mutate(n = rowMeans(select(h3_rna_f,2:3)))

h3_rna_r <- left_join(h3_r1_r_rna, h3_r2_r_rna)
h3_rna_r <- h3_rna_r %>% 
  mutate(n = rowMeans(select(h3_rna_r,2:3)))

# t = 4.5h
h45_r1_f_rna <- get_tsv_input("h45_r1_f","rna_h45_r1")
h45_r1_r_rna <- get_tsv_input("h45_r1_r","rna_h45_r1")
h45_r2_f_rna <- get_tsv_input("h45_r2_f","rna_h45_r2")
h45_r2_r_rna <- get_tsv_input("h45_r2_r","rna_h45_r2")

h45_rna_f <- left_join(h45_r1_f_rna, h45_r2_f_rna)
h45_rna_f <- h45_rna_f %>% 
  mutate(n = rowMeans(select(h45_rna_f,2:3)))

h45_rna_r <- left_join(h45_r1_r_rna, h45_r2_r_rna)
h45_rna_r <- h45_rna_r %>% 
  mutate(n = rowMeans(select(h45_rna_r,2:3)))

# t = 6h
h6_r1_f_rna <- get_tsv_input("h6_r1_f","rna_h6_r1")
h6_r1_r_rna <- get_tsv_input("h6_r1_r","rna_h6_r1")
h6_r2_f_rna <- get_tsv_input("h6_r2_f","rna_h6_r2")
h6_r2_r_rna <- get_tsv_input("h6_r2_r","rna_h6_r2")

h6_rna_f <- left_join(h6_r1_f_rna, h6_r2_f_rna)
h6_rna_f <- h6_rna_f %>% 
  mutate(n = rowMeans(select(h6_rna_f,2:3)))

h6_rna_r <- left_join(h6_r1_r_rna, h6_r2_r_rna)
h6_rna_r <- h6_rna_r %>% 
  mutate(n = rowMeans(select(h6_rna_r,2:3)))


#### 5-nt termseq depths ####

##### 0h #####

# F
j <-  character(length(nrow(h0_rna_f)))
j[1] <- 0
j[2] <- 0
j[nrow(h0_tsv_f)] <- 0
j[nrow(h0_tsv_f)-1] <- 0
for(row in 3:(nrow(h0_tsv_f)-2)){
  k <- (h0_tsv_f$n[row-2]+
        h0_tsv_f$n[row-1]+
        h0_tsv_f$n[row]+
        h0_tsv_f$n[row+1]+
        h0_tsv_f$n[row+2])/5
  j[row] <- k
}

h0_tsv_f <- h0_tsv_f %>%
  mutate(peak := !!j)

# R

j <-  character(length(nrow(h0_tsv_r)))
j[1] <- 0
j[2] <- 0
j[nrow(h0_tsv_r)] <- 0
j[nrow(h0_tsv_r)-1] <- 0
for(row in 3:(nrow(h0_tsv_r)-2)){
  k <- (h0_tsv_r$n[row-2]+
          h0_tsv_r$n[row-1]+
          h0_tsv_r$n[row]+
          h0_tsv_r$n[row+1]+
          h0_tsv_r$n[row+2])/5
  j[row] <- k
}

h0_tsv_r <- h0_tsv_r %>%
  mutate(peak := !!j)

##### 3h #####

# F
j <-  character(length(nrow(h3_tsv_f)))
j[1] <- 0
j[2] <- 0
j[nrow(h3_tsv_f)] <- 0
j[nrow(h3_tsv_f)-1] <- 0
for(row in 3:(nrow(h3_tsv_f)-2)){
  k <- (h3_tsv_f$n[row-2]+
          h3_tsv_f$n[row-1]+
          h3_tsv_f$n[row]+
          h3_tsv_f$n[row+1]+
          h3_tsv_f$n[row+2])/5
  j[row] <- k
}

h3_tsv_f <- h3_tsv_f %>%
  mutate(peak := !!j)

# R

j <-  character(length(nrow(h3_tsv_r)))
j[1] <- 0
j[2] <- 0
j[nrow(h3_tsv_r)] <- 0
j[nrow(h3_tsv_r)-1] <- 0
for(row in 3:(nrow(h3_tsv_r)-2)){
  k <- (h3_tsv_r$n[row-2]+
          h3_tsv_r$n[row-1]+
          h3_tsv_r$n[row]+
          h3_tsv_r$n[row+1]+
          h3_tsv_r$n[row+2])/5
  j[row] <- k
}

h3_tsv_r <- h3_tsv_r %>%
  mutate(peak := !!j)

##### 45h #####

# F
j <-  character(length(nrow(h45_tsv_f)))
j[1] <- 0
j[2] <- 0
j[nrow(h45_tsv_f)] <- 0
j[nrow(h45_tsv_f)-1] <- 0
for(row in 3:(nrow(h45_tsv_f)-2)){
  k <- (h45_tsv_f$n[row-2]+
          h45_tsv_f$n[row-1]+
          h45_tsv_f$n[row]+
          h45_tsv_f$n[row+1]+
          h45_tsv_f$n[row+2])/5
  j[row] <- k
}

h45_tsv_f <- h45_tsv_f %>%
  mutate(peak := !!j)

# R

j <-  character(length(nrow(h45_tsv_r)))
j[1] <- 0
j[2] <- 0
j[nrow(h45_tsv_r)] <- 0
j[nrow(h45_tsv_r)-1] <- 0
for(row in 3:(nrow(h45_tsv_r)-2)){
  k <- (h45_tsv_r$n[row-2]+
          h45_tsv_r$n[row-1]+
          h45_tsv_r$n[row]+
          h45_tsv_r$n[row+1]+
          h45_tsv_r$n[row+2])/5
  j[row] <- k
}

h45_tsv_r <- h45_tsv_r %>%
  mutate(peak := !!j)

##### 6h #####

# F
j <-  character(length(nrow(h6_tsv_f)))
j[1] <- 0
j[2] <- 0
j[nrow(h6_tsv_f)] <- 0
j[nrow(h6_tsv_f)-1] <- 0
for(row in 3:(nrow(h6_tsv_f)-2)){
  k <- (h6_tsv_f$n[row-2]+
          h6_tsv_f$n[row-1]+
          h6_tsv_f$n[row]+
          h6_tsv_f$n[row+1]+
          h6_tsv_f$n[row+2])/5
  j[row] <- k
}

h6_tsv_f <- h6_tsv_f %>%
  mutate(peak := !!j)

# R

j <-  character(length(nrow(h6_tsv_r)))
j[1] <- 0
j[2] <- 0
j[nrow(h6_tsv_r)] <- 0
j[nrow(h6_tsv_r)-1] <- 0
for(row in 3:(nrow(h6_tsv_r)-2)){
  k <- (h6_tsv_r$n[row-2]+
          h6_tsv_r$n[row-1]+
          h6_tsv_r$n[row]+
          h6_tsv_r$n[row+1]+
          h6_tsv_r$n[row+2])/5
  j[row] <- k
}

h6_tsv_r <- h6_tsv_r %>%
  mutate(peak := !!j)

##### add to TTS #####

# Plus
tts_plus <- tts_old %>%
  filter(strand == "+")

f0 <- h0_tsv_f$peak[match(x = tts_plus$TTS.position, table = h0_tsv_f$pos)]
f3 <- h3_tsv_f$peak[match(x = tts_plus$TTS.position, table = h3_tsv_f$pos)]
f45 <- h45_tsv_f$peak[match(x = tts_plus$TTS.position, table = h45_tsv_f$pos)]
f6 <- h6_tsv_f$peak[match(x = tts_plus$TTS.position, table = h6_tsv_f$pos)]

tts_plus <- tts_plus %>%
  add_column(peak_0 = !!f0) %>%
  add_column(peak_3 = !!f3) %>%
  add_column(peak_45 = !!f45) %>%
  add_column(peak_6 = !!f6)

tts_plus <- tts_plus %>%
  mutate(newTTS3 = as.numeric(peak_0)/as.numeric(peak_3),
         newTTS45 = as.numeric(peak_0)/as.numeric(peak_45),
         newTTS6 = as.numeric(peak_0)/as.numeric(peak_6))

# Minus

tts_minus <- tts_old %>%
  filter(strand == "-")

r0 <- h0_tsv_f$peak[match(x = tts_minus$TTS.position, table = h0_tsv_f$pos)]
r3 <- h3_tsv_f$peak[match(x = tts_minus$TTS.position, table = h3_tsv_f$pos)]
r45 <- h45_tsv_f$peak[match(x = tts_minus$TTS.position, table = h45_tsv_f$pos)]
r6 <- h6_tsv_f$peak[match(x = tts_minus$TTS.position, table = h6_tsv_f$pos)]

tts_minus <- tts_minus %>%
  add_column(peak_0 = !!r0) %>%
  add_column(peak_3 = !!r3) %>%
  add_column(peak_45 = !!r45) %>%
  add_column(peak_6 = !!r6)

tts_minus <- tts_minus %>%
  mutate(newTTS3 = as.numeric(peak_0)/as.numeric(peak_3),
         newTTS45 = as.numeric(peak_0)/as.numeric(peak_45),
         newTTS6 = as.numeric(peak_0)/as.numeric(peak_6))

# Merge

tts_all <- rbind(tts_plus,tts_minus) %>%
  arrange(GeneID)

##### TTS Score #####
tts_test <- tts_all %>%
  filter(TTS.position %in% rho$TTS.position)


#tts_test <- tts_test %>%
#  mutate(newTTS3 = log2(newTTS3),
#         newTTS45 = log2(newTTS45),
#         newTTS6 = log2(newTTS6))


old3 <- rho$Rho3[match(x = tts_test$TTS.position, table = rho$TTS.position)]
old45 <- rho$Rho45[match(x = tts_test$TTS.position, table = rho$TTS.position)]
old6 <- rho$Rho6[match(x = tts_test$TTS.position, table = rho$TTS.position)]

tts_final <- tts_test %>%
  add_column(oldTTS3 = !!old3) %>%
  add_column(oldTTS45 = !!old45) %>%
  add_column(oldTTS6 = !!old6)

tts_final <- tts_final %>%
  mutate(oldTTS3 = as.numeric(oldTTS3),
         oldTTS45 = as.numeric(oldTTS45),
         oldTTS6 = as.numeric(oldTTS6)) %>%
  filter(peak_0 !=0 & peak_3 != 0 & peak_45 !=0 & peak_6 !=0)

##### Plot corr #####

ggplot(tts_final) +
  geom_point(aes(x=oldTTS3, y=newTTS3)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,12), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old TTS-score") +
  ylab("new TTS-score") +
  ggtitle("3h TTS-scores")

ggplot(tts_final) +
  geom_point(aes(x=oldTTS45, y=newTTS45)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,8), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old TTS-score") +
  ylab("new TTS-score") +
  ggtitle("4.5h TTS-scores")

ggplot(tts_final) +
  geom_point(aes(x=oldTTS6, y=newTTS6)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,12), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,12), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old TTS-score") +
  ylab("new TTS-score") +
  ggtitle("6h TTS-scores")

#### RT-score calc ####

###### 0h ######

# Plus
tts_h0_f <- tts_plus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+i, table=h0_rna_f$pos)]
  h <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+neg_i, table=h0_rna_f$pos)]
  tts_h0_f <- tts_h0_f %>% add_column(!!k := g, !!m := h)
}

tts_h0_f <- tts_h0_f %>% add_column("0" = h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position, table=h0_rna_f$pos)])

# Minus 

tts_h0_r <- tts_minus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position+neg_i, table=h0_rna_r$pos)]
  h <- h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position+i, table=h0_rna_r$pos)]
  tts_h0_r <- tts_h0_r %>% add_column(!!k := g, !!m := h)
}

tts_h0_r <- tts_h0_r %>% add_column("0" = h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position, table=h0_rna_r$pos)])

# Combined
tts_h0 <- rbind(tts_h0_f,tts_h0_r)

tts_h0 <- tts_h0[,order(as.numeric(names(tts_h0)))] 
tts_h0$depth1 <- rowMeans(subset(tts_h0, select = c(norm_range$min1, norm_range$max1)), na.rm = T)
tts_h0$depth2 <- rowMeans(subset(tts_h0, select = c(norm_range$min2, norm_range$max2)), na.rm = T)
tts_h0 <- tts_h0 %>%
  mutate(depth1 = ifelse(depth1>norm_limit,depth1,norm_limit),
         depth2 = ifelse(depth2>norm_limit,depth2,norm_limit),
         RT0 = depth2/depth1) %>%
  select(GeneID,TTS.position,depth1,depth2,RT0)

###### 3h ######

# Plus
tts_h3_f <- tts_plus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position+i, table=h3_rna_f$pos)]
  h <- h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position+neg_i, table=h3_rna_f$pos)]
  tts_h3_f <- tts_h3_f %>% add_column(!!k := g, !!m := h)
}

tts_h3_f <- tts_h3_f %>% add_column("0" = h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position, table=h3_rna_f$pos)])

# Minus 

tts_h3_r <- tts_minus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position+neg_i, table=h3_rna_r$pos)]
  h <- h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position+i, table=h3_rna_r$pos)]
  tts_h3_r <- tts_h3_r %>% add_column(!!k := g, !!m := h)
}

tts_h3_r <- tts_h3_r %>% add_column("0" = h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position, table=h3_rna_r$pos)])

# Combined
tts_h3 <- rbind(tts_h3_f,tts_h3_r)

tts_h3 <- tts_h3[,order(as.numeric(names(tts_h3)))] 
tts_h3$depth1 <- rowMeans(subset(tts_h3, select = c(norm_range$min1, norm_range$max1)), na.rm = T)
tts_h3$depth2 <- rowMeans(subset(tts_h3, select = c(norm_range$min2, norm_range$max2)), na.rm = T)
tts_h3 <- tts_h3 %>%
  mutate(depth1 = ifelse(depth1>norm_limit,depth1,norm_limit),
         depth2 = ifelse(depth2>norm_limit,depth2,norm_limit),
         RT3 = depth2/depth1) %>%
  select(GeneID,TTS.position,depth1,depth2,RT3)


###### 45h ######

# Plus
tts_h45_f <- tts_plus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position+i, table=h45_rna_f$pos)]
  h <- h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position+neg_i, table=h45_rna_f$pos)]
  tts_h45_f <- tts_h45_f %>% add_column(!!k := g, !!m := h)
}

tts_h45_f <- tts_h45_f %>% add_column("0" = h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position, table=h45_rna_f$pos)])

# Minus 

tts_h45_r <- tts_minus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position+neg_i, table=h45_rna_r$pos)]
  h <- h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position+i, table=h45_rna_r$pos)]
  tts_h45_r <- tts_h45_r %>% add_column(!!k := g, !!m := h)
}

tts_h45_r <- tts_h45_r %>% add_column("0" = h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position, table=h45_rna_r$pos)])

# Combined
tts_h45 <- rbind(tts_h45_f,tts_h45_r)

tts_h45 <- tts_h45[,order(as.numeric(names(tts_h45)))] 
tts_h45$depth1 <- rowMeans(subset(tts_h45, select = c(norm_range$min1, norm_range$max1)), na.rm = T)
tts_h45$depth2 <- rowMeans(subset(tts_h45, select = c(norm_range$min2, norm_range$max2)), na.rm = T)
tts_h45 <- tts_h45 %>%
  mutate(depth1 = ifelse(depth1>norm_limit,depth1,norm_limit),
         depth2 = ifelse(depth2>norm_limit,depth2,norm_limit),
         RT45 = depth2/depth1) %>%
  select(GeneID,TTS.position,depth1,depth2,RT45)


###### 6h ######

# Plus
tts_h6_f <- tts_plus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position+i, table=h6_rna_f$pos)]
  h <- h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position+neg_i, table=h6_rna_f$pos)]
  tts_h6_f <- tts_h6_f %>% add_column(!!k := g, !!m := h)
}

tts_h6_f <- tts_h6_f %>% add_column("0" = h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position, table=h6_rna_f$pos)])

# Minus 

tts_h6_r <- tts_minus %>%
  select(GeneID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position+neg_i, table=h6_rna_r$pos)]
  h <- h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position+i, table=h6_rna_r$pos)]
  tts_h6_r <- tts_h6_r %>% add_column(!!k := g, !!m := h)
}

tts_h6_r <- tts_h6_r %>% add_column("0" = h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position, table=h6_rna_r$pos)])

# Combined
tts_h6 <- rbind(tts_h6_f,tts_h6_r)

tts_h6 <- tts_h6[,order(as.numeric(names(tts_h6)))] 
tts_h6$depth1 <- rowMeans(subset(tts_h6, select = c(norm_range$min1, norm_range$max1)), na.rm = T)
tts_h6$depth2 <- rowMeans(subset(tts_h6, select = c(norm_range$min2, norm_range$max2)), na.rm = T)
tts_h6 <- tts_h6 %>%
  mutate(depth1 = ifelse(depth1>norm_limit,depth1,norm_limit),
         depth2 = ifelse(depth2>norm_limit,depth2,norm_limit),
         RT6 = depth2/depth1) %>%
  select(GeneID,TTS.position,depth1,depth2,RT6)

##### Final RT Scores #####



RT_table <- tibble(GeneID = tts_all$GeneID,
                   RT0 = tts_h0$RT0,
                   RT3 = tts_h3$RT3,
                   RT45 = tts_h45$RT45,
                   RT6 = tts_h6$RT6,
                   depth0 = tts_h0$depth1,
                   depth3 = tts_h3$depth1,
                   depth45 = tts_h45$depth1,
                   depth6 = tts_h6$depth1) %>%
  mutate(newRT3 = RT3/RT0,
         newRT45 = RT45/RT0,
         newRT6 = RT6/RT0) %>%
  select(GeneID,depth0,depth3,depth45,depth6,newRT3,newRT45,newRT6)

tts_all <- left_join(tts_all, RT_table)

##### Add to RT #####

rt_test <- tts_all %>%
  filter(TTS.position %in% nt_100$Position)

old3 <- nt_100$RT0_3[match(x = rt_test$TTS.position, table = nt_100$Position)]
old45 <- nt_100$RT0_45[match(x = rt_test$TTS.position, table = nt_100$Position)]
old6 <- nt_100$RT0_6[match(x = rt_test$TTS.position, table = nt_100$Position)]

rt_final <- rt_test %>%
  add_column(oldRT3 = !!old3) %>%
  add_column(oldRT45 = !!old45) %>%
  add_column(oldRT6 = !!old6)

rt_final <- rt_final %>%
  mutate(oldRT3 = as.numeric(oldRT3),
         oldRT45 = as.numeric(oldRT45),
         oldRT6 = as.numeric(oldRT6)) 

##### Plot corr #####

ggplot(rt_final) +
  geom_point(aes(x=oldRT3, y=newRT3)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old RT-score") +
  ylab("new RT-score") +
  ggtitle("3h RT-scores")

ggplot(rt_final) +
  geom_point(aes(x=oldRT45, y=newRT45)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old RT-score") +
  ylab("new RT-score") +
  ggtitle("45h RT-scores")

ggplot(rt_final) +
  geom_point(aes(x=oldRT6, y=newRT6)) +
  theme_alex() +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,6), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("old RT-score") +
  ylab("new RT-score") +
  ggtitle("6h RT-scores")

#### Do RT and TTS correlate? ####

compare_test <- tts_all %>%
  filter(peak_0 !=0 & peak_3 != 0 & peak_45 !=0 & peak_6 !=0)
  
ggplot(compare_test) +
  geom_point(aes(x=newTTS3, y=newRT3, colour=log2(depth3/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("3h")
  
ggplot(compare_test) +
  geom_point(aes(x=newTTS45, y=newRT45, colour=log2(depth45/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("4.5h")

ggplot(compare_test) +
  geom_point(aes(x=newTTS6, y=newRT6, colour=log2(depth6/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,6), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("6h")

# RTscore vs expression

ggplot(compare_test) +
  geom_point(aes(x=log2(newTTS3), y=log2(depth3/depth0)), size = 0.8, alpha = 0.2) +
  theme_alex() +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,5), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("Expression ratio") +
  ggtitle("3h")

ggplot(compare_test) +
  geom_point(aes(x=log2(newTTS45), y=log2(depth45/depth0)), size = 0.8, alpha = 0.2) +
  theme_alex() +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,5), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("Expression ratio") +
  ggtitle("45h")

ggplot(compare_test) +
  geom_point(aes(x=log2(newTTS6), y=log2(depth6/depth0)), size = 0.8, alpha = 0.2) +
  theme_alex() +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-5,5), expand=c(0,0)) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("Expression ratio") +
  ggtitle("6h")
  
# Is normalising the TTS-score to the expression the answer?!?!?!?!
# (Hint, probably not...)

compare_norm <- compare_test %>%
  mutate(exp3 = depth3/depth0,
         exp45 = depth45/depth0,
         exp6 = depth6/depth0) %>%
  mutate(normTTS3 = newTTS3/exp3,
         normTTS45 = newTTS45/exp45,
         normTTS6 = newTTS6/exp6)

ggplot(compare_norm) +
  geom_point(aes(x=normTTS3, y=newRT3, colour=log2(depth3/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("3h")

ggplot(compare_norm) +
  geom_point(aes(x=normTTS45, y=newRT45, colour=log2(depth45/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("4.5h")

ggplot(compare_norm) +
  geom_point(aes(x=normTTS6, y=newRT6, colour=log2(depth6/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("new TTS-score") +
  ylab("new RT-score") +
  ggtitle("6h")


# Is term-seq peak height correlating with RNA-seq coverage??

compare_limit <- compare_test %>%
#  filter(peak_0 > 4 & peak_3 >4 & peak_45 >4 & peak_6 >4)
  filter(as.numeric(peak_0) > 4)

ggplot(compare_test) +
  geom_point(aes(x=log2(as.numeric(peak_0)), y=log2(depth0), colour=log2(depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("Term-seq peak") +
  ylab("RNAseq depth") +
  ggtitle("0h")

ggplot(compare_test) +
  geom_point(aes(x=log2(as.numeric(peak_3)), y=log2(depth3), colour=log2(depth3/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("Term-seq peak") +
  ylab("RNAseq depth") +
  ggtitle("3h")

ggplot(compare_test) +
  geom_point(aes(x=log2(as.numeric(peak_45)), y=log2(depth45), colour=log2(depth45/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("Term-seq peak") +
  ylab("RNAseq depth") +
  ggtitle("4.5h")

ggplot(compare_test) +
  geom_point(aes(x=log2(as.numeric(peak_6)), y=log2(depth6), colour=log2(depth6/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("Term-seq peak") +
  ylab("RNAseq depth") +
  ggtitle("6h")
  
# Scatterplots at timepoints, yay!
ggplot(compare_limit) +
  geom_point(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3)), colour=log2(depth3/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_smooth(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3))), method = "lm", se = FALSE, colour = "grey10") +
  xlab("term-peak0") +
  ylab("term-peak3") +
  ggtitle("3h")

ggplot(compare_limit ) +
  geom_point(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_45)), colour=log2(depth45/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,10), expand = c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_smooth(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3))), method = "lm", se = FALSE, colour = "grey10") +
  xlab("term-peak0") +
  ylab("term-peak45") +
  ggtitle("4.5h")

ggplot(compare_limit) +
  geom_point(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_6)), colour=log2(depth6/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed") +
  geom_smooth(aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3))), method = "lm", se = FALSE, colour = "grey10") +
  xlab("term-peak0") +
  ylab("term-peak6") +
  ggtitle("6h")

# Red or blue?

red3 <- compare_test %>% filter(log2(depth3/depth0) < -1)
red45 <- compare_test %>% filter(log2(depth45/depth0) < -1)
red6 <- compare_test %>% filter(log2(depth6/depth0) < -1)

blue3 <- compare_test %>% filter(log2(depth3/depth0) > 1)
blue45 <- compare_test %>% filter(log2(depth45/depth0) > 1)
blue6 <- compare_test %>% filter(log2(depth6/depth0) > 1)

ggplot() +
  geom_point(data = red3, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3)), colour=log2(depth3/depth0)), size = 0.8) +
  geom_point(data = blue3, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_3)), colour=log2(depth3/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("term-peak0") +
  ylab("term-peak3") +
  ggtitle("3h")

ggplot() +
  geom_point(data = red45, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_45)), colour=log2(depth45/depth0)), size = 0.8) +
  geom_point(data = blue45, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_45)), colour=log2(depth45/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("term-peak0") +
  ylab("term-peak45") +
  ggtitle("45h")

ggplot() +
  geom_point(data = red6, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_6)), colour=log2(depth6/depth0)), size = 0.8) +
  geom_point(data = blue6, aes(x=log2(as.numeric(peak_0)), y=log2(as.numeric(peak_6)), colour=log2(depth6/depth0)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("term-peak0") +
  ylab("term-peak6") +
  ggtitle("6h")

# In-compare

ggplot(compare_test) +
  geom_point(aes(x=as.numeric(peak_3), y=as.numeric(peak_45), colour=log2(depth45/depth3)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("term-peak3") +
  ylab("term-peak45") +
  ggtitle("3h vs 4.5h")

ggplot(compare_test) +
  geom_point(aes(x=as.numeric(peak_45), y=as.numeric(peak_6), colour=log2(depth6/depth45)), size = 0.8) +
  theme_heatmap2() +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_x_continuous(limits=c(0,3), expand=c(0,0)) +
  scale_colour_gradient2(low = colours$red, mid = "grey80", high = colours$azure, midpoint = 0) +
  geom_abline(slope=1, intercept = 0) +
  xlab("term-peak45") +
  ylab("term-peak6") +
  ggtitle("4.5h vs 6h")
  

