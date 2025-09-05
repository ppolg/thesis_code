###########################################################################

# Additional data analysis for D'Halluin et al., 2023

#   • Plot depth around TTS/PS/others
#     - Also do groupuing, heatmap

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

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

# geometric means, because R has no base function?!?!?
geomean <- function(x) exp(mean(log(x)))

# tsv loader
get_tsv_input <- function(tsv_file, norm_to) {
  input_tsv <- paste(here::here("NUGA_data/"), tsv_file, ".tsv", sep = "")
  column_name <- as.character(tsv_file)
  n <- nreads[[norm_to]]
  tsv <- read.csv(input_tsv, sep='\t') %>%
    dplyr::select(2:3) %>%
    dplyr::rename("pos"=1,!!column_name:=2)
  tsv[[column_name]] = round(((tsv[[column_name]]/n)*1000000),3) # to convert to CPM
  return(tsv)
}

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

# Get the dictionary with nreads
nreads <- hash::hash(
  h0_r1 = 10436615, h0_r1_f = 8672072, h0_r1_r = 1764543,
  h0_r2 = 10319815, h0_r2_f = 8692127, h0_r2_r = 1627688,
  h3_r1 = 9694862, h3_r1_f = 8159988, h3_r1_r = 1534874,
  h3_r2 = 9728650, h3_r2_f = 8195737, h3_r2_r = 1532913,
  h45_r1 = 9410193, h45_r1_f = 7904673, h45_r1_r = 1505520,
  h45_r2 = 10430338, h45_r2_f = 8835943, h45_r2_r = 1594395,
  h6_r1 = 8428220, h6_r1_f = 6923251, h6_r1_r = 1504969,
  h6_r2 = 8825242, h6_r2_f = 7054755, h6_r2_r = 1770487
)

# The range in which to check coverage
nt_range <- 100

filter_limit <- 3

norm_min <- -75
norm_max <- -50

norm_range <- hash::hash(
  nt1 = nt_range+norm_min+1,
  nt2 = nt_range+norm_max+1
)

#RhoTermPredict Scores

preds <- read.csv(here::here("NUGA_Data/TTS_find_table.csv")) %>%
  dplyr::rename("pred100" = "RhoTermPredict..100nt.",
                "pred150" = "RhoTermPredict..150nt.") %>%
  mutate(pred150 = ifelse(pred150 == "Found",T,F)) %>%
  filter(pred150 == T) %>%
  select(ID,Locus,Gene)

########
# Main #
########

#### Load tsv  ####

# All .tsv in NUGA_Data, generated by samtools depth

# t = 0h
h0_r1_f_tsv <- get_tsv_input("h0_r1_f","h0_r1")
h0_r1_r_tsv <- get_tsv_input("h0_r1_r","h0_r1")
h0_r2_f_tsv <- get_tsv_input("h0_r2_f","h0_r2")
h0_r2_r_tsv <- get_tsv_input("h0_r2_r","h0_r2")

h0_tsv_f <- left_join(h0_r1_f_tsv, h0_r2_f_tsv)
h0_tsv_f <- h0_tsv_f %>% 
  mutate(n = rowMeans(select(h0_tsv_f,2:3)))

h0_tsv_r <- left_join(h0_r1_r_tsv, h0_r2_r_tsv)
h0_tsv_r <- h0_tsv_r %>% 
  mutate(n = rowMeans(select(h0_tsv_r,2:3)))

# t = 3h
h3_r1_f_tsv <- get_tsv_input("h3_r1_f","h3_r1")
h3_r1_r_tsv <- get_tsv_input("h3_r1_r","h3_r1")
h3_r2_f_tsv <- get_tsv_input("h3_r2_f","h3_r2")
h3_r2_r_tsv <- get_tsv_input("h3_r2_r","h3_r2")

h3_tsv_f <- left_join(h3_r1_f_tsv, h3_r2_f_tsv)
h3_tsv_f <- h3_tsv_f %>% 
  mutate(n = rowMeans(select(h3_tsv_f,2:3)))

h3_tsv_r <- left_join(h3_r1_r_tsv, h3_r2_r_tsv)
h3_tsv_r <- h3_tsv_r %>% 
  mutate(n = rowMeans(select(h3_tsv_r,2:3)))

# t = 4.5h
h45_r1_f_tsv <- get_tsv_input("h45_r1_f","h45_r1")
h45_r1_r_tsv <- get_tsv_input("h45_r1_r","h45_r1")
h45_r2_f_tsv <- get_tsv_input("h45_r2_f","h45_r2")
h45_r2_r_tsv <- get_tsv_input("h45_r2_r","h45_r2")

h45_tsv_f <- left_join(h45_r1_f_tsv, h45_r2_f_tsv)
h45_tsv_f <- h45_tsv_f %>% 
  mutate(n = rowMeans(select(h45_tsv_f,2:3)))

h45_tsv_r <- left_join(h45_r1_r_tsv, h45_r2_r_tsv)
h45_tsv_r <- h45_tsv_r %>% 
  mutate(n = rowMeans(select(h45_tsv_r,2:3)))

# t = 6h
h6_r1_f_tsv <- get_tsv_input("h6_r1_f","h6_r1")
h6_r1_r_tsv <- get_tsv_input("h6_r1_r","h6_r1")
h6_r2_f_tsv <- get_tsv_input("h6_r2_f","h6_r2")
h6_r2_r_tsv <- get_tsv_input("h6_r2_r","h6_r2")

h6_tsv_f <- left_join(h6_r1_f_tsv, h6_r2_f_tsv)
h6_tsv_f <- h6_tsv_f %>% 
  mutate(n = rowMeans(select(h6_tsv_f,2:3)))

h6_tsv_r <- left_join(h6_r1_r_tsv, h6_r2_r_tsv)
h6_tsv_r <- h6_tsv_r %>% 
  mutate(n = rowMeans(select(h6_tsv_r,2:3)))

# average coverage?
avg_cov <- hash::hash(
  h0_f = as.numeric(mean(h0_tsv_f[,4],na.rm = T)),
  h0_r = as.numeric(mean(h0_tsv_r[,4],na.rm = T)),
  h3_f = as.numeric(mean(h3_tsv_f[,4],na.rm = T)),
  h3_r = as.numeric(mean(h3_tsv_r[,4],na.rm = T)),
  h45_f = as.numeric(mean(h45_tsv_f[,4],na.rm = T)),
  h45_r = as.numeric(mean(h45_tsv_r[,4],na.rm = T)),
  h6_f = as.numeric(mean(h6_tsv_f[,4],na.rm = T)),
  h6_r = as.numeric(mean(h6_tsv_r[,4],na.rm = T))
)

#h0_tsv_f <- h0_tsv_f %>% mutate(n = ifelse(n<avg_cov$h0_f, 1, n))
#h0_tsv_r <- h0_tsv_r %>% mutate(n = ifelse(n<avg_cov$h0_r, 1, n))
#h3_tsv_f <- h3_tsv_f %>% mutate(n = ifelse(n<avg_cov$h3_f, 1, n))
#h3_tsv_r <- h3_tsv_r %>% mutate(n = ifelse(n<avg_cov$h3_r, 1, n))
#h45_tsv_f <- h45_tsv_f %>% mutate(n = ifelse(n<avg_cov$h45_f, 1, n))
#h45_tsv_r <- h45_tsv_r %>% mutate(n = ifelse(n<avg_cov$h45_r, 1, n))
#h6_tsv_f <- h6_tsv_f %>% mutate(n = ifelse(n<avg_cov$h6_f, 1, n))
#h6_tsv_r <- h6_tsv_r %>% mutate(n = ifelse(n<avg_cov$h6_r, 1, n))

###### TTS ######
# Import TTS table from supplementary figure to make it easy for everyone :)
tts_table <- read.csv(here::here("NUGA_Data/TTS_find_table.csv")) %>%
  filter(Gene != "rrl")

# Filter out TTS with 2+ times more expression downstream(!)
tts_naughty <- read.csv(here::here("NUGA_Data/naughty_list.csv"))

tts_table <- tts_table %>%
  filter(ID %nin% tts_naughty$ID)

tts_plus <- tts_table %>% filter(strand == "+")
tts_minus <- tts_table %>% filter(strand == "-")



cond <- read.csv(here::here("NUGA_Data/condTTS_table.csv"))
tts_rho_3h <- read.csv(here::here("NUGA_Data/TTS_3h.csv"))
tts_rho_45h <- read.csv(here::here("NUGA_Data/TTS_45h.csv"))
tts_rho_6h <- read.csv(here::here("NUGA_Data/TTS_6h.csv"))

###### 0h ######

# Plus
tts_h0_f <- tts_plus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+i, table=h0_tsv_f$pos)]
  h <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+neg_i, table=h0_tsv_f$pos)]
  tts_h0_f <- tts_h0_f %>% add_column(!!k := g, !!m := h)
}

tts_h0_f <- tts_h0_f %>% add_column("0" = h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position, table=h0_tsv_f$pos)])

# Minus 

tts_h0_r <- tts_minus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position+neg_i, table=h0_tsv_r$pos)]
  h <- h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position+i, table=h0_tsv_r$pos)]
  tts_h0_r <- tts_h0_r %>% add_column(!!k := g, !!m := h)
}

tts_h0_r <- tts_h0_r %>% add_column("0" = h0_tsv_r$n[base::match(x=tts_h0_r$TTS.position, table=h0_tsv_r$pos)])

# Combined
tts_h0 <- rbind(tts_h0_f,tts_h0_r)
tts_h0_pred <- tts_h0 %>% filter(ID %in% preds$ID)
tts_h0_cond <- tts_h0 %>% filter(ID %in% cond$ID)

#All
tts_h0 <- tts_h0 %>%
  select(3:ncol(tts_h0))

tts_h0 <- tts_h0[,order(as.numeric(names(tts_h0)))] 
tts_h0$depth <- rowMeans(subset(tts_h0, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h0 <- tts_h0 %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h0_nums <- tts_h0 %>% colMeans(., na.rm = T)

h0_table <- tibble(pos = (0-(nt_range)):(nt_range), n0 = h0_nums)

# RhoTermPred
tts_h0_pred <- tts_h0_pred %>%
  select(3:ncol(tts_h0_pred))

tts_h0_pred <- tts_h0_pred[,order(as.numeric(names(tts_h0_pred)))] 
tts_h0_pred$depth <- rowMeans(subset(tts_h0_pred, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h0_pred <- tts_h0_pred %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h0_pred_nums <- tts_h0_pred %>% colMeans(., na.rm = T)

h0_pred_table <- tibble(pos = (0-(nt_range)):(nt_range), n0_pred = h0_pred_nums)

# Cond
tts_h0_cond <- tts_h0_cond %>%
  select(3:ncol(tts_h0_cond))

tts_h0_cond <- tts_h0_cond[,order(as.numeric(names(tts_h0_cond)))] 
tts_h0_cond$depth <- rowMeans(subset(tts_h0_cond, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h0_cond <- tts_h0_cond %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h0_cond_nums <- tts_h0_cond %>% colMeans(., na.rm = T)

h0_cond_table <- tibble(pos = (0-(nt_range)):(nt_range), n0_cond = h0_cond_nums)

##### 0h, long, to prove a point #####

# Plus
long_range <- 200

norm_long <- hash::hash(
  nt1 = long_range+norm_min+1,
  nt2 = long_range+norm_max+1
)

tts_long_f <- tts_plus %>%
  select(ID,TTS.position)

for(i in 1:long_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+i, table=h0_tsv_f$pos)]
  h <- h0_tsv_f$n[base::match(x=tts_h0_f$TTS.position+neg_i, table=h0_tsv_f$pos)]
  tts_long_f <- tts_long_f %>% add_column(!!k := g, !!m := h)
}

tts_long_f <- tts_long_f %>% add_column("0" = h0_tsv_f$n[base::match(x=tts_long_f$TTS.position, table=h0_tsv_f$pos)])

# Minus 

tts_long_r <- tts_minus %>%
  select(ID,TTS.position)

for(i in 1:long_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h0_tsv_r$n[base::match(x=tts_long_r$TTS.position+neg_i, table=h0_tsv_r$pos)]
  h <- h0_tsv_r$n[base::match(x=tts_long_r$TTS.position+i, table=h0_tsv_r$pos)]
  tts_long_r <- tts_long_r %>% add_column(!!k := g, !!m := h)
}

tts_long_r <- tts_long_r %>% add_column("0" = h0_tsv_r$n[base::match(x=tts_long_r$TTS.position, table=h0_tsv_r$pos)])

# Combined
tts_long <- rbind(tts_long_f,tts_long_r)

#All
tts_long <- tts_long %>%
  select(3:ncol(tts_long))

tts_long <- tts_long[,order(as.numeric(names(tts_long)))] 
tts_long$depth <- rowMeans(subset(tts_long, select = c(norm_long$nt1, norm_long$nt2)), na.rm = T)
tts_long <- tts_long %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*long_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*long_range))

long_nums <- tts_long %>% colMeans(., na.rm = T)

long_table <- tibble(pos = (0-(long_range)):(long_range), n0 = long_nums)


###### 3h ######
# Plus
tts_h3_f <- tts_plus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position+i, table=h3_tsv_f$pos)]
  h <- h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position+neg_i, table=h3_tsv_f$pos)]
  tts_h3_f <- tts_h3_f %>% add_column(!!k := g, !!m := h)
}

tts_h3_f <- tts_h3_f %>% add_column("0" = h3_tsv_f$n[base::match(x=tts_h3_f$TTS.position, table=h3_tsv_f$pos)])

# Minus 

tts_h3_r <- tts_minus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position+neg_i, table=h3_tsv_r$pos)]
  h <- h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position+i, table=h3_tsv_r$pos)]
  tts_h3_r <- tts_h3_r %>% add_column(!!k := g, !!m := h)
}

tts_h3_r <- tts_h3_r %>% add_column("0" = h3_tsv_r$n[base::match(x=tts_h3_r$TTS.position, table=h3_tsv_r$pos)])

# Combined
tts_h3 <- rbind(tts_h3_f,tts_h3_r)
tts_h3_fdr <- tts_h3 %>% filter(ID %in% tts_rho_3h$ID)
tts_h3_pred <- tts_h3 %>% filter(ID %in% preds$ID)
tts_h3_cond <- tts_h3 %>% filter(ID %in% cond$ID)


#All
tts_h3 <- tts_h3 %>%
  select(3:ncol(tts_h3))

tts_h3 <- tts_h3[,order(as.numeric(names(tts_h3)))]
tts_h3$depth <- rowMeans(subset(tts_h3, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)

tts_h3 <- tts_h3 %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h3_nums <- tts_h3 %>% colMeans(., na.rm = T)

h3_table <- tibble(pos = (0-(nt_range)):(nt_range), n3_pred = h3_nums)

# FDR

tts_h3_fdr <- tts_h3_fdr %>%
  select(3:ncol(tts_h3_fdr))

tts_h3_fdr <- tts_h3_fdr[,order(as.numeric(names(tts_h3_fdr)))]
tts_h3_fdr$depth <- rowMeans(subset(tts_h3_fdr, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)


tts_h3_fdr <- tts_h3_fdr %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h3_nums_fdr <- tts_h3_fdr %>% colMeans(., na.rm = T)

h3_table_fdr <- tibble(pos = (0-(nt_range)):(nt_range), n3_fdr = h3_nums_fdr)

# RhoTermPred
tts_h3_pred <- tts_h3_pred %>%
  select(3:ncol(tts_h3_pred))

tts_h3_pred <- tts_h3_pred[,order(as.numeric(names(tts_h3_pred)))] 
tts_h3_pred$depth <- rowMeans(subset(tts_h3_pred, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h3_pred <- tts_h3_pred %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h3_pred_nums <- tts_h3_pred %>% colMeans(., na.rm = T)

h3_pred_table <- tibble(pos = (0-(nt_range)):(nt_range), n3 = h3_pred_nums)

# Cond
tts_h3_cond <- tts_h3_cond %>%
  select(3:ncol(tts_h3_cond))

tts_h3_cond <- tts_h3_cond[,order(as.numeric(names(tts_h3_cond)))] 
tts_h3_cond$depth <- rowMeans(subset(tts_h3_cond, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h3_cond <- tts_h3_cond %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h3_cond_nums <- tts_h3_cond %>% colMeans(., na.rm = T)

h3_cond_table <- tibble(pos = (0-(nt_range)):(nt_range), n3_cond = h3_cond_nums)

###### 4.5h ######

# Plus
tts_h45_f <- tts_plus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position+i, table=h45_tsv_f$pos)]
  h <- h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position+neg_i, table=h45_tsv_f$pos)]
  tts_h45_f <- tts_h45_f %>% add_column(!!k := g, !!m := h)
}

tts_h45_f <- tts_h45_f %>% add_column("0" = h45_tsv_f$n[base::match(x=tts_h45_f$TTS.position, table=h45_tsv_f$pos)])

# Minus 

tts_h45_r <- tts_minus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position+neg_i, table=h45_tsv_r$pos)]
  h <- h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position+i, table=h45_tsv_r$pos)]
  tts_h45_r <- tts_h45_r %>% add_column(!!k := g, !!m := h)
}

tts_h45_r <- tts_h45_r %>% add_column("0" = h45_tsv_r$n[base::match(x=tts_h45_r$TTS.position, table=h45_tsv_r$pos)])

# Combined
tts_h45 <- rbind(tts_h45_f,tts_h45_r)
tts_h45_fdr <- tts_h45 %>% filter(ID %in% tts_rho_45h$ID)
tts_h45_pred <- tts_h45 %>% filter(ID %in% preds$ID)
tts_h45_cond <- tts_h45 %>% filter(ID %in% cond$ID)

#All
tts_h45 <- tts_h45 %>%
  select(3:ncol(tts_h45))

tts_h45 <- tts_h45[,order(as.numeric(names(tts_h45)))]
tts_h45$depth <- rowMeans(subset(tts_h45, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)

tts_h45 <- tts_h45 %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h45_nums <- tts_h45 %>% colMeans(., na.rm = T)

h45_table <- tibble(pos = (0-(nt_range)):(nt_range), n45 = h45_nums)

# FDR
tts_h45_fdr <- tts_h45_fdr %>%
  select(3:ncol(tts_h45_fdr))

tts_h45_fdr <- tts_h45_fdr[,order(as.numeric(names(tts_h45_fdr)))]
tts_h45_fdr$depth <- rowMeans(subset(tts_h45_fdr, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)

tts_h45_fdr <- tts_h45_fdr %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h45_nums_fdr <- tts_h45_fdr %>% colMeans(., na.rm = T)

h45_table_fdr <- tibble(pos = (0-(nt_range)):(nt_range), n45_fdr = h45_nums_fdr)

# RhoTermPred
tts_h45_pred <- tts_h45_pred %>%
  select(3:ncol(tts_h45_pred))

tts_h45_pred <- tts_h45_pred[,order(as.numeric(names(tts_h45_pred)))] 
tts_h45_pred$depth <- rowMeans(subset(tts_h45_pred, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h45_pred <- tts_h45_pred %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h45_pred_nums <- tts_h45_pred %>% colMeans(., na.rm = T)

h45_pred_table <- tibble(pos = (0-(nt_range)):(nt_range), n45_pred = h45_pred_nums)

# Cond
tts_h45_cond <- tts_h45_cond %>%
  select(3:ncol(tts_h45_cond))

tts_h45_cond <- tts_h45_cond[,order(as.numeric(names(tts_h45_cond)))] 
tts_h45_cond$depth <- rowMeans(subset(tts_h45_cond, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h45_cond <- tts_h45_cond %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h45_cond_nums <- tts_h45_cond %>% colMeans(., na.rm = T)

h45_cond_table <- tibble(pos = (0-(nt_range)):(nt_range), n45_cond = h45_cond_nums)

###### 6h ######

# Plus
tts_h6_f <- tts_plus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position+i, table=h6_tsv_f$pos)]
  h <- h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position+neg_i, table=h6_tsv_f$pos)]
  tts_h6_f <- tts_h6_f %>% add_column(!!k := g, !!m := h)
}

tts_h6_f <- tts_h6_f %>% add_column("0" = h6_tsv_f$n[base::match(x=tts_h6_f$TTS.position, table=h6_tsv_f$pos)])

# Minus 

tts_h6_r <- tts_minus %>%
  select(ID,TTS.position)

for(i in 1:nt_range){
  neg_i <- 0-i
  k <- paste0(i)
  m <- paste0(neg_i)
  g <- h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position+neg_i, table=h6_tsv_r$pos)]
  h <- h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position+i, table=h6_tsv_r$pos)]
  tts_h6_r <- tts_h6_r %>% add_column(!!k := g, !!m := h)
}

tts_h6_r <- tts_h6_r %>% add_column("0" = h6_tsv_r$n[base::match(x=tts_h6_r$TTS.position, table=h6_tsv_r$pos)])

# Combined
tts_h6 <- rbind(tts_h6_f,tts_h6_r)
tts_h6_fdr <- tts_h6 %>% filter(ID %in% tts_rho_6h$ID)
tts_h6_pred <- tts_h6 %>% filter(ID %in% preds$ID)
tts_h6_cond <- tts_h6 %>% filter(ID %in% cond$ID)

#All
tts_h6 <- tts_h6 %>%
  select(3:ncol(tts_h6))

tts_h6 <- tts_h6[,order(as.numeric(names(tts_h6)))]
tts_h6$depth <- rowMeans(subset(tts_h6, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)

tts_h6 <- tts_h6 %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h6_nums <- tts_h6 %>% colMeans(., na.rm = T)

h6_table <- tibble(pos = (0-(nt_range)):(nt_range), n6 = h6_nums)

# FDR
tts_h6_fdr <- tts_h6_fdr %>%
  select(3:ncol(tts_h6_fdr))

tts_h6_fdr <- tts_h6_fdr[,order(as.numeric(names(tts_h6_fdr)))]
tts_h6_fdr$depth <- rowMeans(subset(tts_h6_fdr, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)

tts_h6_fdr <- tts_h6_fdr %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h6_nums_fdr <- tts_h6_fdr %>% colMeans(., na.rm = T)

h6_table_fdr <- tibble(pos = (0-(nt_range)):(nt_range), n6_fdr = h6_nums_fdr)

# RhoTermPred
tts_h6_pred <- tts_h6_pred %>%
  select(3:ncol(tts_h6_pred))

tts_h6_pred <- tts_h6_pred[,order(as.numeric(names(tts_h6_pred)))] 
tts_h6_pred$depth <- rowMeans(subset(tts_h6_pred, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h6_pred <- tts_h6_pred %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h6_pred_nums <- tts_h6_pred %>% colMeans(., na.rm = T)

h6_pred_table <- tibble(pos = (0-(nt_range)):(nt_range), n6_pred = h6_pred_nums)

# Cond
tts_h6_cond <- tts_h6_cond %>%
  select(3:ncol(tts_h6_cond))

tts_h6_cond <- tts_h6_cond[,order(as.numeric(names(tts_h6_cond)))] 
tts_h6_cond$depth <- rowMeans(subset(tts_h6_cond, select = c(norm_range$nt1, norm_range$nt2)), na.rm = T)
tts_h6_cond <- tts_h6_cond %>%
  mutate(test = ifelse(depth>filter_limit,"Y","N")) %>%
  filter(test == "Y") %>%
  select(1:(2+2*nt_range)) %>%
  mutate(across(.fns=~./`depth`)) %>%
  select(1:(1+2*nt_range))

h6_cond_nums <- tts_h6_cond %>% colMeans(., na.rm = T)

h6_cond_table <- tibble(pos = (0-(nt_range)):(nt_range), n6_cond = h6_cond_nums)

##### Combine/melt #####
all_term_table <- left_join(h0_table,h3_table) %>%
  left_join(.,h45_table) %>%
  left_join(.,h6_table) %>%
  reshape2::melt(id.vars="pos")

fdr_term_table <- left_join(h0_table,h3_table_fdr) %>%
  left_join(.,h45_table_fdr) %>%
  left_join(.,h6_table_fdr) %>%
  reshape2::melt(id.vars="pos")

pred_term_table <- left_join(h0_pred_table,h3_pred_table) %>%
  left_join(.,h45_pred_table) %>%
  left_join(.,h6_pred_table) %>%
  reshape2::melt(id.vars="pos")

cond_term_table <- left_join(h0_cond_table,h3_cond_table) %>%
  left_join(.,h45_cond_table) %>%
  left_join(.,h6_cond_table) %>%
  reshape2::melt(id.vars="pos")

long_term_table <- reshape2::melt(long_table)

########
# Plot #
########

#All
ggplot(data=all_term_table,aes(x=pos,y=log(value,2),colour=variable)) +
  geom_line(size = 1.8, alpha = 0.8) +
  theme_alex() +
  xlab("Relative position to TTS") +
  ylab("Change in coverage (log2)") +
  scale_x_continuous(limits = c(-100,100),expand=c(0,0)) +
  scale_y_continuous(limits = c(-1.25,1.25),expand=c(0,0), breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1")) +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green)) +
  theme(legend.position = "none") +
  geom_vline(size = 0.6, xintercept = 0, alpha = 0.6, linetype = "longdash") +
  ggtitle("Coverage per nt, normalised to expression")

#FDR
ggplot(data=fdr_term_table,aes(x=pos,y=log(value,2),colour=variable)) +
  geom_line(size = 1.8, alpha = 0.8) +
  theme_alex() +
  xlab("Relative position to TTS") +
  ylab("Change in coverage (log2)") +
  scale_x_continuous(limits = c(-100,100),expand=c(0,0)) +
  scale_y_continuous(limits = c(-1.1,1),expand=c(0,0), breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1")) +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green)) +
  theme(legend.position = "none") +
  geom_vline(size = 0.6, xintercept = 0, alpha = 0.6, linetype = "longdash") +
  ggtitle("Coverage per nt, for FDR-based Rho-TTS")

#Long
ggplot(data=long_table,aes(x=pos,y=log(n0,2))) +
  geom_line(size = 1.8, alpha = 1, colour = "grey40") +
  theme_alex() +
  xlab("Relative position to TTS") +
  ylab("Change in coverage (log2)") +
  scale_x_continuous(limits = c(-100,200),expand=c(0,0)) +
  scale_y_continuous(limits = c(-1.5,1.5),expand=c(0,0), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("-1.5","-1","-0.5","0","0.5","1","1.5")) +
  theme(legend.position = "none") +
  geom_vline(size = 0.6, xintercept = 0, alpha = 0.6, linetype = "longdash") +
  ggtitle("Coverage per nt, long (for showing 100nt limit)")

#RTP
ggplot(data=pred_term_table,aes(x=pos,y=log(value,2),colour=variable)) +
  geom_line(size = 1.8, alpha = 0.8) +
  theme_alex() +
  xlab("Relative position to TTS") +
  ylab("Change in coverage (log2)") +
  scale_x_continuous(limits = c(-100,100),expand=c(0,0)) +
  scale_y_continuous(limits = c(-1.25,1.25),expand=c(0,0), breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1")) +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green)) +
  theme(legend.position = "none") +
  geom_vline(size = 0.6, xintercept = 0, alpha = 0.6, linetype = "longdash") +
  ggtitle("Coverage per nt, for RhoTermPredicted")

#Cond
ggplot(data=cond_term_table,aes(x=pos,y=log(value,2),colour=variable)) +
  geom_line(size = 1.8, alpha = 0.8) +
  theme_alex() +
  xlab("Relative position to TTS") +
  ylab("Change in coverage (log2)") +
  scale_x_continuous(limits = c(-100,100),expand=c(0,0)) +
  scale_y_continuous(limits = c(-1.25,1.25),expand=c(0,0), breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1")) +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green)) +
  theme(legend.position = "none") +
  geom_vline(size = 0.6, xintercept = 0, alpha = 0.6, linetype = "longdash") +
  ggtitle("Coverage per nt, for CondTTS")


#### Heatmap ####
zmin <- 0.25
zmax <- 8

###### tts_h0 ######
tts_h0_melt <- tts_h0 
tts_h0_melt[is.na(tts_h0_melt)] = 0
h0_melt <- tts_h0_melt %>%
  mutate(k = c(1:nrow(tts_h0_melt)))

h0_melt$depth <- rowMeans(subset(h0_melt, select = c(126, 176)), na.rm = T)
h0_melt$k <- reorder(h0_melt$k,h0_melt$depth)
h0_melt <- h0_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h0_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$darkgold) +
  ggtitle("0h") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

###### tts_h3 ######
tts_h3_melt <- tts_h3
tts_h3_melt[is.na(tts_h3_melt)] = 0
h3_melt <- tts_h3_melt %>%
  mutate(k = c(1:nrow(tts_h3_melt)))

h3_melt$depth <- rowMeans(subset(h3_melt, select = c(126, 176)), na.rm = T)
h3_melt$k <- reorder(h3_melt$k,h3_melt$depth)
h3_melt <- h3_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h3_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$azure) +
  ggtitle("3h") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

####### tts_h45 ######
tts_h45_melt <- tts_h45
tts_h45_melt[is.na(tts_h45_melt)] = 0
h45_melt <- tts_h45_melt %>%
  mutate(k = c(1:nrow(tts_h45_melt)))

h45_melt$depth <- rowMeans(subset(h45_melt, select = c(126, 176)), na.rm = T)
h45_melt$k <- reorder(h45_melt$k,h45_melt$depth)
h45_melt <- h45_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h45_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$purple) +
  ggtitle("4.5h") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

###### tts_h6 ######
tts_h6_melt <- tts_h6
tts_h6_melt[is.na(tts_h6_melt)] = 0
h6_melt <- tts_h6_melt %>%
  mutate(k = c(1:nrow(tts_h6_melt)))

h6_melt$depth <- rowMeans(subset(h6_melt, select = c(126, 176)), na.rm = T)
h6_melt$k <- reorder(h6_melt$k,h6_melt$depth)
h6_melt <- h6_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h6_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("6h") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### FDR #####

# These will probably look mediocre with a low sample size?

###### 3h ######

tts_h3_fdr_melt <- tts_h3_fdr
tts_h3_fdr_melt[is.na(tts_h3_fdr_melt)] = 0
h3_fdr_melt <- tts_h3_fdr_melt %>%
  mutate(k = c(1:nrow(tts_h3_fdr_melt)))

h3_fdr_melt$depth <- rowMeans(subset(h3_fdr_melt, select = c(126, 176)), na.rm = T)
h3_fdr_melt$k <- reorder(h3_fdr_melt$k,h3_fdr_melt$depth)
h3_fdr_melt <- h3_fdr_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h3_fdr_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("3h, FDR") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Relative position to TTS (nt)")

###### 4.5h ######

tts_h45_fdr_melt <- tts_h45_fdr
tts_h45_fdr_melt[is.na(tts_h45_fdr_melt)] = 0
h45_fdr_melt <- tts_h45_fdr_melt %>%
  mutate(k = c(1:nrow(tts_h45_fdr_melt)))

h45_fdr_melt$depth <- rowMeans(subset(h45_fdr_melt, select = c(126, 176)), na.rm = T)
h45_fdr_melt$k <- reorder(h45_fdr_melt$k,h45_fdr_melt$depth)
h45_fdr_melt <- h45_fdr_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h45_fdr_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("4.5h, FDR") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Relative position to TTS (nt)")

###### 6h ######

tts_h6_fdr_melt <- tts_h6_fdr
tts_h6_fdr_melt[is.na(tts_h6_fdr_melt)] = 0
h6_fdr_melt <- tts_h6_fdr_melt %>%
  mutate(k = c(1:nrow(tts_h6_fdr_melt)))

h6_fdr_melt$depth <- rowMeans(subset(h6_fdr_melt, select = c(126, 176)), na.rm = T)
h6_fdr_melt$k <- reorder(h6_fdr_melt$k,h6_fdr_melt$depth)
h6_fdr_melt <- h6_fdr_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h6_fdr_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("6h, FDR") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Relative position to TTS (nt)")

#### RhoTermPredict heatmaps ####

##### 0h #####

tts_h0_pred_melt <- tts_h0_pred
tts_h0_pred_melt[is.na(tts_h0_pred_melt)] = 0
h0_pred_melt <- tts_h0_pred_melt %>%
  mutate(k = c(1:nrow(tts_h0_pred_melt)))

h0_pred_melt$depth <- rowMeans(subset(h0_pred_melt, select = c(126, 176)), na.rm = T)
h0_pred_melt$k <- reorder(h0_pred_melt$k,h0_pred_melt$depth)
h0_pred_melt <- h0_pred_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h0_pred_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$darkgold) +
  ggtitle("0h for RhoTermPredict") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 3h #####

tts_h3_pred_melt <- tts_h3_pred
tts_h3_pred_melt[is.na(tts_h3_pred_melt)] = 0
h3_pred_melt <- tts_h3_pred_melt %>%
  mutate(k = c(1:nrow(tts_h3_pred_melt)))

h3_pred_melt$depth <- rowMeans(subset(h3_pred_melt, select = c(126, 176)), na.rm = T)
h3_pred_melt$k <- reorder(h3_pred_melt$k,h3_pred_melt$depth)
h3_pred_melt <- h3_pred_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h3_pred_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$azure) +
  ggtitle("3h for RhoTermPredict") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 45h #####

tts_h45_pred_melt <- tts_h45_pred
tts_h45_pred_melt[is.na(tts_h45_pred_melt)] = 0
h45_pred_melt <- tts_h45_pred_melt %>%
  mutate(k = c(1:nrow(tts_h45_pred_melt)))

h45_pred_melt$depth <- rowMeans(subset(h45_pred_melt, select = c(126, 176)), na.rm = T)
h45_pred_melt$k <- reorder(h45_pred_melt$k,h45_pred_melt$depth)
h45_pred_melt <- h45_pred_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h45_pred_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$purple) +
  ggtitle("4.5h for RhoTermPredict") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 6h #####

tts_h6_pred_melt <- tts_h6_pred
tts_h6_pred_melt[is.na(tts_h6_pred_melt)] = 0
h6_pred_melt <- tts_h6_pred_melt %>%
  mutate(k = c(1:nrow(tts_h6_pred_melt)))

h6_pred_melt$depth <- rowMeans(subset(h6_pred_melt, select = c(126, 176)), na.rm = T)
h6_pred_melt$k <- reorder(h6_pred_melt$k,h6_pred_melt$depth)
h6_pred_melt <- h6_pred_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h6_pred_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("6h for RhoTermPredict") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")





##### CondTTS #####
##### 0h #####

tts_h0_cond_melt <- tts_h0_cond
tts_h0_cond_melt[is.na(tts_h0_cond_melt)] = 0
h0_cond_melt <- tts_h0_cond_melt %>%
  mutate(k = c(1:nrow(tts_h0_cond_melt)))

h0_cond_melt$depth <- rowMeans(subset(h0_cond_melt, select = c(126, 176)), na.rm = T)
h0_cond_melt$k <- reorder(h0_cond_melt$k,h0_cond_melt$depth)
h0_cond_melt <- h0_cond_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h0_cond_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$darkgold) +
  ggtitle("0h for CondTTS") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 3h #####

tts_h3_cond_melt <- tts_h3_cond
tts_h3_cond_melt[is.na(tts_h3_cond_melt)] = 0
h3_cond_melt <- tts_h3_cond_melt %>%
  mutate(k = c(1:nrow(tts_h3_cond_melt)))

h3_cond_melt$depth <- rowMeans(subset(h3_cond_melt, select = c(126, 176)), na.rm = T)
h3_cond_melt$k <- reorder(h3_cond_melt$k,h3_cond_melt$depth)
h3_cond_melt <- h3_cond_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h3_cond_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$azure) +
  ggtitle("3h for CondTTS") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 45h #####

tts_h45_cond_melt <- tts_h45_cond
tts_h45_cond_melt[is.na(tts_h45_cond_melt)] = 0
h45_cond_melt <- tts_h45_cond_melt %>%
  mutate(k = c(1:nrow(tts_h45_cond_melt)))

h45_cond_melt$depth <- rowMeans(subset(h45_cond_melt, select = c(126, 176)), na.rm = T)
h45_cond_melt$k <- reorder(h45_cond_melt$k,h45_cond_melt$depth)
h45_cond_melt <- h45_cond_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h45_cond_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$purple) +
  ggtitle("4.5h for CondTTS") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")

##### 6h #####

tts_h6_cond_melt <- tts_h6_cond
tts_h6_cond_melt[is.na(tts_h6_cond_melt)] = 0
h6_cond_melt <- tts_h6_cond_melt %>%
  mutate(k = c(1:nrow(tts_h6_cond_melt)))

h6_cond_melt$depth <- rowMeans(subset(h6_cond_melt, select = c(126, 176)), na.rm = T)
h6_cond_melt$k <- reorder(h6_cond_melt$k,h6_cond_melt$depth)
h6_cond_melt <- h6_cond_melt %>%
  select(-depth) %>%
  melt(id.vars="k",variable.name = "m") %>%
  mutate(k = as.numeric(k),
         m = as.numeric(m) - 101,
         value = ifelse(value > zmax,zmax,value),
         value = ifelse(value < zmin,zmin,value))


ggplot(h6_cond_melt, aes(x=m,y=k,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low="white",high=colours$green) +
  ggtitle("6h for CondTTS") +
  theme_heatmap() +
  scale_x_continuous(limits=c(-100,100), expand=c(0,0), labels = c(-100,50,"TTS",50,100)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend_heatmap()) +
  xlab("Relative position to TTS (nt)")




