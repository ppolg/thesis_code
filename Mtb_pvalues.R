###########################################################################

# Additional data analysis for D'Halluin et al., 2022
#   • Note - most of this is copy-paste and not solidly defined functions,
#   • due to lack of time (and energy). Apologies in advance if you check it!

#   • p-value corrections for RT-score

###########################################################################

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","reshape2","DESeq2")
invisible(lapply(packages, require, character.only = TRUE))

theme_heatmap2 <- function(base_size=20) {
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
               label.theme = element_text(size = 20,
                                          family = "Arial",
                                          colour = "black",
                                          face = "bold"),
               ticks = F,
               nbin = 300,
               barwidth = 2,
               barheight = 14)
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

nt_100 <- read.csv(here::here("NUGA_Data/100nt_table.csv")) %>%
  dplyr::select(ID,Locus,Class,RT.score.0H.3H,RT.score.0H.4.5H,RT.score.0H.6H,Fisher.0.3H,Fisher.0.4.5H,Fisher.0.6H) %>%
  dplyr::rename(RT3 = RT.score.0H.3H,
                RT45 = RT.score.0H.4.5H,
                RT6 = RT.score.0H.6H,
                p3 = Fisher.0.3H,
                p45 = Fisher.0.4.5H,
                p6 = Fisher.0.6H)

nt_100 <- nt_100[-which(nt_100$ID == ""), ]

nt_100 <- nt_100 %>%
  dplyr::mutate(p3_FDR = p.adjust(p3, method = "BH"),
                p45_FDR = p.adjust(p45, method = "BH"),
                p6_FDR = p.adjust(p6, method = "BH"))



tts <- read.csv(here::here("NUGA_Data/tts_score.csv"))

rho <- read.csv(here::here("NUGA_Data/rho_table.csv")) %>%
  dplyr::mutate(LTF = ifelse(Leading == "Found", "Leading",
                             ifelse(Trailing == "Found", "Trailing", "Final"))) %>%
  dplyr::select(ID,TTS.position,Strand,Locus,Category,Class,LTF,Rho.Log2.0H.3H.,Rho.Log2.0H.4.5H.,Rho.Log2.0H.6H.) %>%
  dplyr::rename("Rho3" = Rho.Log2.0H.3H.,
                "Rho45" = Rho.Log2.0H.4.5H.,
                "Rho6" = Rho.Log2.0H.6H.) %>%
  mutate("p3" = tts$X3H.FDR,
         "p45" = tts$X4.5H.FDR,
         "p6" = tts$X6H.FDR) %>%
  filter(!is.na(as.numeric(Rho3)),
         !is.na(as.numeric(Rho45)),
         !is.na(as.numeric(Rho6)))




rt_3h <- nt_100 %>% dplyr::filter(RT3 > 1.1 & p3_FDR < 0.05)
rt_45h <- nt_100 %>% dplyr::filter(RT45 > 1.1 & p45_FDR < 0.05)
rt_6h <- nt_100 %>% dplyr::filter(RT6 > 1.1 & p6_FDR < 0.05)

rt_3h_tine <- nt_100 %>% dplyr::filter(RT3 > 1.1 & p3_FDR < 0.1)
rt_45h_tine <- nt_100 %>% dplyr::filter(RT45 > 1.1 & p45_FDR < 0.1)
rt_6h_tine <- nt_100 %>% dplyr::filter(RT6 > 1.1 & p6_FDR < 0.1)

tts_3h <- rho %>% dplyr::filter(Rho3 > 1.1 & p3 < 0.05)
tts_45h <- rho %>% dplyr::filter(Rho45 > 1.1 & p45 < 0.05)
tts_6h <- rho %>% dplyr::filter(Rho6 > 1.1 & p6 < 0.05)


pmin <- 0.005
pmax <- 0.2
pdiv <- 0.005

RTmin <- 1
RTmax <- 2
RTdiv <- 0.05

TTSmin <- 1
TTSmax <- 2
TTSdiv <- 0.05


p_values <- seq(pmin,pmax,by=pdiv)
RT_values <- seq(RTmin,RTmax,by=RTdiv)
TTS_values <- seq(TTSmin,TTSmax,by=TTSdiv)


#### 3h, no FDR ####
rt_3h_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT3 > n & p3 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_3h_table <- rt_3h_table %>%
    add_column(!!i := k)
}

h3_melt <- reshape2::melt(rt_3h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))
  
h3_melt <- h3_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)

#### 45h, no FDR ####
rt_45h_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT45 > n & p45 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_45h_table <- rt_45h_table %>%
    add_column(!!i := k)
}

h45_melt <- reshape2::melt(rt_45h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

h45_melt <- h45_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)

#### 6h, no FDR ####
rt_6h_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT6 > n & p6 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_6h_table <- rt_6h_table %>%
    add_column(!!i := k)
}

h6_melt <- reshape2::melt(rt_6h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

h6_melt <- h6_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)

#### 3h, FDR ####
rt_3h_fdr_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT3 > n & p3_FDR < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_3h_fdr_table <- rt_3h_fdr_table %>%
    add_column(!!i := k)
}

h3_fdr_melt <- reshape2::melt(rt_3h_fdr_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

h3_fdr_melt <- h3_fdr_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)

#### 45h, FDR ####
rt_45h_fdr_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT45 > n & p45_FDR < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_45h_fdr_table <- rt_45h_fdr_table %>%
    add_column(!!i := k)
}

h45_fdr_melt <- reshape2::melt(rt_45h_fdr_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

h45_fdr_melt <- h45_fdr_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)

#### 6h, FDR ####
rt_6h_fdr_table <- tibble(p = p_values)

for(RT in 1:length(RT_values)){
  workspace <- tibble()
  n <- RT_values[RT]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- nt_100 %>% dplyr::filter(RT6 > n & p6_FDR < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  rt_6h_fdr_table <- rt_6h_fdr_table %>%
    add_column(!!i := k)
}

h6_fdr_melt <- reshape2::melt(rt_6h_fdr_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

h6_fdr_melt <- h6_fdr_melt %>%  mutate(variable = (variable-1)*RTdiv + RTmin)


#### 3h, TTS score ####
tts_3h_table <- tibble(p = p_values)

for(TTS in 1:length(TTS_values)){
  workspace <- tibble()
  n <- TTS_values[TTS]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- rho %>% dplyr::filter(Rho3 > n & p3 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  tts_3h_table <- tts_3h_table %>%
    add_column(!!i := k)
}

tts_h3_melt <- reshape2::melt(tts_3h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

tts_h3_melt <- tts_h3_melt %>%  mutate(variable = (variable-1)*TTSdiv + TTSmin)

#### 45h, TTS score ####
tts_45h_table <- tibble(p = p_values)

for(TTS in 1:length(TTS_values)){
  workspace <- tibble()
  n <- TTS_values[TTS]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- rho %>% dplyr::filter(Rho45 > n & p45 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  tts_45h_table <- tts_45h_table %>%
    add_column(!!i := k)
}

tts_h45_melt <- reshape2::melt(tts_45h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

tts_h45_melt <- tts_h45_melt %>%  mutate(variable = (variable-1)*TTSdiv + TTSmin)

#### 6h, TTS score ####
tts_6h_table <- tibble(p = p_values)

for(TTS in 1:length(TTS_values)){
  workspace <- tibble()
  n <- TTS_values[TTS]
  i <- paste0(n)
  k <- c()
  for(p_val in 1:length(p_values)){
    l <- p_values[p_val]
    workspace <- rho %>% dplyr::filter(Rho6 > n & p6 < l)
    m <- nrow(workspace)
    k <- c(k,m)
  }
  tts_6h_table <- tts_6h_table %>%
    add_column(!!i := k)
}

tts_h6_melt <- reshape2::melt(tts_6h_table,id.vars="p") %>%
  mutate(variable = as.numeric(variable), value=as.numeric(value))

tts_h6_melt <- tts_h6_melt %>%  mutate(variable = (variable-1)*TTSdiv + TTSmin)

#### PLOT! ####

#RT - all
ggplot(h3_melt,aes(x=p, y=variable,fill=log(value,2 ))) +
  geom_tile() +
  scale_fill_gradient(low=colours$blue,high=colours$red,limits=c(NA,14)) +
  ggtitle("3h, no FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = rel(0.8)),
        legend.key.size= unit(0.6, "cm"),) +
  ylab("RT cutoff")

ggplot(h45_melt,aes(x=p, y=variable,fill=log(value,2))) +
  geom_tile() +
  scale_fill_gradient(low=colours$blue,high=colours$red, limits = c(NA,14)) +
  ggtitle("4.5h, no FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = rel(0.8)),
        legend.key.size= unit(0.6, "cm"),) +
  ylab("RT cutoff")

ggplot(h6_melt,aes(x=p, y=variable,fill=log2(value) )) +
  geom_tile() +
  scale_fill_gradient(low=colours$blue,high=colours$red, limits = c(8,14)) +
  ggtitle("6h, no FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = rel(0.8)),
        legend.key.size= unit(0.6, "cm"),) +
  ylab("RT cutoff")


#RT - FDR

ggplot(h3_fdr_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$alex_B1,high=colours$azure, limits=c(NA,110)) +
  ggtitle("3h, FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("RT cutoff")

ggplot(h45_fdr_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$alex_B1,high=colours$azure) +
  ggtitle("4.5h, FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("RT cutoff")

ggplot(h6_fdr_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$alex_B1,high=colours$azure) +
  ggtitle("6h, FDR") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("RT cutoff")

#TTSscore

ggplot(tts_h3_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$ashenskin,high=colours$alex_R3_edge, limits=c(0,1000)) +
  ggtitle("3h, TTS Scores") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("TTS cutoff")

ggplot(tts_h45_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$ashenskin,high=colours$alex_R3_edge, limits=c(0,1200)) +
  ggtitle("45h, TTS Scores") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("TTS cutoff")

ggplot(tts_h6_melt,aes(x=p, y=variable,fill=value)) +
  geom_tile() +
  scale_fill_gradient(low=colours$ashenskin,high=colours$alex_R3_edge, limits=c(600,1800)) +
  ggtitle("6h, TTS Scores") +
  theme_heatmap2() +
  scale_x_continuous(expand=c(0,0),labels = c("0","0.05","0.1","0.15","0.2")) +
  scale_y_continuous(expand=c(0,0),breaks=c(1,1.25,1.5,1.75,2),labels = c("1","1.25","1.5","1.75","2")) +
  guides(fill = guide_legend_heatmap()) +
  xlab("q-value cutoff") +
  ylab("TTS cutoff")

fwrite(rt_3h, file = here::here("NUGA_Out/RT_3h.csv"))
fwrite(rt_45h, file = here::here("NUGA_Out/RT_45h.csv"))
fwrite(rt_6h, file = here::here("NUGA_Out/RT_6h.csv"))

fwrite(tts_3h, file = here::here("NUGA_Out/TTS_3h.csv"))
fwrite(tts_45h, file = here::here("NUGA_Out/TTS_45h.csv"))
fwrite(tts_6h, file = here::here("NUGA_Out/TTS_6h.csv"))

  

