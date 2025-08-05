###########################################################################

# Data analysis for D'Halluin et al., 2023
 
#   • Plot RT and Rho score distribution as requestsed
#     - Also do subsets of types, maybe RT in antisense is different?

###########################################################################

########
# Libs #
########

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats")
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
  alex_R3_edge = "#520101",
  alex_B1 = "#D7E0F7",
  alex_B2 = "#87A9D9",
  alex_B3 = "#5577AA",
  alex_B1_edge = "#BBD0FA"
)

# ggplot2 theme for figures
theme_alex <- function(base_size=20) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5, vjust=3),
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

# For filling 0's for better loess
k=round(seq.int(0,10, by=0.01),2) # Becasue seq.int without round gives weird results where n != n??
replacement_numbers <- tibble("RT" = k, n = 0 )
replacement_Rho <- tibble("Rho" = k, n = 0 )


#### RT/Rho ####


##### 100 nt (Fig 4C) #####

###### Full dataset ######
nt_100 <- read.csv(here::here("NUGA_Data/100nt_table.csv")) %>%
  dplyr::select(ID,Position,Strand,Locus,Class,RT.score.0H.3H,RT.score.0H.4.5H,RT.score.0H.6H) %>%
  dplyr::rename("RT0_3" = RT.score.0H.3H,
                "RT0_45" = RT.score.0H.4.5H,
                "RT0_6" = RT.score.0H.6H) %>%
  filter(!is.na(as.numeric(RT0_3)),
         !is.na(as.numeric(RT0_45)),
         !is.na(as.numeric(RT0_6))) %>%
  mutate(across(6:8, as.numeric),
         across(6:8, round, 2))

nt_100_sum3 <- nt_100 %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))


nt_100_3 <- rbind(nt_100_sum3,replacement_numbers %>% filter(RT %nin% nt_100_sum3$RT)) %>%
  arrange(RT)


nt_100_sum45 <- nt_100 %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

nt_100_45 <- rbind(nt_100_sum45,replacement_numbers %>% filter(RT %nin% nt_100_sum45$RT)) %>%
  arrange(RT)

nt_100_sum6 <- nt_100 %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

nt_100_6 <- rbind(nt_100_sum6,replacement_numbers %>% filter(RT %nin% nt_100_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = nt_100_3,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = nt_100_45,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B2", fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = nt_100_6,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B3", fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(-1,2), expand = c(0,0), breaks = c(-1,0,1,2)) +
  scale_y_continuous(limits = c(0,50), expand = c(0,0)) +
  geom_vline(xintercept = 0, alpha = 0.8) +
  xlab("RT Score (log2, all TTS)") +
  ylab("Number of TTS") +
  ggtitle("Fig 4C replacement")



###### Only I ######
nt_100_I <- nt_100 %>%
  dplyr::filter(Class == "I")

nt_100_I_sum3 <- nt_100_I %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))

nt_100_I_3 <- rbind(nt_100_I_sum3,replacement_numbers %>% filter(RT %nin% nt_100_I_sum3$RT)) %>%
  arrange(RT)

nt_100_I_sum45 <- nt_100_I %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

nt_100_I_45 <- rbind(nt_100_I_sum45,replacement_numbers %>% filter(RT %nin% nt_100_I_sum45$RT)) %>%
  arrange(RT)

nt_100_I_sum6 <- nt_100_I %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

nt_100_I_6 <- rbind(nt_100_I_sum6,replacement_numbers %>% filter(RT %nin% nt_100_I_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = nt_100_I_3,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = nt_100_I_45,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B2",fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = nt_100_I_6,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.05))),
              colour=colours$"alex_B3",fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,35), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("RT ratio") +
  ylab("Count") +
  ggtitle("Fig 4C - Internal")

###### Only O ######
nt_100_O <- nt_100 %>%
  dplyr::filter(Class == "O")

nt_100_O_sum3 <- nt_100_O %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))

nt_100_O_3 <- rbind(nt_100_O_sum3,replacement_numbers %>% filter(RT %nin% nt_100_O_sum3$RT)) %>%
  arrange(RT)

nt_100_O_sum45 <- nt_100_O %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

nt_100_O_45 <- rbind(nt_100_O_sum45,replacement_numbers %>% filter(RT %nin% nt_100_O_sum45$RT)) %>%
  arrange(RT)

nt_100_O_sum6 <- nt_100_O %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

nt_100_O_6 <- rbind(nt_100_O_sum6,replacement_numbers %>% filter(RT %nin% nt_100_O_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = nt_100_O_3,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = nt_100_O_45,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B2",fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = nt_100_O_6,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B3",fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("RT ratio") +
  ylab("Count") +
  ggtitle("Fig 4C - Orphan")

###### Only A ######
nt_100_A <- nt_100 %>%
  dplyr::filter(Class == "A")

nt_100_A_sum3 <- nt_100_A %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))

nt_100_A_3 <- rbind(nt_100_A_sum3,replacement_numbers %>% filter(RT %nin% nt_100_A_sum3$RT)) %>%
  arrange(RT)

nt_100_A_sum45 <- nt_100_A %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

nt_100_A_45 <- rbind(nt_100_A_sum45,replacement_numbers %>% filter(RT %nin% nt_100_A_sum45$RT)) %>%
  arrange(RT)

nt_100_A_sum6 <- nt_100_A %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

nt_100_A_6 <- rbind(nt_100_A_sum6,replacement_numbers %>% filter(RT %nin% nt_100_A_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = nt_100_A_3,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = nt_100_A_45,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B2",fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = nt_100_A_6,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B3",fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("RT ratio") +
  ylab("Count") +
  ggtitle("Fig 4C - Antisense")

###### Only P ######
nt_100_P <- nt_100 %>%
  dplyr::filter(Class == "P")

nt_100_P_sum3 <- nt_100_P %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))

nt_100_P_3 <- rbind(nt_100_P_sum3,replacement_numbers %>% filter(RT %nin% nt_100_P_sum3$RT)) %>%
  arrange(RT)

nt_100_P_sum45 <- nt_100_P %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

nt_100_P_45 <- rbind(nt_100_P_sum45,replacement_numbers %>% filter(RT %nin% nt_100_P_sum45$RT)) %>%
  arrange(RT)

nt_100_P_sum6 <- nt_100_P %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

nt_100_P_6 <- rbind(nt_100_P_sum6,replacement_numbers %>% filter(RT %nin% nt_100_P_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = nt_100_P_3,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = nt_100_P_45,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B2",fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = nt_100_P_6,
              aes(x=RT, y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.1))),
              colour=colours$"alex_B3",fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("RT ratio") +
  ylab("Count") +
  ggtitle("Fig 4C - Primary")



##### CondTTS (Fig 4F) #####

cond_tts <- read.csv(here::here("NUGA_Data/condTTS_table.csv")) %>%
  dplyr::select(ID,Strand,Locus,Class,RT.Score.0H.3H,RT.Score.0H.4.5H,RT.score.0H.6H) %>%
  dplyr::rename("RT0_3" = RT.Score.0H.3H,
                "RT0_45" = RT.Score.0H.4.5H,
                "RT0_6" = RT.score.0H.6H) %>%
  filter(!is.na(as.numeric(RT0_3)),
         !is.na(as.numeric(RT0_45)),
         !is.na(as.numeric(RT0_6))) %>%
  mutate(across(5:7, as.numeric),
         across(5:7, round, 2))

cond_tts_sum3 <- cond_tts %>%
  group_by(RT0_3) %>%
  dplyr::summarise(n = length(RT0_3)) %>%
  dplyr::rename("RT" = 'RT0_3') %>%
  dplyr::mutate(RT = as.numeric(RT))

cond_tts_3 <- rbind(cond_tts_sum3,replacement_numbers %>% filter(RT %nin% cond_tts_sum3$RT)) %>%
  arrange(RT)

cond_tts_sum45 <- cond_tts %>%
  group_by(RT0_45) %>%
  dplyr::summarise(n = length(RT0_45)) %>%
  dplyr::rename("RT" = 'RT0_45')

cond_tts_45 <- rbind(cond_tts_sum45,replacement_numbers %>% filter(RT %nin% cond_tts_sum45$RT)) %>%
  arrange(RT)

cond_tts_sum6 <- cond_tts %>%
  group_by(RT0_6) %>%
  dplyr::summarise(n = length(RT0_6)) %>%
  dplyr::rename("RT" = 'RT0_6')

cond_tts_6 <- rbind(cond_tts_sum6,replacement_numbers %>% filter(RT %nin% cond_tts_sum6$RT)) %>%
  arrange(RT)

ggplot() +
  theme_alex() +
  geom_ribbon(data = cond_tts_3,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.08))),
              colour=colours$"alex_B1_edge", fill=colours$"alex_B1", alpha=0.6) +
  geom_ribbon(data = cond_tts_45,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.08))),
              colour=colours$"alex_B2",fill=colours$"alex_B2", alpha=0.6) +
  geom_ribbon(data = cond_tts_6,
              aes(x=log2(RT), y=n, ymin=0, ymax=predict(loess(n ~ RT, span = 0.08))),
              colour=colours$"alex_B3",fill=colours$"alex_B3", alpha=0.6) +
  scale_x_continuous(limits = c(-1,2), expand = c(0,0), breaks = c(-1,0,1,2)) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  geom_vline(xintercept = 0, alpha = 0.8) +
  xlab("RT score (log2, CondTTS") +
  ylab("Number of TTS") +
  ggtitle("Fig 4F replacement")



##### Plot peak distribution (Fig 4B and 4E) #####

###### All - Fig 4B ######
rho <- read.csv(here::here("NUGA_Data/rho_table.csv")) %>%
  dplyr::mutate(LTF = ifelse(Leading == "Found", "Leading",
                             ifelse(Trailing == "Found", "Trailing", "Final"))) %>%
  dplyr::select(ID,TTS.position,Strand,Locus,Category,Class,LTF,Rho.Log2.0H.3H.,Rho.Log2.0H.4.5H.,Rho.Log2.0H.6H.) %>%
  dplyr::rename("Rho3" = Rho.Log2.0H.3H.,
                "Rho45" = Rho.Log2.0H.4.5H.,
                "Rho6" = Rho.Log2.0H.6H.) %>%
  filter(!is.na(as.numeric(Rho3)),
         !is.na(as.numeric(Rho45)),
         !is.na(as.numeric(Rho6))) %>%
  mutate(across(8:10, as.numeric),
         across(8:10, round, 2))

rho_sum3 <- rho %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_3 <- rbind(rho_sum3,replacement_Rho %>% filter(Rho %nin% rho_sum3$Rho)) %>%
  arrange(Rho)

rho_sum45 <- rho %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_45 <- rbind(rho_sum45,replacement_Rho %>% filter(Rho %nin% rho_sum45$Rho)) %>%
  arrange(Rho)

rho_sum6 <- rho %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_6 <- rbind(rho_sum6,replacement_Rho %>% filter(Rho %nin% rho_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_3,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.05))),
              colour=colours$"alex_R3_edge", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_45,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.07))),
              colour=colours$"alex_R3_edge",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_6,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.07))),
              colour=colours$"alex_R3_edge",fill=colours$"alex_R3",  alpha=0.6) +
  scale_x_continuous(limits = c(-1,2), expand = c(0,0), breaks = c(-2,-1,0,1,2)) +
  scale_y_continuous(limits = c(0,50), expand = c(0,0)) +
  geom_vline(xintercept = 0, alpha = 0.8) +
  xlab("TTS Score (log2, all TTS)") +
  ylab("Number of TTS") +
  ggtitle("Fig 4A")

  ###### Only I ######

rho_I <- rho %>%
  dplyr::filter(Class == "I")

rho_I_sum3 <- rho_I %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_I_3 <- rbind(rho_I_sum3,replacement_Rho %>% filter(Rho %nin% rho_I_sum3$Rho)) %>%
  arrange(Rho)

rho_I_sum45 <- rho_I %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_I_45 <- rbind(rho_I_sum45,replacement_Rho %>% filter(Rho %nin% rho_I_sum45$Rho)) %>%
  arrange(Rho)

rho_I_sum6 <- rho_I %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_I_6 <- rbind(rho_I_sum6,replacement_Rho %>% filter(Rho %nin% rho_I_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_I_3,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.05))),
              colour=colours$"alex_R1", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_I_45,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.05))),
              colour=colours$"alex_R2",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_I_6,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.05))),
              colour=colours$"alex_R3",fill=colours$"alex_R3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("log2 TTS Score") +
  ylab("Count") +
  ggtitle("Fig 4B - Internal")


###### Only O ######
rho_O <- rho %>%
  dplyr::filter(Class == "O")

rho_O_sum3 <- rho_O %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_O_3 <- rbind(rho_O_sum3,replacement_Rho %>% filter(Rho %nin% rho_O_sum3$Rho)) %>%
  arrange(Rho)

rho_O_sum45 <- rho_O %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_O_45 <- rbind(rho_O_sum45,replacement_Rho %>% filter(Rho %nin% rho_O_sum45$Rho)) %>%
  arrange(Rho)

rho_O_sum6 <- rho_O %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_O_6 <- rbind(rho_O_sum6,replacement_Rho %>% filter(Rho %nin% rho_O_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_O_3,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R1", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_O_45,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R2",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_O_6,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R3",fill=colours$"alex_R3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("log2 TTS Score") +
  ylab("Count") +
  ggtitle("Fig 4B - Orphan")



###### Only A ######
rho_A <- rho %>%
  dplyr::filter(Class == "A")

rho_A_sum3 <- rho_A %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_A_3 <- rbind(rho_A_sum3,replacement_Rho %>% filter(Rho %nin% rho_A_sum3$Rho)) %>%
  arrange(Rho)

rho_A_sum45 <- rho_A %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_A_45 <- rbind(rho_A_sum45,replacement_Rho %>% filter(Rho %nin% rho_A_sum45$Rho)) %>%
  arrange(Rho)

rho_A_sum6 <- rho_A %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_A_6 <- rbind(rho_A_sum6,replacement_Rho %>% filter(Rho %nin% rho_A_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_A_3,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R1", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_A_45,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R2",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_A_6,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R3",fill=colours$"alex_R3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("log2 TTS Score") +
  ylab("Count") +
  ggtitle("Fig 4B - Antisense")

###### Only P ######

rho_P <- rho %>%
  dplyr::filter(Class == "P")

rho_P_sum3 <- rho_P %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_P_3 <- rbind(rho_P_sum3,replacement_Rho %>% filter(Rho %nin% rho_P_sum3$Rho)) %>%
  arrange(Rho)

rho_P_sum45 <- rho_P %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_P_45 <- rbind(rho_P_sum45,replacement_Rho %>% filter(Rho %nin% rho_P_sum45$Rho)) %>%
  arrange(Rho)

rho_P_sum6 <- rho_P %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_P_6 <- rbind(rho_P_sum6,replacement_Rho %>% filter(Rho %nin% rho_P_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_P_3,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R1", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_P_45,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R2",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_P_6,
              aes(x=Rho, y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.1))),
              colour=colours$"alex_R3",fill=colours$"alex_R3", alpha=0.6) +
  scale_x_continuous(limits = c(0.25,2.5), expand = c(0,0), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  geom_vline(xintercept = 1, alpha = 0.8) +
  xlab("log2 TTS Score") +
  ylab("Count") +
  ggtitle("Fig 4B - Primary")

###### Conditional - Fig 4E ######
rho_cond <- rho %>%
  dplyr::filter(Category == "Cond_TTS")

rho_cond_sum3 <- rho_cond %>%
  group_by(Rho3) %>%
  dplyr::summarise(n = length(Rho3)) %>%
  dplyr::rename("Rho" = 'Rho3') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_cond_3 <- rbind(rho_cond_sum3,replacement_Rho %>% filter(Rho %nin% rho_cond_sum3$Rho)) %>%
  arrange(Rho)

rho_cond_sum45 <- rho_cond %>%
  group_by(Rho45) %>%
  dplyr::summarise(n = length(Rho45)) %>%
  dplyr::rename("Rho" = 'Rho45') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_cond_45 <- rbind(rho_cond_sum45,replacement_Rho %>% filter(Rho %nin% rho_cond_sum45$Rho)) %>%
  arrange(Rho)

rho_cond_sum6 <- rho_cond %>%
  group_by(Rho6) %>%
  dplyr::summarise(n = length(Rho6)) %>%
  dplyr::rename("Rho" = 'Rho6') %>%
  dplyr::mutate(Rho = as.numeric(Rho))

rho_cond_6 <- rbind(rho_cond_sum6,replacement_Rho %>% filter(Rho %nin% rho_cond_sum6$Rho)) %>%
  arrange(Rho)

ggplot() +
  theme_alex() +
  geom_ribbon(data = rho_cond_3,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.06))),
              colour=colours$"alex_R3_edge", fill=colours$"alex_R1", alpha=0.6) +
  geom_ribbon(data = rho_cond_45,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.08))),
              colour=colours$"alex_R3_edge",fill=colours$"alex_R2", alpha=0.6) +
  geom_ribbon(data = rho_cond_6,
              aes(x=log2(Rho), y=n, ymin=0, ymax=predict(loess(n ~ Rho, span = 0.08))),
              colour=colours$"alex_R3_edge",fill=colours$"alex_R3", alpha=0.6) +
  scale_x_continuous(limits = c(-1,2), expand = c(0,0), breaks = c(-1,0,1,2)) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  geom_vline(xintercept = 0, alpha = 0.8) +
  xlab("TTS Score (log2, CondTTS)") +
  ylab("Number of TTS") +
  ggtitle("Fig 4E")
