###########################################################################

# Additional data analysis for D'Halluin et al., 2023

#   • Plot distance between PS and TTS, as requested

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

# Function to read the bed(esque) files the TermSeq peak caller made
read_alex_bed <- function(file,strand = "+") {
  path <- paste(here::here("NUGA_data/"), file, ".bed", sep = "")
  bed <- as.tibble(read.table(path)) %>%
    select(2:3) %>%
    dplyr::rename(start = 1,
                  end = 2) %>%
    mutate(ID = row_number(),
           main = floor((start+end)/2)) %>%
    select(ID,main,start,end)
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
  blue = "#4A749E",
  purple = "#612882",
  alex_R1 = "#F0BFC1",
  alex_R2 = "#D46859",
  alex_R3 = "#871E12",
  alex_B1 = "#D7E0F7",
  alex_B2 = "#87A9D9",
  alex_B3 = "#5577AA",
  alex_B1_edge = "#BBD0FA",
  alex_grey = "#E5E5E5"
)

# ggplot2 theme for figures
theme_alex <- function(base_size=20) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5, vjust = 3),
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

########
# Main #
########

# Get TTS

nt_100 <- read.csv(here::here("NUGA_Data/100nt_table.csv")) %>%
  dplyr::select(ID,Position,Strand,Locus,Class) 

tts_plus <- nt_100 %>%
  dplyr::filter(Strand == "+") # Filter plus strand

tts_minus <- nt_100 %>%
  dplyr::filter(Strand != "+") # Filter minus strand

# Get PS

ps<- read.csv(here::here("NUGA_Data/PS_table.csv"))

ps_plus <- ps %>%
  dplyr::filter(Strand == "+") # Filter plus strand

ps_minus <- ps %>%
  dplyr::filter(Strand != "+") # Filter minus strand

# Get Termseq peaks

term_minus <- read_alex_bed("Peaks_ExpoMin")
term_plus <- read_alex_bed("Peaks_ExpoPlus")

term_minus_new <- read_alex_bed("PeaksCalling_min")
term_plus_new <- read_alex_bed("PeaksCalling_plus")

# Define range:
# This is the MAX range between TTS and PS we check for the curve.
nt_range <- 200

#### Old ####

# NB: This is now decrepit, only keeping it for completeness' sake

# Plus

#for(i in 0:nt_range){
#  k <- paste0(i)
#  term_plus <- term_plus %>% 
#  m <- ifelse((term_plus$main + i) %in% ps_plus$Position,1,0)
#    add_column(!!k := m)
#}


#for(i in 6:(5+nt_range)){
#  term_plus[,i] <- term_plus[,i] + term_plus[,i-1]
#}

#term_plus_final <- term_plus %>%
#  mutate(across (5:(5+nt_range), ~ ifelse(.x == 0,1,0)))

#term_plus_sum <- term_plus_final %>%
#  select(5:(5+nt_range)) %>%
#  summarise(across(where(is.numeric),sum))

#term_plus_melt <- reshape2::melt(term_plus_sum) %>%
#  mutate( var = as.numeric(variable)-1)

# Minus

#for(i in 0:nt_range){
#  k <- paste0(i)
#  m <- ifelse((term_minus$main - i) %in% ps_minus$Position,1,0)
#  term_minus <- term_minus %>% 
#    add_column(!!k := m)
#}


#for(i in 6:(5+nt_range)){
#  term_minus[,i] <- term_minus[,i] + term_minus[,i-1]
#}

#term_minus_final <- term_minus %>%
#  mutate(across (5:(5+nt_range), ~ ifelse(.x == 0,1,0)))

#term_minus_sum <- term_minus_final %>%
#  select(5:(5+nt_range)) %>%
#  summarise(across(where(is.numeric),sum))

#term_minus_melt <- reshape2::melt(term_minus_sum) %>%
#  mutate( var = as.numeric(variable)-1)

# Combine
#term_both_final <- rbind(term_plus_final,term_minus_final)

#term_both_sum <- term_both_final %>%
#  select(5:(5+nt_range)) %>%
#  summarise(across(where(is.numeric),sum))

#term_both_melt <- reshape2::melt(term_both_sum) %>%
#  mutate( var = as.numeric(variable)-1)

#### New ####
# Plus

for(i in 0:nt_range){
  k <- paste0(i)
  m <- ifelse((term_plus_new$main + i) %in% ps_plus$Position,1,0)
  term_plus_new <- term_plus_new %>% 
    add_column(!!k := m)
}


for(i in 6:(5+nt_range)){
  term_plus_new[,i] <- term_plus_new[,i] + term_plus_new[,i-1]
}

term_plus_new_final <- term_plus_new %>%
  mutate(across (5:(5+nt_range), ~ ifelse(.x == 0,1,0)))

term_plus_new_sum <- term_plus_new_final %>%
  select(5:(5+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

term_plus_new_melt <- reshape2::melt(term_plus_new_sum) %>%
  mutate( var = as.numeric(variable)-1)

# Minus

for(i in 0:nt_range){
  k <- paste0(i)
  m <- ifelse((term_minus_new$main - i) %in% ps_minus$Position,1,0)
  term_minus_new <- term_minus_new %>% 
    add_column(!!k := m)
}


for(i in 6:(5+nt_range)){
  term_minus_new[,i] <- term_minus_new[,i] + term_minus_new[,i-1]
}

term_minus_new_final <- term_minus_new %>%
  mutate(across (5:(5+nt_range), ~ ifelse(.x == 0,1,0)))

term_minus_new_sum <- term_minus_new_final %>%
  select(5:(5+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

term_minus_new_melt <- reshape2::melt(term_minus_new_sum) %>%
  mutate( var = as.numeric(variable)-1)

# Combine
term_both_new_final <- rbind(term_plus_new_final,term_minus_new_final)

term_both_new_sum <- term_both_new_final %>%
  select(5:(5+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

term_both_new_melt <- reshape2::melt(term_both_new_sum) %>%
  mutate( var = as.numeric(variable)-1)

#### Permutation test ####
test <- tibble(main = sample.int(4411532, size = nrow(term_both_new_final), replace = T),
               strand = ifelse(sample.int(2, size = nrow(term_both_new_final), replace = T) == 2,"+","-")) %>% 
  arrange(main)

test_plus <- test %>% filter(strand == "+")
test_minus <- test %>% filter(strand == "-")

# Plus

for(i in 0:nt_range){
  k <- paste0(i)
  m <- ifelse((test_plus$main + i) %in% ps_plus$Position,1,0)
  test_plus <- test_plus %>% 
    add_column(!!k := m)
}


for(i in 4:(3+nt_range)){
  test_plus[,i] <- test_plus[,i] + test_plus[,i-1]
}

test_plus_final <- test_plus %>%
  mutate(across (3:(3+nt_range), ~ ifelse(.x == 0,1,0)))

test_plus_sum <- test_plus_final %>%
  select(3:(3+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

test_plus_melt <- reshape2::melt(test_plus_sum) %>%
  mutate( var = as.numeric(variable)-1)

# Minus

for(i in 0:nt_range){
  k <- paste0(i)
  m <- ifelse((test_minus$main - i) %in% ps_minus$Position,1,0)
  test_minus <- test_minus %>% 
    add_column(!!k := m)
}

for(i in 4:(3+nt_range)){
  test_minus[,i] <- test_minus[,i] + test_minus[,i-1]
}

test_minus_final <- test_minus %>%
  mutate(across (3:(3+nt_range), ~ ifelse(.x == 0,1,0)))

test_minus_sum <- test_minus_final %>%
  select(3:(3+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

test_minus_melt <- reshape2::melt(test_minus_sum) %>%
  mutate( var = as.numeric(variable)-1)

# Combine
test_both_final <- rbind(test_plus_final,test_minus_final)

test_both_sum <- test_both_final %>%
  select(3:(3+nt_range)) %>%
  summarise(across(where(is.numeric),sum))

test_both_melt <- reshape2::melt(test_both_sum) %>%
  mutate( var = as.numeric(variable)-1)




#p <- wilcox.test(test_both_melt$value,term_both_melt$value, paired = F)$p.value
pnew <- wilcox.test(test_both_melt$value,term_both_new_melt$value, paired = F)$p.value

#### PLOT! ####
# This is the old one.
#ggplot() +
#  geom_line(data = term_both_melt,aes(x=var,y=value),size = 1.5, colour = colours$brown) +
#  geom_line(data = test_both_melt,aes(x=var,y=value),size = 1.5, colour = colours$grey) +
#  scale_y_log10(limits = c(1000,10000), breaks=c(1000,10000)) +
#  scale_y_continuous(expand=c(0,0),limits = c(0,7000)) +
#  scale_x_continuous(expand = c(0,0), breaks = c(0,50,100,150,200)) +
#  theme_alex() +
#  geom_vline(xintercept = 50, alpha = 0.8, linetype="dashed") +
#  xlab("Distance between Term peak and 5' mono-p") +
#  ylab("n of found TTS") +
#  ggtitle("Found TTS based on distance")


# Current Fig 1C
ggplot() +
  geom_line(data = term_both_new_melt,aes(x=var,y=value),size = 4, colour = "black") +
  geom_line(data = test_both_melt,aes(x=var,y=value),size = 4, colour = colours$grey, linetype = "dashed") +
  #scale_y_log10(limits = c(1000,10000), breaks=c(1000,10000)) +
    scale_y_continuous(expand=c(0,0),limits = c(0,7000)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,50,100,150,200)) +
  theme_alex() +
  geom_vline(xintercept = 50, alpha = 0.8, linetype="dashed") +
  xlab("Distance between Term peak and 5' mono-p") +
  ylab("n of found TTS") +
  ggtitle("Found TTS based on distance")

