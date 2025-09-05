
# PhD Final Stretch: Analysis for multiple organisms
#   • Quick qRT ddCT

###########################################################################

########
# Func #
########

theme_mycopore <- function(base_size=10) {
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

########
# Main #
########

# Get Ct from csv (DA2)

ct_table <- read.csv(paste0(here::here(),"/qRT/ct_sense_as.csv"), comment.char = "#") %>%
  filter(Quantity == "-") %>%
  select(1,2,6,7,8)

ct_table[is.na(ct_table)] <- 0

ct_table_2 <- read.csv(paste0(here::here(),"/qRT/ct_16S_RofZ.csv"), comment.char = "#") %>%
  filter(Quantity == "-") %>%
  select(1,2,6,7,8)

ct_table_2[is.na(ct_table_2)] <- 0

ct_table_sum <- rbind(ct_table,ct_table_2) %>%
  filter(Sample != "NTC") %>%
  mutate(time = str_extract(Sample,"[0-9]+")) %>%
  mutate(induced = ifelse(grepl("no",Sample,fixed = T),F,T)) %>%
  mutate(RT = ifelse(grepl("+",Sample,fixed = T),"+","-"))

ct_table_house<- ct_table_sum %>%
  filter(induced == T & RT == "+" & Target == "16S") %>%
  arrange(Target,as.numeric(time))

ct_table_work <- ct_table_sum %>%
  filter(induced == T & RT == "+" & Target != "16S") %>%
  arrange(Target,as.numeric(time))

# ct_table_house_control <- ct_table_sum %>%
#   filter(induced == F & RT == "+" & Target == "16S") %>%
#   arrange(Target,as.numeric(time))
# 
# ct_table_work_control <- ct_table_sum %>%
#   filter(induced == F & RT == "+" & Target != "16S") %>%
#   arrange(Target,as.numeric(time))

ct_table_work["Cq_house"] <- ct_table_house$Cq.Mean[1:4]
ct_table_work <- ct_table_work %>%
  mutate(deltaCQ = Cq.Mean-Cq_house)

for(i in 1:nrow(ct_table_work)){
  h <- ct_table_work$deltaCQ[i]
  k <- ceiling(i/4) + (ceiling(i/4)-1)*3
  j <- ct_table_work$deltaCQ[k]
  print(paste0("iteration: ", i, " original: ", h, " deducting: ",j, " result: ", h-j))
  ct_table_work$ddCq[i] <- (h-j)
}

ct_table_work <- ct_table_work %>%
  mutate(l2fc = log2(2^-ddCq)) #also equals -ddCq, but it's late so I am a dum dum :(

ct_table_long <- ct_table_work %>%
  select(2,6,12) %>%
  filter(time != 0) %>%
  mutate(time = factor(as.numeric(time))) %>%
  arrange(time) %>%
  reshape2::melt() %>%
  mutate(SE = c(0.46866974,
                0.10839800,
                0.06940776,
                0.28966054,
                0.05943836,
                0.04483783,
                0.53471308,
                0.17064600,
                0.07932335))
  
dodge_test <- position_dodge(width=0.9)

fig_l2fc <- ct_table_long %>%
  ggplot() +
  geom_bar(aes(x=time,y=value,fill=Target),position=dodge_test,stat="identity") +
  geom_errorbar(aes(x=time, ymin=value-SE, ymax=value+SE, fill=Target),
                position=dodge_test,width=0.5, colour="grey40") +
  xlab("time after induction (h)") +
  ylab("log2-fold change") +
  scale_fill_manual(values=c(colours$alex_R2,colours$alex_B2,"grey")) +
  theme_mycopore() +
  theme(legend.title  = element_blank())

print(fig_l2fc)
