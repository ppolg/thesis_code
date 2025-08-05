# PhD Final Stretch: Check uORF NUGA PPE

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

uORF_count <- read.csv(here::here("NUGA_Out/uORF_test_all.csv"))
annot <- readRDS(here::here("NUGA_Out/species/Mtb/Mtb_full.rds"))


########
# Main #
########

gff_annot <- annot[["subset"]][["full"]]

k <- gff_annot$function_group[match(x=uORF_count$Rv_name,table=gff_annot$Rv_name)]
uORF_count$function_group <- k

# fix Rv2612c, filter non-unique so multiple uORFs don't distort
uORF_count <- uORF_count %>%
  dplyr::mutate(function_group = ifelse(Rv_name == "Rv2612c","lipid metabolism",function_group)) %>%
  filter(! base::duplicated(Rv_name))

uORF_NUGA <- uORF_count %>%
  filter(n_overlap_main == 4)


# in original?
gff_sum <- gff_annot %>%
  group_by(function_group) %>%
  summarise(n_gff = length(function_group)) %>%
  filter(function_group != "")

# in uORFs?
uORF_sum <- uORF_count %>%
  group_by(function_group) %>%
  summarise(n_all = length(function_group)) %>%
  add_row(function_group = "stable RNAs", n_all = 0) %>%
  add_row(function_group = "unknown", n_all = 0)

uORF_NUGA_sum <- uORF_NUGA %>%
  group_by(function_group) %>%
  summarise(n_NUGA = length(function_group)) %>%
  add_row(function_group = "stable RNAs", n_NUGA = 0) %>%
  add_row(function_group = "unknown", n_NUGA = 0)

# merge
sum_merged <- merge(gff_sum, uORF_sum, by="function_group")
sum_merged <- merge(sum_merged,uORF_NUGA_sum, by="function_group")

sum_merged <- sum_merged %>%
  mutate(r_gff = n_gff/nrow(gff_annot),
         r_all = n_all/nrow(uORF_count),
         r_NUGA = n_NUGA/nrow(uORF_NUGA))

# melt and plot
melt_merged <- sum_merged %>%
  select(1, 5:7) %>%
  dplyr::rename("All ORFs" = "r_gff",
         "ORFs with a uORF" = "r_all",
         "ORFs with a NUGA ORF" = "r_NUGA") %>%
  reshape2::melt()

fig_uORF_func <- ggplot(melt_merged,aes(x=value,y=function_group,fill=variable,colour=variable)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_mycopore() +
  xlab("portion of uORFs") +
  ylab("") +
  ggtitle("Functional distribution of ORFs with uORF") +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c(colours$purple,colours$blue,colours$green)) +
  scale_fill_manual(values=c(colours$purple,colours$blue,colours$green)) +
  theme(legend.title = element_blank())

print(fig_uORF_func)
