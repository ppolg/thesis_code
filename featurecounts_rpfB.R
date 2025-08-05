# FeatureCounts


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

# Get list of bam files

gff <- get_gff_input("Msm3")
directory <- paste(here::here("Brindhaseqcounts/bam/sorted"))
files <- grep("sorted.bam",list.files(directory),value=TRUE)
files <- files[!grepl("bai",files)]
names <- sub(".sorted.bam","\\1", files)



# Featurecounts loop
for(i in 1:length(files)){
  file <- files[i]
  name <- names[i]
  
  # Run featureCounts
  k <- featureCounts(allowMultiOverlap = T,
                     files = paste(directory,"/",file,sep=""),
                     annot.ext = paste(directory,"/Mtb2.gff",sep=""),
                     isGTFAnnotationFile = T,
                     GTF.featureType = "gene",
                     GTF.attrType = "ID",
                     GTF.attrType.extra = "Name",
                     nthreads = 8,
                     isPairedEnd = T)
  
  # Make df from output
  l <- tibble(k[["annotation"]]) %>%
    mutate(counts = k[["counts"]])
  
  # write counts matrix in csv
  outfile = paste0(directory,"/",name,".counts")
  fwrite(l, file = outfile)
}
