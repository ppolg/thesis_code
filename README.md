# Complexities of transcriptional and post-transcriptional control in mycobacteria

Code associated with my PhD thesis, regarding analysis of short overlaps in mycobacteria.

Previously published code, associated with analysis of Term-seq is available at https://github.com/ppolg/Mtb_termseq

Code assoicated with analysis of nanopore sequencing is available at https://github.com/ppolg/mycopore






NUGA_multispecies_2:

•	fix_codon_table(): reads in the raw download of codon table usage from DNA HIVE and reformats it to be usable in downstream steps. This is conditionally called in the final wrapper or the initial read function to allow for usage of files from DNA HIVE directly. Returns a codon usage table as a data.frame object.

•	get_gff_input(): gets a GFF input from a folder and conditionally filters it. Only used as part of the following function.

•	read_files(): reads a FASTA, GFF and codon usage table files from a subfolder of the variable read_path, passed as an argument. GFF3 (or other formats) are further processed through regex string searches on its “attributes” column to find gene name, locus name and biotype. Codons are classified as acidic/basic/rare/start based on their encoded amino acid or frequency in coding sequences (with a threshold of 5 codons per thousand for rare codons). Returns a list of the three files.

•	annotate_overlaps(): appends the GFF file from a list from read_files(). Genome annotation data is split by strand and ordered by start to find the next ORF in sequence to annotate overlapping ORFs. The function also detects annotated 16S rRNA sequences to find the best fit for a pre-determined consensus anti Shine-Dalgarno sequence. This is piped to the Perl package ‘free2bind’ to determine predicted binding energy of the anti-Shine-Dalgarno to the upstream region of ORFs. The function returns the annotated GFF file and QC figures in addition to all the read files, as a list.

•	subset_overlaps(): reads the list from the previous function and creates subsets of the annotated GFF file. Overlap lengths are counted based on coordinates and overlapping sequences are subset based on overlap length, the presence or absence of a clear RBS and whether they are overlapped at the 5’ or 3’ of the ORFs. The function returns a list of these subsets.

•	draw_graphs(): reads in the subsetted list to plot figures (following some pre-processing). Created figures include distribution of the overlap lengths, the sequence distribution of NTGA overlaps and start codons as pie charts, and G and RBS distance distributions for the ORFs with each start codon. The created figures are both printed and returned as a list.

•	do_last10(): takes an annotated GFF file (from annotate_overlaps()) and the FASTA and codon rarity tables from read_files(), to perform a set of amino acid analyses on the last n (defined as a global variable) codons of ORFs overlapped at the 3’ end. Last n codons of every ORF, based on the GFF coordinates, are grabbed, and are translated to amino acids and ORFs ending in non-stop codons are filtered out, assuming an annotation error. Amino acids at every position as flagged as acidic, basic or rare, depending on the input tasks. Prolines and multi-proline stretches are also counted. The data is subsetted for only NUGA/URRUG overlaps, as well as based on the presence or absence of detected RBSs and for CGA-ending ORFs. Distributions in total and on a per-position base are then plotted, and amino acid distributions are visualised as logos.

•	do_first5(): functions similarly to do_last10(), but only does a pre-defined subset of the tasks (not taking a list of tasks as an argument), to analyse the first n (by standard 5, defined as a global variable) amino acids of ORFs overlapped at the 5’.

•	do_SD_QC(): takes the subsetted GFF files from subset_overlaps() and the codon analysis file from do_last10() to assess the quality of finding RBSs on the subsets. Subsets for NUGA, URRUG and CGA-ending overlapped ORFs are selected, as well as subsets based on the presence or absence of detected RBSs. Previously calculated G values of binding and distances to predicted binding position of RBSs are collated into new tibbles and subsequently plotted as bar charts per subsetted group. The function returns the new tables and figures.

•	analyse_nt(): checks the nucleotide context at the overlaps of ORFs to check for patterns and nucleotide enrichment. First the function calculates the GC% of the FASTA file from the read_files() function. Following this, it runs two for loops (one for each strand) for a number of iterations defined by a global variable corresponding to the number of nucleotides to analyse on each side of the overlaps. Following the counting of nucleotides, the function plots the distribution of the nucleotides per position as a combined point and line graph for each subset and returns all calculated tables and figures.

•	check_multi_P_seq(): Analyses the proline codon distribution of the last codons of 1- and 4-nucleotide overlaps (read from the output of the do_last10() function). For 1nt, 4nt subsets both with and without detected RBS (as defined in the subset_overlaps() function), counts the number of the four proline-encoding codons (CCT, CCG, CCA, CCC) and plots their relative abundance in each subset onto a bar chart that is subsequently printed and returned. 

•	calc_NUGA_RER(): calculates the “relative enrichment ratio” (RER) of both URRUG and NUGA overlaps, which indicates the enrichment of the last codons of ORFs (just prior the stop codon) compared to all ORFs. RER is normalised to sequence exclusions (such as the absence of CUGA overlaps) and is expressed as a log2-fold enrichment. The function also performs hypergeometric statistical tests with BH FDR corrections and plots the enrichment distributions as bar charts and volcano plots.

•	show_stop(): finds the stop codons and creates figures of their distribution, both as a bar and as a pie chart.

•	perform_GOA(): performs gene ontology analysis based on the COG database. Reads the subset table and external files from the COG database. Finds gene functional groups based on the database and calculates and plots distribution of the functional groups for 1nt and 4nt overlaps both as frequencies and as log2-fold changes compared to all genes. Returns results as list.

•	save_all(): Reads all previous returned outputs and saves all tables and figures in a subfolder based on the species shorthand. 

•	analyse_overlaps(): Main function that runs all the previous functions sequentially. Also saves the full list as an RDS file for future reading and utilisation.
