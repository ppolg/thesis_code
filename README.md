# Complexities of transcriptional and post-transcriptional control in mycobacteria

## Repository information
Code associated with my PhD thesis, primarily regarding analysis of short overlaps in mycobacteria.

## Availability of other code
Previously published code, associated with analysis of Term-seq is available at https://github.com/ppolg/Mtb_termseq

Code assoicated with analysis of nanopore sequencing is available at https://github.com/ppolg/mycopore

The FreeAlign.pm perl module is Copyright (C) 2004, Joshua Starmer, used as part of the pipeline under the GNU General Public License.

## Code:

### NUGA_multispecies_2:

Baseline overlapping ORF analysis pipeline. Output per species is available in the NUGA_Out folder

•	_fix_codon_table()_: reads in the raw download of codon table usage from DNA HIVE and reformats it to be usable in downstream steps. This is conditionally called in the final wrapper or the initial read function to allow for usage of files from DNA HIVE directly. Returns a codon usage table as a data.frame object.

•	_get_gff_input()_: gets a GFF input from a folder and conditionally filters it. Only used as part of the following function.

•	_read_files()_: reads a FASTA, GFF and codon usage table files from a subfolder of the variable read_path, passed as an argument. GFF3 (or other formats) are further processed through regex string searches on its “attributes” column to find gene name, locus name and biotype. Codons are classified as acidic/basic/rare/start based on their encoded amino acid or frequency in coding sequences (with a threshold of 5 codons per thousand for rare codons). Returns a list of the three files.

•	_annotate_overlaps()_: appends the GFF file from a list from read_files(). Genome annotation data is split by strand and ordered by start to find the next ORF in sequence to annotate overlapping ORFs. The function also detects annotated 16S rRNA sequences to find the best fit for a pre-determined consensus anti Shine-Dalgarno sequence. This is piped to the Perl package ‘free2bind’ to determine predicted binding energy of the anti-Shine-Dalgarno to the upstream region of ORFs. The function returns the annotated GFF file and QC figures in addition to all the read files, as a list.

•	_subset_overlaps()_: reads the list from the previous function and creates subsets of the annotated GFF file. Overlap lengths are counted based on coordinates and overlapping sequences are subset based on overlap length, the presence or absence of a clear RBS and whether they are overlapped at the 5’ or 3’ of the ORFs. The function returns a list of these subsets.

•	_draw_graphs()_: reads in the subsetted list to plot figures (following some pre-processing). Created figures include distribution of the overlap lengths, the sequence distribution of NTGA overlaps and start codons as pie charts, and deltaG and RBS distance distributions for the ORFs with each start codon. The created figures are both printed and returned as a list.

•	_do_last10()_: takes an annotated GFF file (from _annotate_overlaps()_) and the FASTA and codon rarity tables from _read_files()_, to perform a set of amino acid analyses on the last n (defined as a global variable) codons of ORFs overlapped at the 3’ end. Last n codons of every ORF, based on the GFF coordinates, are grabbed, and are translated to amino acids and ORFs ending in non-stop codons are filtered out, assuming an annotation error. Amino acids at every position as flagged as acidic, basic or rare, depending on the input tasks. Prolines and multi-proline stretches are also counted. The data is subsetted for only NUGA/URRUG overlaps, as well as based on the presence or absence of detected RBSs and for CGA-ending ORFs. Distributions in total and on a per-position base are then plotted, and amino acid distributions are visualised as logos.

•	_do_first5()_: functions similarly to _do_last10()_, but only does a pre-defined subset of the tasks (not taking a list of tasks as an argument), to analyse the first n (by standard 5, defined as a global variable) amino acids of ORFs overlapped at the 5’.

•	_do_SD_QC()_: takes the subsetted GFF files from _subset_overlaps()_ and the codon analysis file from _do_last10()_ to assess the quality of finding RBSs on the subsets. Subsets for NUGA, URRUG and CGA-ending overlapped ORFs are selected, as well as subsets based on the presence or absence of detected RBSs. Previously calculated deltaG values of binding and distances to predicted binding position of RBSs are collated into new tibbles and subsequently plotted as bar charts per subsetted group. The function returns the new tables and figures.

•	_analyse_nt()_: checks the nucleotide context at the overlaps of ORFs to check for patterns and nucleotide enrichment. First the function calculates the GC% of the FASTA file from the _read_files()_ function. Following this, it runs two for loops (one for each strand) for a number of iterations defined by a global variable corresponding to the number of nucleotides to analyse on each side of the overlaps. Following the counting of nucleotides, the function plots the distribution of the nucleotides per position as a combined point and line graph for each subset and returns all calculated tables and figures.

•	_check_multi_P_seq()_: Analyses the proline codon distribution of the last codons of 1- and 4-nucleotide overlaps (read from the output of the do_last10() function). For 1nt, 4nt subsets both with and without detected RBS (as defined in the _subset_overlaps()_ function), counts the number of the four proline-encoding codons (CCT, CCG, CCA, CCC) and plots their relative abundance in each subset onto a bar chart that is subsequently printed and returned. 

•	_calc_NUGA_RER()_: calculates the “relative enrichment ratio” (RER) of both URRUG and NUGA overlaps, which indicates the enrichment of the last codons of ORFs (just prior the stop codon) compared to all ORFs. It is calculated as the ratio of specific codons in the last position for overlap _versus_ all ORFs.
RER is also normalised to sequence exclusions (such as the absence of CUGA overlaps) and is expressed as a non-log enrichment. The function also performs hypergeometric statistical tests with BH FDR corrections and plots the enrichment distributions as bar charts and volcano plots.

•	_show_stop()_: finds the stop codons and creates figures of their distribution, both as a bar and as a pie chart.

•	_perform_GOA()_: performs gene ontology analysis based on the COG database. Reads the subset table and external files from the COG database. Finds gene functional groups based on the database and calculates and plots distribution of the functional groups for 1nt and 4nt overlaps both as frequencies and as log2-fold changes compared to all genes. Returns results as list.

•	_save_all()_: Reads all previous returned outputs and saves all tables and figures in a subfolder based on the species shorthand. 

•	_analyse_overlaps()_: Main function that runs all the previous functions sequentially. Also saves the full list as an RDS file for future reading and utilisation.


### NUGA_compare2:

•	_read_full()_: Reads an RDS file from NUGA_multispecies_2. Also calculates GC content for downstream comparative analysis

•	_read_all()_: Reads multiple RDS files based on a list by calling _read_full()_, return a table of tables

•	_process_tables()_:Iterate through read RDS files (_i.e._ the table of tables from _read_all()_) to collate specific values into separate, merged tables. Essentially refactors previous datasets around comparable characteristics, to allow for plotting 

•	_plot_merged()_: Plots figures from _process_tables()_

•	_calc_NUGA_percents()_: Calculates NUGA percentages for overlaps (i.e. the % of overlaps that fit the subsets defined by _NUGA_multies_2.R_. Also makes comparative figures and prints them


### find_SD:
Compares and contrasts different methodologies of detecting and picking up RBSs. Tests include consecutive purines with/without mismatch, and different pairings to anti-SD. Also contains random controls. Plots QC figures including deltaG and distance to RBS.


### qRT_analysis:
Contains code related to brief analysis of qRT-PCR of RofZ KD. 


### mycopore_init:
Contains themes, colours, plotting settings, shorthands etc.




