# Selection-Assay
An in vitro selection assay was devlepoved to design small gRNAs that can promote efficient A to I conversion using human ADAR2. 
This repository contains raw fastqc data for each round of the selection assay, as well as the script used to analyse the data and the raw result
Fastq files are labeled Xnt-Y-RZ where X is either 21 or 31 nt, Y is the round being analyzed and Z is one of the two pair-end read files. Round 0 corresponds to the unedited substrates and the Round C corresponds to the unrandomized substrates
The Diversity.py script is use to determine at each round the number of unique sequences encountered, the number of reads in which each unique guide sequence was present, the distribution of mismatches in the guide region, the distribution of errors outside the guide region and the base distribution at each position in the guide region. The results of running the script are in the 21nt and 31nt_raw_results files.
The Enrichment_rate.py script is use to determine how much a particular guide sequence was enriched from one round to the next. The results of running the script are in the 21nt and 31nt_Rate_analysis_9_rounds files.
