1) codeml.ctl is the codeml control file used as input for all PAML runs. 

2) compare_means_paml.R is an R script for visualizing the bootstrap comparison of means analysis for dN, dS, and dN/dS between introgressed and non-introgressed genes. 

3) extract.py parses the output of the PAML tests conducted in paml_introgressed_genes.py and paml_species_tree_genes.py. It also looks for genes with stop codons in the sequence alignments and removes these from consideration. 

4) paml_80.R is an R script that creates a plot for comparison of introgressed and non-introgressed metrics at the 80% posterior probability threshold. The metrics compared are dN, dS, dN/dS, percent coding, protein-coding genes / Mb, and positively selected genes / Mb. 

5) paml_90.R is an R script that creates a plot for comparison of introgressed and non-introgressed metrics at the 90% posterior probability threshold. The metrics compared are dN, dS, dN/dS, percent coding, protein-coding genes / Mb, and positively selected genes / Mb. 

6) paml_compare_means.py gets the distribution of dN, dS, and dN/dS values for both introgressed and non-introgressed genes and performs a bootstrap comparison of means. Genes with dN/dS > 1.5 are dropped from the analysis because they signal errors in reference alignment or genotyping. The datasets are also filtered to remove values that are more than 4 standard deviations from the mean (z-score > 4).  

7) paml_introgressed_genes.py finds all genes with more than half of their bases introgressed, filters the genes based on coverage depth metrics, builds multiple sequence alignment for those that pass filter, and runs codeml of PAML on each alignment. 

8) paml_species_tree_genes.py finds all genes with more than half of their bases confidently called for the species tree, filters the genes based on coverage depth metrics, builds multiple sequence alignment for those that pass filter, and runs codeml of PAML on each alignment. 