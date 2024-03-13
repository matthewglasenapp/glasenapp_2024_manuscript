1) count_genes_introgression_tracts.py counts the number of protein coding genes and coding bases that overlap introgression tracts. The 10kb introgression tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results. 

2) count_genes_species_tree_tracts.py counts the number of protein coding genes and coding bases that overlap species tree tracts. The 10kb species tree tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results. 

3) count_psg_introgression_tracts.py counts the number of genes with a history of positive selection overlapping introgression tracts. The list of genes with a history of positive selection can be found at data/spur_genome_metadata/psg.csv. The 10kb introgression tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results. 

4) count_psg_species_tracts.py counts the number of genes with a history of positive selection overlapping species tree tracts. The list of genes with a history of positive selection can be found at data/spur_genome_metadata/psg.csv. The 10kb species tree tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results.

5) create_null_intervals.py randomly permutes the 10 kb introgression tracts into regions confidently called for the species tree, creating 1,000 replicate interval files. The 10kb species tree tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results.

6) get_gene_density_windows.py divides the genome into windows of a specified size and calculate the propotion of the bases in that window that are part of a protein-coding gene. 

7) get_scaffold_gene_density.py calculates the proportion of each scaffold that is coding. 

8) resample_introgression_tracts.py bootstrap resamples the 10kb introgression tracts with replacement to create 1,000 pseudoreplicate interval files. These files are used as the input to count_genes_introgression_tracts.py and count_psg_introgression_tracts.py to get a distribution of the number of overlapping protein-coding genes, genes with a history of positive selection, and coding bases. The 10kb introgression tract interval file can be found at data/phylonet_hmm/phylonet_hmm_results. 