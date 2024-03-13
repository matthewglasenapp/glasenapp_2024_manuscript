1) aggregate_output.py aggregates the output from running PhyloNet-HMM 100 times on each scaffold. 

2) move_bestrun_files.py moves the bestrun files for each PhyloNet-HMM file for safe keeping. 

3) process_hmm.py and process_hmm_species_tree.py process the output of PhyloNet-HMM to identify the introgression tracts and species tree tracts present in the genome. 

4) verify_tracts.py is a sanity check for process_hmm.py. 