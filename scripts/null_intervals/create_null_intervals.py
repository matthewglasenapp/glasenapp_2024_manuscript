import os
from joblib import Parallel, delayed
import csv 
from statistics import mean

cores = 24

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/null_intervals/"

# Specify input files for bedtools shuffle
number_replicates = 1000
introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"
scaffold_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/scaffolds.tsv"
species_tree_regions = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm_species_tree/ten_kb_tracts.bed"

def get_replicate_file_lst():
	replicate_file_lst = ["replicate_" + str(i) for i in range(1,number_replicates + 1)]
	return replicate_file_lst

def shuffle(file):
	shuffle = "bedtools shuffle -i {} -g {} -incl {} -noOverlapping > {}{}".format(introgression_tract_file, scaffold_file, species_tree_regions, root_dir, file)
	os.system(shuffle)

def main():
	replicate_file_lst = get_replicate_file_lst()
	Parallel(n_jobs=cores)(delayed(shuffle)(replicate) for replicate in replicate_file_lst)

if __name__ == "__main__":
	main()
