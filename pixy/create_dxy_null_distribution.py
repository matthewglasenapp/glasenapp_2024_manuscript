import os
from joblib import Parallel, delayed
import csv 
from statistics import mean

cores = 24

root_dir = "/hb/scratch/mglasena/dxy/pixy/species_tree_tracts/"

# Specify input files for bedtools shuffle
number_replicates = 1000
introgression_tract_file = "ten_kb_tracts.bed"
scaffold_file = "scaffolds.txt"
species_tree_regions = "ten_kb_tracts_st.bed"

# Specify additional input files for pixy
vcf_file = "/hb/scratch/mglasena/dxy/vcf_files/combined_vcf_files/filtered_genotype_calls_individual_genotypes.g.vcf.gz"
popfile = "popfile.txt"

# Directory for replicate species tree interval files 
replicate_dir = root_dir + "replicate_interval_files/"
make_replicate_dir = "mkdir -p {}".format(replicate_dir)
os.system(make_replicate_dir)

# Directory for replicate species tree interval files 
pixy_dir = root_dir + "pixy_files/"
make_pixy_dir = "mkdir -p {}".format(pixy_dir)
os.system(make_pixy_dir)

def get_replicate_file_lst():
	replicate_file_lst = ["replicate_" + str(i) for i in range(1,number_replicates + 1)]
	return replicate_file_lst

def shuffle(file):
	shuffle = "bedtools shuffle -i {} -g {} -incl {} -noOverlapping > {}{}".format(introgression_tract_file, scaffold_file, species_tree_regions, replicate_dir, file)
	os.system(shuffle)

def run_pixy(file):
	input_file = replicate_dir + file
	run_pixy = "pixy --stats dxy --vcf {} --populations {} --bed_file {} --output_folder {} --output_prefix {}".format(vcf_file, popfile, input_file, pixy_dir, file)
	os.system(run_pixy)

def get_dist():
	mean_lst = []
	for file in os.listdir(pixy_dir):
		mean_dxy = mean([float(row.split("\t")[5]) for row in open(pixy_dir + file,"r").read().splitlines()[1:] if row.split("\t")[5] != "NA"])
		mean_lst.append(mean_dxy)

	output_csv = open("species_tree_tracts_mean_dxy.csv","w")
	writer = csv.writer(output_csv)
	writer.writerow(mean_lst)
	output_csv.close()
	
def main():
	replicate_file_lst = get_replicate_file_lst()
	Parallel(n_jobs=cores)(delayed(shuffle)(replicate) for replicate in replicate_file_lst)
	Parallel(n_jobs=cores)(delayed(run_pixy)(replicate) for replicate in replicate_file_lst)
	get_dist()

if __name__ == "__main__":
	main()