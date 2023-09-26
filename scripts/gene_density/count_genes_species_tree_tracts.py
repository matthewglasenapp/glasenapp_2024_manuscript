"""
Given a directory of replicate interval files:
1. Intersect each replicate interval file with protein_coding_genes.bed and unique_CDS.bed
2. Using the intersect bed file, count the number of overlapping protein coding genes that had >50% of bases overlapping
3. Output a csv file of the number of overlapping protein coding genes 
4. Using the intersect file from unique_exions.bed, count the number of overlapping coding bases for each replicate interval file
"""

from joblib import Parallel, delayed
import os
import csv 

# Number of cores available for Parallel computing
cores = 24

# Percentage of overlapping bases required for a gene to be counted.
percent_gene_overlap_threshold = 0.5

# Directory containing replicate interval (BED) files 
replicate_interval_file_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/null_intervals/"

# BED file containing all protein-coding genes in the S. purpuratus reference genome assembly
gene_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/protein_coding_genes.bed"

# BED file containing all unique protein-coding CDS in the S. purpuratus reference genome assembly 
CDS_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/nonoverlapping_unique_CDS.bed"

# Initialize dictionary to record the number of overlapping genes per replicate interval file
count_genes_dict = dict()

# Initialize dictionary to record the number of overlapping coding bases per replicate interval file 
count_bases_dict = dict()

# Create output directory for the output files resulting from running bedtools intersect with the replicate interval files and protein_coding_genes.bed
output_directory_genes = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/gene_density/species_tree_tracts/bedtools_intersect_dir_genes/"
make_output_directory_genes = "mkdir -p {}".format(output_directory_genes)
os.system(make_output_directory_genes)

# Create output directory for the output files resulting from running bedtools intersect with the replicate interval files and unique_CDS.bed
output_directory_CDS = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/gene_density/species_tree_tracts/bedtools_intersect_dir_CDS/"
make_output_directory_CDS = "mkdir -p {}".format(output_directory_CDS)
os.system(make_output_directory_CDS)

# Use bedtools intersect to create file bed file containing the overlap between an interval file and a feature file
def intersect(tract_file, overlap_file, output_directory):
	outfile = output_directory + tract_file.split("/")[-1]
	os.system("bedtools intersect -a " + tract_file + " -b " + overlap_file + " -wo > " + outfile)

def count_genes(intersect_file):
	overlaps = open(intersect_file,"r").read().splitlines()

	overlap_dict = dict()

	for line in overlaps:
		if line.split("\t")[7] in overlap_dict:
			overlap_dict[line.split("\t")[7]][1] += int(line.split("\t")[-1])
		else:
			overlap_dict[line.split("\t")[7]] = [((int(line.split("\t")[6])) - (int(line.split("\t")[5]))), int(line.split("\t")[-1])]

	filtered_overlap_dict = {key:value for key,value in overlap_dict.items() if value[1]/value[0] >= percent_gene_overlap_threshold}

	number_overlapping_genes = len(filtered_overlap_dict)

	return number_overlapping_genes

def add_to_count_genes_dict(intersect_file):
	replicate = intersect_file.split("/")[-1]
	count_genes_dict[replicate] = count_genes(intersect_file)

def count_coding_bases(intersect_file):
	overlaps = open(intersect_file,"r").read().splitlines()

	overlap_dict = dict()

	for line in overlaps:
		if line.split("\t")[7] in overlap_dict:
			overlap_dict[line.split("\t")[7]] += int(line.split("\t")[-1])
		else:
			overlap_dict[line.split("\t")[7]] = int(line.split("\t")[-1])

	total_coding_bases = sum(overlap_dict.values())

	return total_coding_bases

def add_to_count_bases_dict(intersect_file):
	replicate = intersect_file.split("/")[-1]
	count_bases_dict[replicate] = count_coding_bases(intersect_file)

def write_to_csv(dict, out_file):
	count_list = list(dict.values())

	output_csv = open(out_file,"w")
	writer = csv.writer(output_csv)
	writer.writerow(count_list)
	output_csv.close()

def main():
	# Get a list of replicate interval file paths
	replicate_file_path_lst = os.listdir(replicate_interval_file_dir)
	replicate_files = [replicate_interval_file_dir + item for item in replicate_file_path_lst]

	# Intersect each replicate interval file with protein_coding_genes.bed
	Parallel(n_jobs=cores)(delayed(intersect)(replicate, gene_file, output_directory_genes) for replicate in replicate_files)

	# Get a list of the output files from running bedtools intersect with protein_coding_genes.bed
	intersect_file_path_lst_genes = os.listdir(output_directory_genes)
	intersect_files_genes = [output_directory_genes + item for item in intersect_file_path_lst_genes]

	# Populate count_genes_dict with each replicate interval file and the number of overlapping protein-coding genes 
	for file in intersect_files_genes:
		add_to_count_genes_dict(file)

	# Write distribution of the number of overlapping protein-coding genes to csv 
	write_to_csv(count_genes_dict, "species_tree_tract_gene_count_90.csv")

	# Intersect each replicate interval file with protein_coding_genes.bed
	Parallel(n_jobs=cores)(delayed(intersect)(replicate, CDS_file, output_directory_CDS) for replicate in replicate_files)

	# Get a list of the output files from running bedtools intersect with unique_CDS.bed
	intersect_file_path_lst_CDS = os.listdir(output_directory_CDS)
	intersect_files_CDS = [output_directory_CDS + item for item in intersect_file_path_lst_CDS]

	# Populate count_CDS_dict with each replicate interval file and the number of overlapping protein-coding bases
	for file in intersect_files_CDS:
		add_to_count_bases_dict(file)

	# Write distribution of the number of overlapping protein-coding bases to csv 
	write_to_csv(count_bases_dict, "species_tree_tract_base_count_90.csv")

if __name__ == "__main__":
	main()
