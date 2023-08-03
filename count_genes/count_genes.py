from joblib import Parallel, delayed
import os
import csv 

cores = 24

percent_introgressed_threshold = 0.5

replicate_interval_file_dir = "/hb/scratch/mglasena/dxy/pixy/process_hmm_90/species_tree_tracts/replicate_interval_files/"

output_directory = "/hb/scratch/mglasena/dxy/pixy/process_hmm_90/species_tree_tracts/bedtools_intersect_dir/"
make_output_directory = "mkdir -p {}".format(output_directory)
os.system(make_output_directory)

gene_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/protein_coding_genes.bed"

count_genes_dict = dict()

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file):
	outfile = output_directory + tract_file.split("/")[-1] + ".bed"
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_file + " -wo > " + outfile)

def count_genes(intersect_file):
	overlaps = open(intersect_file,"r").read().splitlines()

	overlap_dict = dict()

	for line in overlaps:
		if line.split("\t")[7] in overlap_dict:
			overlap_dict[line.split("\t")[7]][1] += int(line.split("\t")[14])
		else:
			overlap_dict[line.split("\t")[7]] = [((int(line.split("\t")[6])) - (int(line.split("\t")[5]))), int(line.split("\t")[14])]

	filtered_overlap_dict = {key:value for key,value in overlap_dict.items() if value[1]/value[0] >= percent_introgressed_threshold}

	number_introgressed = len(filtered_overlap_dict)

	return number_introgressed

def add_to_count_genes_dict(intersect_file):
	replicate = intersect_file.split("/")[-1]
	count_genes_dict[replicate] = count_genes(intersect_file)

def write_gene_counts_to_csv():
	count_list = list(count_genes_dict.values())

	output_csv = open("test.csv","w")
	writer = csv.writer(output_csv)
	writer.writerow(count_list)
	output_csv.close()

def main():
	#replicate_file_path_lst = os.listdir(replicate_interval_file_dir)
	#replicate_files = [replicate_interval_file_dir + item for item in replicate_file_path_lst]

	#Parallel(n_jobs=cores)(delayed(intersect_genes)(replicate) for replicate in replicate_files)

	intersect_file_path_lst = os.listdir(output_directory)
	intersect_files = [output_directory + item for item in intersect_file_path_lst]

	for file in intersect_files:
		add_to_count_genes_dict(file)

	write_gene_counts_to_csv()

if __name__ == "__main__":
	main()
