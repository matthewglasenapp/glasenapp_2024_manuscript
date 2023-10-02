import os
import csv 

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/gene_density/introgression_tracts/"

output_directory_genes = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/gene_density/introgression_tracts/bedtools_intersect_dir_genes/"

all_genes_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts/all_genes.csv"

psg_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/psg.csv"

LOC_SPU_dict = dict()

# Initialize dictionary to record the number of overlapping genes per replicate interval file
count_genes_dict = dict()

percent_gene_overlap_threshold = 0.5

def make_LOC_SPU_dict():
	all_genes_data = open(all_genes_file,"r").read().splitlines()
	all_genes_data = [record.split(",") for record in all_genes_data[1:]]
	for record in all_genes_data:
		if record[0] not in LOC_SPU_dict:
			if "SPU" in record[3]:
				LOC_SPU_dict[record[0]] = record[3].split("|")
			else:
				LOC_SPU_dict[record[0]] = "missing"
		else:
			print("Danger! Duplicate entry: {}".format(record[0]))

def count_overlapping_psg(intersect_file):
	overlaps = open(intersect_file,"r").read().splitlines()

	overlap_dict = dict()

	for line in overlaps:
		if line.split("\t")[7] in overlap_dict:
			overlap_dict[line.split("\t")[7]][1] += int(line.split("\t")[14])
		else:
			overlap_dict[line.split("\t")[7]] = [((int(line.split("\t")[6])) - (int(line.split("\t")[5]))), int(line.split("\t")[14])]

	filtered_overlap_dict = {key.split("gene-")[1]:value for key,value in overlap_dict.items() if value[1]/value[0] >= percent_gene_overlap_threshold}

	# Cross-reference for PSGs
	psg_list = list(csv.reader(open(psg_file,"r"), delimiter = ","))
	psg_list = [n[0] for n in psg_list[4:1012]]

	counter = 0
	for record in filtered_overlap_dict:
		if record in LOC_SPU_dict:
			for value in LOC_SPU_dict[record]:
				if value in psg_list:
					counter += 1
					break
	
	return counter

def add_to_count_genes_dict(intersect_file):
	replicate = intersect_file.split("/")[-1]
	count_genes_dict[replicate] = count_overlapping_psg(intersect_file)

def write_to_csv(dict, out_file):
	count_list = list(dict.values())

	output_csv = open(out_file,"w")
	writer = csv.writer(output_csv)
	writer.writerow(count_list)
	output_csv.close()

def main():
	os.chdir(root_dir)

	make_LOC_SPU_dict()
	
	# Get a list of the output files from running bedtools intersect with protein_coding_genes.bed
	intersect_file_path_lst_genes = os.listdir(output_directory_genes)
	intersect_files_genes = [output_directory_genes + item for item in intersect_file_path_lst_genes]

	# Populate count_genes_dict with each replicate interval file and the number of overlapping protein-coding genes 
	for file in intersect_files_genes:
		add_to_count_genes_dict(file)

	# Write distribution of the number of overlapping protein-coding genes to csv 
	write_to_csv(count_genes_dict, "dist_overlapping_psg.csv")

if __name__ == "__main__":
	main()
