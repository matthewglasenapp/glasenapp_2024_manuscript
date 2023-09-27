### Goal: Create 1,000 replicate datasets of ___ genes, the same number that were introgressed ###

import os
import csv 
import random
from statistics import mean

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/species_tree_genes/"

paml_output_file = root_dir + "dNdS.tsv"

# Number of bootstrap replicates
num_replicates = 1000

paml_lst = [item.split("\t") for item in open(paml_output_file,"r").read().splitlines()[1:]]

# Get the sample size of introgressed genes with paml metrics
sample_size_introgressed = len(open(paml_output_file,"r").read().splitlines()[1:])

replicate_dict = dict()

paml_output_dict = dict()

def create_paml_output_dict():
	for gene in paml_lst:
		paml_output_dict[gene[0]] = [gene[1], gene[2], gene[3]]

def get_random_gene():
	random_number = random.randint(0,len(paml_lst) - 1)
	random_gene = paml_lst[random_number][0]
	return random_gene

def create_replicate():
	resampled_paml_lst = []
	for i in range(1, sample_size_introgressed + 1):
		resampled_paml_lst.append(get_random_gene())
	return resampled_paml_lst

def create_replicate_dict():
	for i in range(1, num_replicates + 1):
		key = "replicate_" + str(i)
		value = create_replicate()
		replicate_dict[key] = value

def get_distributions():
	distribution_dict = dict()

	mean_dS_dist = []
	mean_dN_dist = []
	mean_dNdS_dist = []
	
	for key,value in replicate_dict.items():
		dS_list = []
		dN_list = []
		dNdS_list = []

		for gene in value:
			dS_list.append(float((paml_output_dict[gene][0])))
			dN_list.append(float(paml_output_dict[gene][1]))
			dNdS_list.append(float(paml_output_dict[gene][2]))

		mean_dS_dist.append(mean(dS_list))
		mean_dN_dist.append(mean(dN_list))
		mean_dNdS_dist.append(mean(dNdS_list))

	distribution_dict["mean_dS"] = mean_dS_dist
	distribution_dict["mean_dN"] = mean_dN_dist
	distribution_dict["mean_dNdS_dist"] = mean_dNdS_dist

	return distribution_dict

def write_csv(output_file_name, dist):
	output_file = output_file_name + ".csv"
	output_csv = open(output_file, "w")
	writer = csv.writer(output_csv)
	value_lst = dist
	writer.writerow(value_lst)
	output_csv.close()

def main():
	os.chdir(root_dir)
	create_paml_output_dict()
	create_replicate_dict()
	
	distribution_dict = get_distributions()

	for key,value in distribution_dict.items():
		write_csv(key,value)

if __name__ == "__main__":
	main()
