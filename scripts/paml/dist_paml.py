import os
import csv 
from random import randint
import statistics

paml_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/"
introgressed_genes_dnds_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/introgressed_genes/"
non_introgressed_genes_dnds_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/species_tree_genes/"

distribution_dict = dict()

def get_dnds_distributions():
	introgressed_genes_dnds_file = introgressed_genes_dnds_dir + "dNdS.tsv"
	introgressed_genes_dN = [float(item.split("\t")[1]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:]]
	introgressed_genes_dS = [float(item.split("\t")[2]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:]]
	introgressed_genes_dNdS = [float(item.split("\t")[3]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:]]

	number_introgressed = len(introgressed_genes_dN)
	non_introgressed_genes_dnds_file = non_introgressed_genes_dnds_dir + "dNdS.tsv"
	non_introgressed_genes_dN = [float(item.split("\t")[1]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:]]
	non_introgressed_genes_dS = [float(item.split("\t")[2]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:]]
	non_introgressed_genes_dNdS = [float(item.split("\t")[3]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:]]

	subset_non_introgressed_genes_dN = []
	subset_non_introgressed_genes_dS = []
	subset_non_introgressed_genes_dNdS = []
	
	number_non_introgressed = len(non_introgressed_genes_dN)
	
	for i in range(1, number_introgressed + 1):
		random_number = randint(0, number_non_introgressed - 1)
		subset_non_introgressed_genes_dN.append(non_introgressed_genes_dN[random_number])
		subset_non_introgressed_genes_dS.append(non_introgressed_genes_dS[random_number])
		subset_non_introgressed_genes_dNdS.append(non_introgressed_genes_dNdS[random_number])

	distribution_dict["dN"] = [introgressed_genes_dN, subset_non_introgressed_genes_dN]
	distribution_dict["dS"] = [introgressed_genes_dS, subset_non_introgressed_genes_dS]
	distribution_dict["dNdS"] = [introgressed_genes_dNdS, subset_non_introgressed_genes_dNdS]

def write_csv(output_file_name, dist):
	if "non_introgressed" in output_file_name:
		output_dir = non_introgressed_genes_dnds_dir
	else:
		output_dir = introgressed_genes_dnds_dir
	output_file = output_dir + output_file_name + ".csv"
	output_csv = open(output_file, "w")
	writer = csv.writer(output_csv)
	value_lst = dist
	writer.writerow(value_lst)
	output_csv.close()

def run_approximate_permutation_test(metric, dist_lst, num_replicates):
	introgressed_lst = dist_lst[0]
	non_introgressed_lst = dist_lst[1]
	pooled_values = introgressed_lst + non_introgressed_lst
	num_pooled_values = len(introgressed_lst)
	
	difference_dist = []

	for i in range(1,num_replicates + 1):
		pseudoreplicate_one = []
		pseudoreplicate_two = []

		for i in range(0, num_pooled_values):
			random_number_one = randint(0, num_pooled_values - 1)
			value_one = pooled_values[random_number_one]
			pseudoreplicate_one.append(value_one)
			
			random_number_two = randint(0, num_pooled_values - 1)
			value_two = pooled_values[random_number_two]
			pseudoreplicate_two.append(value_two)
		
		mean_pseudoreplicate_one = statistics.mean(pseudoreplicate_one)
		mean_pseudoreplicate_two = statistics.mean(pseudoreplicate_two)
		difference = mean_pseudoreplicate_two - mean_pseudoreplicate_one
		difference_dist.append(difference)

	dif_csv = open(paml_output_dir + metric + "_difference_distribution_90.csv","w")
	dif_writer = csv.writer(dif_csv)
	dif_writer.writerow(difference_dist)
	dif_csv.close()

	test_statistic = statistics.mean(introgressed_lst) - statistics.mean(non_introgressed_lst)
	print("Mean Introgressed {}: {}".format(metric, statistics.mean(introgressed_lst)))
	print("Mean Non-Introgressed {}: {}".format(metric, statistics.mean(non_introgressed_lst)))
	print("True Difference: {}".format(test_statistic))

def main():
	get_dnds_distributions()

	for key,value in distribution_dict.items():
		introgressed_dist = value[0]
		write_csv(key + "_introgressed_90", introgressed_dist)
		non_introgressed_dist = value[1]
		write_csv(key + "_non_introgressed_90", non_introgressed_dist)

		run_approximate_permutation_test(key, value, 1000)

if __name__ == "__main__":
	main()