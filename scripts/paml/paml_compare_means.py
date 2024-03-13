import os
import csv 
import random
from random import randint
import statistics
from scipy.stats import zscore

paml_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/"
introgressed_genes_dnds_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/introgressed_genes/"
non_introgressed_genes_dnds_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/species_tree_genes/"

distribution_dict = dict()
dif_dist_dict = dict()

z_score_threshold = 4

def get_dnds_distributions():
	# Exclude genes where dNdS > 1.5. 
	introgressed_genes_dnds_file = introgressed_genes_dnds_dir + "dNdS.tsv"
	introgressed_genes_dN = [float(item.split("\t")[1]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	#print(introgressed_genes_dN)
	introgressed_genes_dS = [float(item.split("\t")[2]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	introgressed_genes_dNdS = [float(item.split("\t")[3]) for item in open(introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	#print(introgressed_genes_dNdS)
	#print(sorted(introgressed_genes_dNdS))

	# Remove severe outliers more than 4 standard deviations above the mean
	zscores_introgressed_genes_dN = zscore(introgressed_genes_dN).tolist()
	#print(zscores_introgressed_genes_dN)
	zscores_introgressed_genes_dS = zscore(introgressed_genes_dS).tolist()
	zscores_introgressed_genes_dNdS = zscore(introgressed_genes_dNdS).tolist()
	#print(zscores_introgressed_genes_dNdS)
	
	zscores_indices_dN = [index for index, item in enumerate(zscores_introgressed_genes_dN) if item > z_score_threshold]
	#print(zscores_indices_dN)
	zscores_indices_dS = [index for index, item in enumerate(zscores_introgressed_genes_dS) if item > z_score_threshold]
	#print(zscores_indices_dS)
	zscores_indices_dNdS = [index for index, item in enumerate(zscores_introgressed_genes_dNdS) if item > z_score_threshold]
	#print(zscores_indices_dNdS)
	
	indices_to_remove = sorted(list(set(zscores_indices_dN + zscores_indices_dS + zscores_indices_dNdS)))
	#print("Indices to remove: {}".format(indices_to_remove))
	#print(introgressed_genes_dNdS)
	
	introgressed_genes_dN_filtered = [introgressed_genes_dN[i] for i in range(len(introgressed_genes_dN)) if i not in indices_to_remove]
	introgressed_genes_dS_filtered = [introgressed_genes_dS[i] for i in range(len(introgressed_genes_dS)) if i not in indices_to_remove]
	introgressed_genes_dNdS_filtered = [introgressed_genes_dNdS[i] for i in range(len(introgressed_genes_dNdS)) if i not in indices_to_remove]

	number_introgressed = len(introgressed_genes_dN_filtered)
	print("Number introgressed: {}".format(number_introgressed))
	#print("Max introgressed dNdS: {}".format(max(introgressed_genes_dNdS_filtered)))
	
	# Get non-introgressed set of the same length
	non_introgressed_genes_dnds_file = non_introgressed_genes_dnds_dir + "dNdS.tsv"
	non_introgressed_genes_dN = [float(item.split("\t")[1]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	non_introgressed_genes_dS = [float(item.split("\t")[2]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	non_introgressed_genes_dNdS = [float(item.split("\t")[3]) for item in open(non_introgressed_genes_dnds_file,"r").read().splitlines()[1:] if float(item.split("\t")[3]) <= 1.5]
	#print("Total number non-introgressed to draw from: {}".format(len(non_introgressed_genes_dNdS)))

	subset_non_introgressed_genes_dN = []
	subset_non_introgressed_genes_dS = []
	subset_non_introgressed_genes_dNdS = []
	
	# Start with extra values so that after removing the outliers, there are still at least as many as number_introgressed
	number_non_introgressed = len(non_introgressed_genes_dN)
	random_indices = sorted(random.sample(range(0,number_non_introgressed),number_introgressed + 20))
	for random_number in random_indices:
		subset_non_introgressed_genes_dN.append(non_introgressed_genes_dN[random_number])
		subset_non_introgressed_genes_dS.append(non_introgressed_genes_dS[random_number])
		subset_non_introgressed_genes_dNdS.append(non_introgressed_genes_dNdS[random_number])

	zscores_non_introgressed_genes_dN = zscore(subset_non_introgressed_genes_dN).tolist()
	zscores_non_introgressed_genes_dS = zscore(subset_non_introgressed_genes_dS).tolist()
	zscores_non_introgressed_genes_dNdS = zscore(subset_non_introgressed_genes_dNdS).tolist()

	zscores_indices_non_introgressed_dN = [index for index, item in enumerate(zscores_non_introgressed_genes_dN) if item > z_score_threshold]
	zscores_indices_non_introgressed_dS = [index for index, item in enumerate(zscores_non_introgressed_genes_dS) if item > z_score_threshold]
	zscores_indices_non_introgressed_dNdS = [index for index, item in enumerate(zscores_non_introgressed_genes_dNdS) if item > z_score_threshold]
	
	indices_to_remove_non_introgressed = sorted(list(set(zscores_indices_non_introgressed_dN + zscores_indices_non_introgressed_dS + zscores_indices_non_introgressed_dNdS)))
	
	subset_non_introgressed_genes_dN_filtered = [subset_non_introgressed_genes_dN[i] for i in range(len(subset_non_introgressed_genes_dN)) if i not in indices_to_remove_non_introgressed]
	subset_non_introgressed_genes_dS_filtered = [subset_non_introgressed_genes_dS[i] for i in range(len(subset_non_introgressed_genes_dS)) if i not in indices_to_remove_non_introgressed]
	subset_non_introgressed_genes_dNdS_filtered = [subset_non_introgressed_genes_dNdS[i] for i in range(len(subset_non_introgressed_genes_dNdS)) if i not in indices_to_remove_non_introgressed]

	random_indices_2 = sorted(random.sample(range(0,len(subset_non_introgressed_genes_dNdS_filtered)),number_introgressed))
	final_non_introgressed_genes_dN = []
	final_non_introgressed_genes_dS = []
	final_non_introgressed_genes_dNdS = []
	
	for random_number in random_indices_2:
		final_non_introgressed_genes_dN.append(subset_non_introgressed_genes_dN_filtered[random_number])
		final_non_introgressed_genes_dS.append(subset_non_introgressed_genes_dS_filtered[random_number])
		final_non_introgressed_genes_dNdS.append(subset_non_introgressed_genes_dNdS_filtered[random_number])
	
	#print(final_non_introgressed_genes_dNdS)
	print("Number non-introgressed: {}".format(len(final_non_introgressed_genes_dNdS)))
	#print("Max non-introgressed dNdS: {}".format(max(final_non_introgressed_genes_dNdS)))
	distribution_dict["dN"] = [introgressed_genes_dN_filtered, final_non_introgressed_genes_dN]
	distribution_dict["dS"] = [introgressed_genes_dS_filtered, final_non_introgressed_genes_dS]
	distribution_dict["dNdS"] = [introgressed_genes_dNdS_filtered, final_non_introgressed_genes_dNdS]

def write_csv(output_file_name):
	output_file = paml_output_dir + output_file_name + ".csv"
	output_csv = open(output_file, "w")
	writer = csv.writer(output_csv)
	header = ["dist_type", "dN", "dS", "dNdS"]
	writer.writerow(header)

	dist_length = len(distribution_dict["dN"][0])
	for i in range(0, dist_length):
		dist_type = "introgressed"
		for key,value in distribution_dict.items():
			dN = distribution_dict["dN"][0][i]
			dS = distribution_dict["dS"][0][i]
			dNdS = distribution_dict["dNdS"][0][i]
		data = [dist_type, dN, dS, dNdS]
		writer.writerow(data)
	
	# Add for non_introgressed genes 
	for i in range(0, dist_length):
		dist_type = "non_introgressed"
		for key,value in distribution_dict.items():
			dN = distribution_dict["dN"][1][i]
			dS = distribution_dict["dS"][1][i]
			dNdS = distribution_dict["dNdS"][1][i]
		data = [dist_type, dN, dS, dNdS]
		writer.writerow(data)
	
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

		for j in range(0, num_pooled_values):
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

	dif_dist_dict[metric] = difference_dist


	test_statistic = statistics.mean(introgressed_lst) - statistics.mean(non_introgressed_lst)
	print("Mean Introgressed {}: {}".format(metric, statistics.mean(introgressed_lst)))
	print("Mean Non-Introgressed {}: {}".format(metric, statistics.mean(non_introgressed_lst)))
	print("True Difference: {}".format(test_statistic))
	#num_values_larger = len([item for item in difference_dist if item > test_statistic])
	#prop_values_larger = num_values_larger / len(difference_dist)
	#p_value = 1 - prop_values_larger
	proportion_values_larger = sum(value <= test_statistic for value in difference_dist) / len(difference_dist)
	print("p value: {}".format(proportion_values_larger))
	#print(proportion_values_larger)

def write_difference_dist_csv():
	dif_csv = open(paml_output_dir + "difference_distribution_90.csv","w")
	writer = csv.writer(dif_csv)
	header = ["dN", "dS", "dNdS"]
	writer.writerow(header)
	dN_dist = dif_dist_dict["dN"]
	dS_dist = dif_dist_dict["dS"]
	dNdS_dist = dif_dist_dict["dNdS"]
	for i in range(0,len(dN_dist)):
		data = [dN_dist[i], dS_dist[i], dNdS_dist[i]]
		writer.writerow(data)

	dif_csv.close()


def main():
	get_dnds_distributions()

	write_csv("paml_dist")

	for key,value in distribution_dict.items():
		run_approximate_permutation_test(key, value, 100000)

	#dNdS_dist = distribution_dict["dNdS"]
	#run_approximate_permutation_test("dNdS", dNdS_dist, 100000)

	write_difference_dist_csv()

if __name__ == "__main__":
	main()