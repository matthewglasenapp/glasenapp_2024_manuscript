import csv 
import random
from statistics import mean

# Number of bootstrap replicates
num_replicates = 1000

paml_output_file = "dNdS.tsv"

paml_lst = [item.split("\t") for item in open(paml_output_file,"r").read().splitlines()[1:]]

pseudoreplicate_dict = dict()

paml_output_dict = dict()

def create_paml_output_dict():
	for gene in paml_lst:
		paml_output_dict[gene[0]] = [gene[1], gene[2], gene[3]]

def get_random_gene():
	random_number = random.randint(0,len(paml_lst) - 1)
	random_gene = paml_lst[random_number][0]
	return random_gene

def create_pseudoreplicate():
	resampled_paml_lst = []
	for i in range(1, len(paml_lst) + 1):
		resampled_paml_lst.append(get_random_gene())
	return resampled_paml_lst

def create_pseudoreplicate_dict():
	for i in range(1, num_replicates + 1):
		key = "replicate_" + str(i)
		value = create_pseudoreplicate()
		pseudoreplicate_dict[key] = value

def get_distributions():
	distribution_dict = dict()

	mean_dN_dist = []
	mean_dS_dist = []
	mean_dNdS_dist = []
	
	for key,value in pseudoreplicate_dict.items():
		dN_list = []
		dS_list = []
		dNdS_list = []

		for gene in value:
			dN_list.append(float((paml_output_dict[gene][0])))
			dS_list.append(float(paml_output_dict[gene][1]))
			dNdS_list.append(float(paml_output_dict[gene][2]))

		mean_dN_dist.append(mean(dN_list))
		mean_dS_dist.append(mean(dS_list))
		mean_dNdS_dist.append(mean(dNdS_list))

	distribution_dict["mean_dN"] = mean_dN_dist
	distribution_dict["mean_dS"] = mean_dS_dist
	distribution_dict["mean_dNdS"] = mean_dNdS_dist

	return distribution_dict

def write_csv(output_file_name, dist):
	output_file = output_file_name + ".csv"
	output_csv = open(output_file, "w")
	writer = csv.writer(output_csv)
	value_lst = dist
	writer.writerow(value_lst)
	output_csv.close()

def main():
	create_paml_output_dict()
	create_pseudoreplicate_dict()
	
	distribution_dict = get_distributions()

	for key,value in distribution_dict.items():
		write_csv(key,value)

if __name__ == "__main__":
	main()