import csv 
import random
from statistics import mean

# Number of bootstrap replicates
num_replicates = 1000

pixy_output_file = "introgression_tracts_dxy.txt"

dxy_lst = [float(row.split("\t")[5]) for row in open(pixy_output_file,"r").read().splitlines()[1:] if row.split("\t")[5] != "NA"]

pseudoreplicate_dict = dict()

def get_random_value():
	random_number = random.randint(0,len(dxy_lst) - 1)
	random_value = dxy_lst[random_number]
	return random_value

def create_pseudoreplicate():
	resampled_dxy_lst = []
	for i in range(1, len(dxy_lst) + 1):
		resampled_dxy_lst.append(get_random_value())
	return resampled_dxy_lst

def create_pseudoreplicate_dict():
	for i in range(1, num_replicates + 1):
		key = "replicate_" + str(i)
		value = create_pseudoreplicate()
		pseudoreplicate_dict[key] = value

def get_dist_mean_dxy():
	mean_dxy_dist = []
	for value in pseudoreplicate_dict.values():
		mean_dxy_dist.append(mean(value))

	output_csv = open("introgression_tracts_mean_dxy_dist.csv","w")
	writer = csv.writer(output_csv)
	writer.writerow(mean_dxy_dist)
	output_csv.close()

def main():
	create_pseudoreplicate_dict()
	get_dist_mean_dxy()

if __name__ == "__main__":
	main()
