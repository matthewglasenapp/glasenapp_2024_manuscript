import os
import statistics
from random import randint
import csv

null_interval_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/null_intervals/"

introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"
scaffold_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/scaffolds.tsv"
species_tree_regions = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm_species_tree/ten_kb_tracts.bed"
pixy_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/pixy_dxy/"
it_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/pixy_dxy/introgression_tracts/"
st_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/pixy_dxy/species_tree_tracts/"

#Pixy input_files
vcf_file="/hb/scratch/mglasena/phylonet_hmm/hpul_sfra_invariant_sites_vcf/filtered_genotype_calls_individual_genotypes.g.vcf.gz"
pop_file="/hb/home/mglasena/dissertation/scripts/phylonet_hmm/pixy_dxy/popfile.txt"
cores="24"

# Bedtools Shuffle
def shuffle_introgression_tracts(input_file, output_file):
	shuffle = "bedtools shuffle -i {} -g {} -incl {} -noOverlapping > {}".format(input_file, scaffold_file, species_tree_regions, output_file)
	os.system(shuffle)

# Run pixy
def run_pixy(input_file, output_folder, output_prefix):
	run_pixy = "pixy --stats dxy --vcf {} --populations {} --bed_file {} --output_folder {} --output_prefix {} --n_cores {}".format(vcf_file, pop_file, input_file, output_folder, output_prefix, cores)
	os.system(run_pixy)

# If the dXY of some regions is "NA" after running pixy
# Drop that region and re-run shuffle to get another random region of equivalent length
# Re-run pixy
def check_dxy_dist_counts():
	# Get introgression tract dXY values
	introgression_tract_dxy_file = it_output_dir + "introgression_tracts_dxy.txt"
	introgression_tract_dxy_values = [float(item.split("\t")[5]) for item in open(introgression_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	# Get the species tree tract dXY values
	species_tree_tract_dxy_file = st_output_dir + "null_intervals_dxy.txt"
	species_tree_tract_dxy_values = [float(item.split("\t")[5]) for item in open(species_tree_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	while len(species_tree_tract_dxy_values) < len(introgression_tract_dxy_values):
		with open("redo.bed", "w") as f:
			for line in open(species_tree_tract_dxy_file,"r").read().splitlines()[1:]:
				if line.split("\t")[5] == 'NA':
					chrom = line.split("\t")[2]
					start = line.split("\t")[3]
					stop = line.split("\t")[4]
					f.write(chrom + "\t" + start + "\t" + stop + "\n")
		
		shuffle_introgression_tracts("redo.bed", "redo_null_interval_tracts.bed")

		original_lines = ["\t".join(item.split("\t")[0:3]) for item in open(null_interval_dir + "null_intervals.bed","r").read().splitlines()]
		redo_lines = open("redo.bed").read().splitlines()
		new_lines = open("redo_null_interval_tracts.bed", "r").read().splitlines()
		with open(null_interval_dir + "null_intervals.bed","w") as f:
			for line in original_lines:
				if not line in redo_lines:
					f.write(line + "\n")
			for line in new_lines:
				f.write(line + "\n")

		run_pixy(null_interval_dir + "null_intervals.bed", st_output_dir, "null_intervals")

		species_tree_tract_dxy_file = st_output_dir + "null_intervals_dxy.txt"
		species_tree_tract_dxy_values = [float(item.split("\t")[5]) for item in open(species_tree_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

		os.system('rm redo.bed')
		os.system('rm redo_null_interval_tracts.bed')

	# Write tract dxy distributions
	it_csv = open(it_output_dir + "introgression_tracts_dxy_dist_90.csv","w")
	it_writer = csv.writer(it_csv)
	it_writer.writerow(introgression_tract_dxy_values)
	it_csv.close()

	st_csv = open(st_output_dir + "species_tree_tracts_dxy_dist_90.csv","w")
	st_writer = csv.writer(st_csv)
	st_writer.writerow(species_tree_tract_dxy_values)
	st_csv.close()

def run_approximate_permutation_test(num_replicates):
	# Get introgression tract dXY values 
	introgression_tract_dxy_file = it_output_dir + "introgression_tracts_dxy.txt"
	introgression_tract_dxy_values = [float(item.split("\t")[5]) for item in open(introgression_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	# Get the species tree tract dXY values
	species_tree_tract_dxy_file = st_output_dir + "null_intervals_dxy.txt"
	species_tree_tract_dxy_values = [float(item.split("\t")[5]) for item in open(species_tree_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	# Pool all dXY values (both introgression tract and species tree tract)
	pooled_dxy = introgression_tract_dxy_values + species_tree_tract_dxy_values

	# Bootstrap resample the pooled values
	num_pooled_tracts = len(pooled_dxy)
	difference_dist = []

	for i in range(1,num_replicates + 1):
		pseudoreplicate_one = []
		pseudoreplicate_two = []

		for i in range(0, num_pooled_tracts):
			random_number_one = randint(0, num_pooled_tracts - 1)
			value_one = pooled_dxy[random_number_one]
			pseudoreplicate_one.append(value_one)
			
			random_number_two = randint(0, num_pooled_tracts - 1)
			value_two = pooled_dxy[random_number_two]
			pseudoreplicate_two.append(value_two)
		
		mean_pseudoreplicate_one = statistics.mean(pseudoreplicate_one)
		mean_pseudoreplicate_two = statistics.mean(pseudoreplicate_two)
		difference = mean_pseudoreplicate_two - mean_pseudoreplicate_one
		difference_dist.append(difference)

	dif_csv = open(pixy_output_dir + "difference_distribution_90.csv","w")
	dif_writer = csv.writer(dif_csv)
	dif_writer.writerow(difference_dist)
	dif_csv.close()

	test_statistic = statistics.mean(introgression_tract_dxy_values) - statistics.mean(species_tree_tract_dxy_values)
	print("Mean Introgression Tract dXY: {}".format(statistics.mean(introgression_tract_dxy_values)))
	print("Mean Species Tree Tract dXY: {}".format(statistics.mean(species_tree_tract_dxy_values)))
	print("True Difference: {}".format(test_statistic))

# Sanity Check!
def run_approximation_test_alternate(num_replicates):
	# Get introgression tract dXY values 
	introgression_tract_dxy_file = output_dir + "introgression_tracts_dxy.txt"
	introgression_tract_dxy_values = [float(item.split("\t")[5]) for item in open(introgression_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	# Get the species tree tract dXY values
	species_tree_tract_dxy_file = null_interval_dir + "null_intervals_dxy.txt"
	species_tree_tract_dxy_values = [float(item.split("\t")[5]) for item in open(species_tree_tract_dxy_file,"r").read().splitlines()[1:] if item.split("\t")[5] != 'NA']

	# Pool all dXY values (both introgression tract and species tree tract)
	pooled_dxy = introgression_tract_dxy_values + species_tree_tract_dxy_values

	# Bootstrap resample the pooled values
	num_pooled_tracts = len(pooled_dxy)

	group_a = []
	group_b = []

	for i in range(1,num_replicates + 1):
		pseudoreplicate_lst = []
		for i in range(0, num_pooled_tracts):
			random_number = randint(0, num_pooled_tracts - 1)
			value = pooled_dxy[random_number]
			pseudoreplicate_lst.append(value)
		group_a.append(statistics.mean(pseudoreplicate_lst))

	for i in range(1,num_replicates + 1):
		pseudoreplicate_lst = []
		for i in range(0, num_pooled_tracts):
			random_number = randint(0, num_pooled_tracts - 1)
			value = pooled_dxy[random_number]
			pseudoreplicate_lst.append(value)
		group_b.append(statistics.mean(pseudoreplicate_lst))

	dif_dist = []
	for i in range(0,num_replicates):
		value_1 = group_a[i]
		value_2 = group_b[i]
		dif_dist.append(value_2 - value_1)

	dif_csv = open("dif_alternate.csv","w")
	dif_writer = csv.writer(dif_csv)
	dif_writer.writerow(sorted(dif_dist))
	dif_csv.close()

def main():
	run_pixy(introgression_tract_file, it_output_dir, "introgression_tracts")
	shuffle_introgression_tracts(introgression_tract_file, null_interval_dir + "null_intervals.bed")
	run_pixy(null_interval_dir + "null_intervals.bed", st_output_dir, "null_intervals")
	check_dxy_dist_counts()
	run_approximate_permutation_test(1000)
	#run_approximation_test_alternate(1000)


if __name__ == "__main__":
	main()
