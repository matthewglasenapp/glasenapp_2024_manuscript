#!/usr/bin/env python

import os 
import json
import operator
from operator import itemgetter
import statistics
from statistics import mean
import csv
from joblib import Parallel, delayed
import gzip

# Number of threads to use for coverage depth analysis with mosdepth
threads = 4

# Specify the posterior probability cutoff (1-100) for declaring a site as species_tree
posterior_probability_threshold = 90

# Specify the directory for output files 
output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm_species_tree/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

# Specify the directory with the output from aggregate_output.py.
# This directory contains json files of species_tree probability arrays for each scaffold. 
probability_file_dir = "/hb/scratch/mglasena/phylonet_hmm/probability_files_species_tree/"

# Specify the directory containing the scaffold alignments of all sites (including monomorphic sites)
all_sites_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments_all_sites/"

# Specify the directory contianing the scaffold alignments and coordinate files used in the PhyloNet-HMM analysis
phylonet_hmm_alignment_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments/"

# Specify the tsv file with scaffold names and length in base pairs
scaffold_info_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/scaffolds.tsv"

# Specify the file identifying any 100 kb gaps present in the nexus alignments 
gaps_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/25kb_gaps.bed"

# Specify the paths to the BAM files for each species used for assessing coverage depth
# Must be in alphabetical order by species name 
bam_file_paths_list = [
"/hb/home/mglasena/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

# Initialize a list for species_tree tract lists in format of [[start_coordinate, stop_coordinate, length in bp], []]
combined_results = []

# Initialize a dictionary containing scaffold metrics in the format of: 
# {scaffold: [length in base pairs, length spanned by first and last coordinates, number of sites in the all sites alignment, number nexus sites in PhyloNet-HMM alignment, number posterior probabilities > threshold, number of species_tree tracts, number of 10kb species_tree tracts]}
scaffold_dict = {}

# Initialize a dictionary for coverage depth by species of all species_tree tracts identified by PhyloNet-HMM
coverage_dict = dict()

# Initialize a dictionary for coverage depth by species for species_tree tracts passing coverage depth filters 
filtered_coverage_dict = dict()

# Return a list of lists of each scaffold posterior probability file matched with its respective coordinate file
def get_file_paths_pairs_list():
	
	# Create file of paths of global coordinate files for each scaffold. Save this file as "coorindate_file list" 
	os.system('find ' + phylonet_hmm_alignment_dir + ' -name "coordinates" -type f > coordinate_file_list')
		
	# Create file of paths to posterior probability files for each scaffold
	os.system('find ' + probability_file_dir + ' -type f > output_file_list')
	
	# Get zipped list matching the probability file path to the global coordinate file path for each scaffold. [[Scaffold_1.json path, Scaffold_1 coordinates path], []]
	with open("output_file_list","r") as f1:
		output_file_path_list = f1.read().splitlines()
	
	with open("coordinate_file_list","r") as f2:
		coordinate_file_path_list = f2.read().splitlines()
	
	sorted_zipped_list = [list(n) for n in zip(sorted(output_file_path_list),sorted(coordinate_file_path_list))]

	os.system('rm coordinate_file_list')
	os.system('rm output_file_list')
	
	return sorted_zipped_list

# Locate species_tree tracts and add to combined_results
def find_tracts_scaffold(json_file_path, coordinate_file_path):
	coordinates = open(coordinate_file_path,"r").read().splitlines()
	species_tree_probabilities = [probability * 100 for probability in json.load(open(json_file_path,"r"))]

	# Initialize a list of lists of species_tree tracts in the format of [[start index, stop index], []]
	index_tract_list = []
	
	# Initialize a list of lists of species_tree tracts in the format of [[start coordinate, stop coordinate, length, []]
	coordinate_tract_list = []
	
	# Find indices of species_tree tracts and add them to index_tract_list in the form of [Start index, Stop index]
	start = 0
	stop = 0
	i = 0
	while i < len(species_tree_probabilities):
		if species_tree_probabilities[i] >= posterior_probability_threshold:
			start = i
			stop = i
			while i < len(species_tree_probabilities) and species_tree_probabilities[i] >= posterior_probability_threshold:
				stop = i
				i = i+1
			index_tract_list.append([start, stop])
		i=i+1
	
	# Append each species_tree tract from index_tract_list to coordinate_tract_list in format of [start_coordinate, stop_coordinate, length in base pairs] 
	for tract in index_tract_list:
		# Format chr:coordinate. These coordinates are from the vcf file which is 1-based. Need to convert to zero-based
		start = coordinates[tract[0]]
		stop = coordinates[tract[1]]
		chrm = start.split(":")[0]

		# Convert to zero-based
		# Only need to decrement start by 1 because in bed files, the start position is 0-based. 
		start_coordinate = int(start.split(":")[1]) - 1
		stop_coordinate = int(stop.split(":")[1])
		
		start = chrm + ":" + str(start_coordinate)

		length = stop_coordinate - start_coordinate 
		
		if length > 1:
			coordinate_tract_list.append([start, stop, length])
	
	# Sort coordinate_tract_list by index 2 of each list (length in bp) in order of highest to lowest 
	coordinate_tract_list.sort(reverse=True, key=itemgetter(2))

	combined_results.append(coordinate_tract_list)

# Write each species_tree tract to a bed file in the format of scaffold, start, stop, name (scaffold:start_stop)
def write_tracts_to_bed(file_name, tract_list):
	total_number_tracts = len(tract_list)
	print("{} species_tree tracts found and written to tracts.bed".format(total_number_tracts))
	print("\n")
	with open(file_name,"w") as bed_file:
		for tract in tract_list:
			scaffold = str(tract).split(":")[0].split("'")[1]
			start = str(tract).split(":")[1].split("'")[0]
			stop = str(tract[1]).split(":")[1].split("'")[0]
			bed_file.write(scaffold + "\t" + start + "\t" + stop + "\t" + scaffold + ":" + start + "_" + stop + "\n")

# Intersect species_tree tract file with bed file containing positions of 100kb gaps 
def intersect_tract_file_with_100kb_gaps(tract_file):
	os.system("bedtools intersect -a " + tract_file + " -b " + gaps_file + " -wo > gap_overlap.bed")

# Trim or filter species_tree tracts if they overlap with a 100kb missing data gap 
# Write adjusted species_tree tracts to tracts_gap_filter.bed
def filter_tracts_for_gaps(tract_file):
	gap_overlaps = open("gap_overlap.bed", "r").read().splitlines()
	number_gap_overlapped_tracts = len(gap_overlaps)
	print("{} species_tree tract(s) overlapped with a 100 kb gap".format(number_gap_overlapped_tracts))
	print("\n")

	gap_overlapped_tracts = dict()
	for item in gap_overlaps:
		tract = item.split("\t")[3]
		tract_start = item.split("\t")[1]
		tract_stop = item.split("\t")[2]
		gap_start = item.split("\t")[5]
		gap_stop = item.split("\t")[6]
		gap_name = item.split("\t")[7]
		gap_overlapped_tracts[tract] = [tract_start, tract_stop, gap_start, gap_stop, gap_name]
		
	tracts = open(tract_file,"r").read().splitlines()
	
	tracts_written_to_bed = 0 
	tracts_adjusted = 0
	tracts_filtered = 0
	tracts_added = 0
	
	with open("tracts_gap_filter.bed","w") as f2:
		for line in tracts:
			tract = line.split("\t")[3]
			if tract not in gap_overlapped_tracts:
				f2.write(line + "\n")
				tracts_written_to_bed += 1
			
			else:
				print("The tract {} overlapped with the 100 kb gap {}".format(tract, gap_overlapped_tracts[tract][4]))
				scaffold = tract.split(":")[0]
				tract_start = gap_overlapped_tracts[tract][0]
				tract_stop = gap_overlapped_tracts[tract][1]
				gap_start = gap_overlapped_tracts[tract][2]
				gap_stop = gap_overlapped_tracts[tract][3]

				new_tract = False
				new_tract_2 = False

				# When converting a gap start to a tract stop, have to add 1 to the coordinate value 
				# When converting a gap stop to a tract start, have to subtract 1 to the coordinate value 
				if gap_start > tract_start and gap_stop < tract_stop:
					new_tract = [tract_start, str(int(gap_start) + 1)]
					new_tract_2 = [str(int(gap_stop) - 1), tract_stop]
					print("The gap {} splits the tract {} in two.".format(gap_overlapped_tracts[tract][4], tract))
					print("Splitting tract {} into two tracts to remove the gap {}.".format(tract, gap_overlapped_tracts[tract][4]))
					print("The new tracts are {} and {}".format(scaffold + ":" + new_tract[0] + "_" + new_tract[1], scaffold + ":" + new_tract_2[0] + "_" + new_tract_2[1]))
					print("\n")
					tracts_adjusted += 1
					tracts_added += 1

				elif gap_start > tract_start and gap_stop >= tract_stop:
					new_tract = [tract_start, str(int(gap_start) + 1)]
					print("The gap {} overlaps with the end of the tract {}".format(gap_overlapped_tracts[tract][4], tract))
					print("Trimming the end of tract {} from {} to {}".format(tract, tract_stop, new_tract[1]))
					print("The new tract is {}".format(scaffold + ":" + new_tract[0] + "_" + new_tract[1]))
					print("\n")
					tracts_adjusted += 1

				elif gap_start <= tract_start and gap_stop < tract_stop:
					new_tract = [str(int(gap_stop) - 1), tract_stop]
					print("The gap {} overlaps with the beginning of the tract {}".format(gap_overlapped_tracts[tract][4], tract))
					print("Trimming the beggining of tract{} from {} to {}".format(tract, tract_start, new_tract[0]))
					print("The new tract is {}".format(scaffold + ":" + new_tract[0] + "_" + new_tract[1]))
					print("\n")
					tracts_adjusted += 1

				elif gap_start <= tract_start and gap_stop >= tract_stop:
					print("The gap {} spans the entire tract {}".format(tract, gap_overlapped_tracts[tract][4]))
					print("Filtering {} from tracts".format(tract))
					print("\n")
					tracts_filtered += 1

				if new_tract:
					f2.write(scaffold + "\t" + new_tract[0] + "\t" + new_tract[1] + "\t" + scaffold + ":" + new_tract[0] + "_" + new_tract[1] + "\n")
					tracts_written_to_bed += 1

				if new_tract_2:
					f2.write(scaffold + "\t" + new_tract_2[0] + "\t" + new_tract_2[1] + "\t" + scaffold + ":" + new_tract_2[0] + "_" + new_tract_2[1] + "\n")
					tracts_written_to_bed += 1

	print("Started with {} species_tree tracts".format(len(tracts)))
	print("Tracts adjusted: {}".format(tracts_adjusted))
	print("Tracts filtered: {}".format(tracts_filtered))
	print("Tracts added: {}".format(tracts_added))
	print("Net change in tracts: {}".format(tracts_added - tracts_filtered))
	print("\n")
	print("{} tracts written to tracts_gap_filter.bed".format(tracts_written_to_bed))
	print("\n")

# Use mosdepth to get the mean coverage depth of each species_tree tract for each species 
def run_mosdepth(tract_file, bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by " + tract_file + " --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(threads, prefix, bam_file)
	os.system(mosdepth)

# Create a list of the relevant mosdepth output file paths in alphabetic order by species name 
def get_mosdepth_file_paths_pairs_list():
	find_regions_files = "find {} -type f -name '*.regions*' | grep -v 'csi' > mosdepth_regions_files".format(output_dir)
	os.system(find_regions_files)

	find_threshold_files = "find {} -type f -name '*.thresholds*' | grep -v 'csi' > mosdepth_threshold_files".format(output_dir)
	os.system(find_threshold_files)

	regions_files_list = open("mosdepth_regions_files","r").read().splitlines()
	threshold_files_list = open("mosdepth_threshold_files","r").read().splitlines()

	sorted_zipped_list = [list(n) for n in zip(sorted(regions_files_list),sorted(threshold_files_list))]

	os.system("rm mosdepth_regions_files")
	os.system("rm mosdepth_threshold_files")

	print(sorted_zipped_list)

	return sorted_zipped_list

# Populate the coverage_dict dictionary in the format of {tract_name: [coverage depth]}
def parse_mosdepth(regions_file):
	depths = gzip.open(regions_file,"rt").read().splitlines()
	
	for record in depths:
		tract_name = record.split("\t")[3]
		coverage = float(record.split("\t")[4])

		if tract_name not in coverage_dict:
			coverage_dict[tract_name] = [coverage]
		else:
			coverage_dict[tract_name].append(coverage)

def append_thresholds_by_species_to_coverage_dict(thresholds_file):
	thresholds = gzip.open(thresholds_file, "rt").read().splitlines()[1:]

	for record in thresholds:
		tract_name = record.split("\t")[3]
		length = int(record.split("\t")[2]) - int(record.split("\t")[1])
		prop_1x =  float(int(record.split("\t")[4]) / length)
		prop_10x = float(int(record.split("\t")[5]) / length)

		coverage_dict[tract_name].append(prop_1x)
		coverage_dict[tract_name].append(prop_10x)

# Removes all mosdepth outputfiles from output directory
def clean_up_mosdepth_output():
	os.system("rm *.dist.txt")
	os.system("rm *.summary.txt")
	os.system("rm *.gz")
	os.system("rm *.csi")

# Filter the coverage_dict for species_tree tracts where one sample has lower than 5x mean coverage depth or higher than 100x mean coverage depth. 
# Add species_tree tracts passing filter to filtered_coverage_dict. 
# Write tracts passing filter to tracts_pf.bed. Create tract_coverage.tsv with coverage depth for each species by species_tree tract.
def filter_tracts_by_depth():
	tracts_filtered = 0
	tracts_pass_filter = 0 
	tracts = open("tracts_gap_filter.bed","r").read().splitlines()

	for key in list(coverage_dict):
		depths = coverage_dict[key][0:4]
		if min(depths) < 5 or max(depths) >= 100:
			tracts_filtered += 1
			continue
		else:
			filtered_coverage_dict[key] = coverage_dict[key]

	with open("tracts_pf.bed","w") as f3:
		for tract in tracts:
			if tract.split("\t")[3] in filtered_coverage_dict:
				tracts_pass_filter += 1
				f3.write(tract + "\n")

	with open("tract_coverage.tsv","w") as f4:
		for key,value in filtered_coverage_dict.items():
			f4.write(key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\t" + str(value[4]) + "\t" + str(value[5]) + "\t" + str(value[6]) + "\t" + str(value[7]) + "\t" + str(value[8]) + "\t" + str(value[9]) + "\t" + str(value[10]) + "\t" + str(value[11]) + "\n")

	Sdro, Sfra, Spal, Hpul = ([] for i in range(4))

	for key,value in filtered_coverage_dict.items():
		Sdro.append(value[0])
		Sfra.append(value[1])
		Spal.append(value[2])
		Hpul.append(value[3])

	print("{} species_tree tracts were filtered due to low or too high coverage depth".format(tracts_filtered))
	print("{} species_tree tracts passed coverage depth filters and were written to tracts_pf.bed".format(tracts_pass_filter))
	print("Tracts passing filter and their coverage depth metrics were written to tract_coverage.tsv.")
	print("\n")

	print("Sdro mean coverage of species_tree tracts: {}".format(mean(Sdro)))
	print("Sfra mean coverage of species_tree tracts: {}".format(mean(Sfra)))
	print("Spal mean coverage of species_tree tracts: {}".format(mean(Spal)))
	print("Hpul mean coverage of species_tree tracts: {}".format(mean(Hpul)))
	print("\n")

def get_stats_on_tracts(tract_file):
	tracts = open(tract_file,"r").read().splitlines()

	length_dist = [int(tract.split("\t")[2]) - int(tract.split("\t")[1]) for tract in tracts]
	ten_kb_length_dist = [item for item in length_dist if item >= 10000]

	# Total number of species_tree tracts 
	num_tracts = len(length_dist)
	print("Number of species_tree tracts: {}".format(num_tracts))

	# Combined length of all species_tree tracts 
	total_length_species_tree_tracts = sum(length_dist)
	print("Combined length of the species_tree tracts: {}".format(total_length_species_tree_tracts))

	# Calculate mean, median and stdev of tract lengths 
	print("Mean species_tree tract length: {}".format(statistics.mean(length_dist)))
	print("Median species_tree tract length: {}".format(statistics.median(length_dist)))
	print("Standard Deviation for all species_tree tracts: {}".format(statistics.stdev(length_dist)))

	# Number of species_tree tracts 10 kb in length or greater 
	num_ten_kb_tracts = len(ten_kb_length_dist)
	print("Number of 10kb species_tree tracts: {}".format(num_ten_kb_tracts))

	# Combined length of 10kb species_tree tracts 
	total_length_ten_kb_species_tree_tracts = sum(ten_kb_length_dist)
	print("Combined length of the 10kb species_tree tracts: {}".format(total_length_ten_kb_species_tree_tracts))

	# Calculate mean, median, stdev of tracts larger than 10kb
	print("Mean tract length of species_tree tracts greater than 10 kb: {}".format(statistics.mean(ten_kb_length_dist)))
	print("Median tract length of species_tree tracts longer than 10 kb: {}".format(statistics.median(ten_kb_length_dist)))
	print("Standard deviation of species_tree tracts greater than 10kb: {}".format(statistics.stdev(ten_kb_length_dist)))

# Populate scaffold_dict dictionary in the format of {scaffold: [length in base pairs, length spanned by first and last coordinates, number of sites in the all sites alignment, number nexus sites in PhyloNet-HMM alignment, number posterior probabilities > threshold, number of species_tree tracts, number of 10kb species_tree tracts]}
def organize_tracts_by_scaffold(tract_file):
	scaffold_list = open(scaffold_info_file,"r").read().splitlines()

	for scaffold in scaffold_list:
		scaffold, length = scaffold.split("\t")
		scaffold_dict[scaffold] = [int(length)]
		first_coordinate = int(open(phylonet_hmm_alignment_dir + scaffold + "/" + "coordinates","r").read().splitlines()[0].split(":")[1]) - 1
		last_coordinate = int(open(phylonet_hmm_alignment_dir + scaffold + "/" + "coordinates","r").read().splitlines()[-1].split(":")[1])
		total_scaffold_alignment_length = last_coordinate - first_coordinate
		scaffold_dict[scaffold].append(total_scaffold_alignment_length)

	for scaffold in os.listdir(all_sites_dir):
		number_sites_all_sites = len(open(all_sites_dir + scaffold + "/" + "coordinates","r").read().splitlines())
		scaffold_dict[scaffold].append(number_sites_all_sites)

	for scaffold in os.listdir(probability_file_dir):
		scaffold_name = scaffold.split(".json")[0]
		number_sites = len(json.load(open(probability_file_dir + scaffold,"r")))
		number_sites_above_threshold = len([probability for probability in json.load(open(probability_file_dir + scaffold,"r")) if probability >= (posterior_probability_threshold/100)])
		scaffold_dict[scaffold_name].append(number_sites)
		scaffold_dict[scaffold_name].append(number_sites_above_threshold)
	
	# Initialize counters for number tracts and number ten_kb_tracts
	for scaffold in scaffold_dict:
		scaffold_dict[scaffold].append(0)
		scaffold_dict[scaffold].append(0)

	tracts = open(tract_file,"r").read().splitlines()
	for tract in tracts:
		scaffold = tract.split("\t")[0]
		scaffold_dict[scaffold][5] += 1
		length = int(tract.split("\t")[2]) - int(tract.split("\t")[1])
		if length >= 10000:
			scaffold_dict[scaffold][6] += 1

# Write the populated scaffold_dict dictionary to a csv file called species_tree_by_scaffold.csv
def write_scaffold_dict_to_csv():
	csv_file = open("species_tree_by_scaffold.csv","w")
	writer = csv.writer(csv_file)	
	header = ["scaffold", "length (bp)", "length spanned by first and last sites (bp)", "number of sites in all_sites alignment", "number nexus sites in phmm alignment", "number posterior probabilities > threshold", "number_tracts", "number_ten_kb_tracts"]
	writer.writerow(header)
	for scaffold in scaffold_dict:
		data = [scaffold, scaffold_dict[scaffold][0], scaffold_dict[scaffold][1],scaffold_dict[scaffold][2],scaffold_dict[scaffold][3],scaffold_dict[scaffold][4],scaffold_dict[scaffold][5],scaffold_dict[scaffold][6]]
		writer.writerow(data)
	csv_file.close()

# Write the distribution of species_tree tract lengths to a csv file called tract_length_dist.csv. Write the distribution of species_tree tract lengths for species_tree tracts 10kb and longer to tract_length_dist_10kb.csv.
def write_tract_length_dist():
	tracts = open("tracts_pf.bed").read().splitlines()
	tract_length_dist_all = [int(tract.split("\t")[2]) - int(tract.split("\t")[1]) for tract in tracts]
	tract_length_dist_10kb = [int(tract.split("\t")[2]) - int(tract.split("\t")[1]) for tract in tracts if int(tract.split("\t")[2]) - int(tract.split("\t")[1]) >= 10000]
	
	csv_file_all_tracts = open("tract_length_dist.csv","w")
	writer = csv.writer(csv_file_all_tracts)
	writer.writerow(tract_length_dist_all)
	csv_file_all_tracts.close()

	csv_file_ten_kb_tracts = open("tract_length_dist_10kb.csv","w")
	writer = csv.writer(csv_file_ten_kb_tracts)
	writer.writerow(tract_length_dist_10kb)
	csv_file_ten_kb_tracts.close()

# Write species_tree tracts 10kb and longer to ten_kb_tracts.bed
def write_ten_kb_tracts():
	tracts = open("tracts_pf.bed").read().splitlines()
	
	with open("ten_kb_tracts.bed", "w") as f2:
		for tract in tracts:
			if int(tract.split("\t")[2]) - int(tract.split("\t")[1]) >= 10000:
				f2.write(tract + "\n")
def main():
	os.chdir(output_dir)

	# Get a list of lists of [posterior probability file, scaffold coordinate file] for each scaffold
	files_by_scaffold_list = get_file_paths_pairs_list()

	# Process PhyloNet-HMM posterior probabilities from each scaffold
	for scaffold_file_pair in files_by_scaffold_list:
		json_file_path = scaffold_file_pair[0]
		coordinate_file_path = scaffold_file_pair[1]
		find_tracts_scaffold(json_file_path, coordinate_file_path)

	# Flatten the combined_results list of lists into a single list in format of [[start_coordinate, stop_coordinate, length in bp],[]]
	flattened_combined_results = [item for sublist in combined_results for item in sublist]
	
	# Sort flat_combined_results list of all species_tree tracts by tract length in base pairs in descending order. 
	flattened_combined_results.sort(reverse=True, key=itemgetter(2))

	# Write each species_tree tract to a bed file in the format of scaffold, start, stop, name (scaffold:start_stop)
	write_tracts_to_bed("tracts.bed", flattened_combined_results)

	# Use bedtools intersect to fine overlap between species_tree tracts and the 100 kb gaps in the scaffold alignments identified by find_gaps.py
	intersect_tract_file_with_100kb_gaps("tracts.bed")

	# Trim or filter species_tree tracts if they overlap with a 100kb missing data gap 
	# Write adjusted species_tree tracts to tracts_gap_filter.bed
	filter_tracts_for_gaps("tracts.bed")
	
	#Use mosdepth and tracts_gap_filter.bed to get the mean coverage depth of each species_tree tract for each species
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)("tracts_gap_filter.bed", bam_file) for bam_file in bam_file_paths_list)

	# Create a list of the relevant mosdepth output file paths in alphabetic order by species name 
	mosdepth_output_file_list = get_mosdepth_file_paths_pairs_list()

	# Populates the coverage_dict dictionary in the format of {tract_name: [coverage depth Sdro, coverage depth Sfra, coverage depth Spal, coverage depth Hpul]}
	for file_pair in mosdepth_output_file_list:
		regions_file = file_pair[0]
		parse_mosdepth(regions_file)

	for file_pair in mosdepth_output_file_list:
		thresholds_file = file_pair[1]
		append_thresholds_by_species_to_coverage_dict(thresholds_file)

	# Removes all mosdepth outputfiles from output directory
	clean_up_mosdepth_output()

	# Filter the coverage_dict for species_tree tracts where one sample has lower than 5x mean coverage depth or higher than 100x mean coverage depth. Add species_tree tracts passing filter to filtered_coverage_dict. Write tracts passing filter to tracts_pf.bed. Create tract_coverage.tsv with coverage depth for each species by species_tree tract.
	filter_tracts_by_depth()

	# Use mosdepth and tracts_pf.bed to get the mean coverage depth of each species_tree tracts by species for tracts passing filters that will be used in downstream analyses
	# Mosdepth is ran a second time to get metrics for only tracts passing filter 
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)("tracts_pf.bed", bam_file) for bam_file in bam_file_paths_list)

	# Print summary stats for species_tree tracts passing filter 
	get_stats_on_tracts("tracts_pf.bed")
	
	# Populate scaffold_dict dictionary in the format of {scaffold: [length in base pairs, length spanned by first and last coordinates, number of sites in the all sites alignment, number nexus sites in PhyloNet-HMM alignment, number posterior probabilities > threshold, number of species_tree tracts, number of 10kb species_tree tracts]}
	organize_tracts_by_scaffold("tracts_pf.bed")
	
	# Write the populated scaffold_dict dictionary to a csv file called species_tree_by_scaffold.csv
	write_scaffold_dict_to_csv()

	# Write the distribution of species_tree tract lengths to a csv file called tract_length_dist.csv. Write the distribution of species_tree tract lengths for species_tree tracts 10kb and longer to tract_length_dist_10kb.csv.
	write_tract_length_dist()
	
	# Write species_tree tracts 10kb and longer to ten_kb_tracts.bed
	write_ten_kb_tracts()

if __name__ == "__main__":
        main()
