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

# Directory where vcf2phylip was run
original_vcf2phylip_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/"

# Root directory where phylonet_hmm was run
root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm/"

# Directory with json files of introgression probability arrays from cat.py
probability_file_dir = "/hb/scratch/mglasena/phylonet_hmm/probability_files/"

# Directory containing the scaffold alignments of all sites 
all_sites_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/nexus_alignments_all_sites/"

# Directory contianing the scaffold alignments used in phylonet_hmm
phylonet_hmm_alignment_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments/"

# Tsv file with scaffold names and length in base pairs
scaffold_info_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/scaffolds.tsv"

# File containing 100 kb gaps in the nexus alignments 
gaps_file = "/hb/scratch/mglasena/phylonet_hmm/100kb_gaps.bed"

# Reference alignment BAM files for assessing coverage dpeth
bam_file_paths_list = [
"/hb/groups/pogson_group/dissertation/data/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

# List of coordinate_tract_lists for each scaffold in format of [[start_coordinate, stop_coordinate, length in bp], []]
combined_results = []

# Dictionary containing scaffolds in the format {scaffold: [length in base pairs, length spanned by first and last coordinates, number of sites in all_sites alignment, number nexus sites in phmm alignment, number posterior probabilities > 90, number_tracts, number_ten_kb_tracts]}
scaffold_dict = {}

# Dictionary for coverage of all introgression tracts identified by phylonet-hmm
coverage_dict = dict()

# Dictionary for coverage of filtered introgression tracts
filtered_coverage_dict = dict()

posterior_probability_threshold = 90

# Return list of phylonet_hmm output files matched with their respective coordinate files from the scaffold alignments
def get_file_paths_pairs_list():
	
	# Create file of paths of global coordinate files for each scaffold. Save this file as "coorindate_file list" 
	os.system('find ' + original_vcf2phylip_dir + ' -name "coordinates" -type f > coordinate_file_list')
		
	# Create file of paths of rawOutput.json files for each scaffold produced by phylonet_hmm. We only want the rawOutput.json file from the "bestrun" folder. Save this file as output_file_list in hmm directory
	os.system('find ' + probability_file_dir + ' -type f > output_file_list')
	
	# Get zipped list matching rawOutput.json file path to global coordinate file path for each scaffold. [[Scaffold_1_rawOutput.json path, Scaffold_1_coordinates], []]
	with open("output_file_list","r") as f1:
		output_file_path_list = f1.read().splitlines()
	
	with open("coordinate_file_list","r") as f2:
		coordinate_file_path_list = f2.read().splitlines()
	
	sorted_zipped_list = [list(n) for n in zip(sorted(output_file_path_list),sorted(coordinate_file_path_list))]

	os.system('rm coordinate_file_list')
	os.system('rm output_file_list')
	
	return sorted_zipped_list

# Locate introgression tracts and add to combined_results
def find_tracts_scaffold(json_file_path, coordinate_file_path):
	coordinates = open(coordinate_file_path,"r").read().splitlines()
	introgression_probabilities = [probability * 100 for probability in json.load(open(json_file_path,"r"))]

	# Create list of lists of introgression tracts in format [[start_index,stop_index,length], [start_index,stop_index,length]]
	index_tract_list = []
	coordinate_tract_list = []
	start = 0
	stop = 0
	# Find indexes of introgression tracts and add them to index_tract_list in the form of [Start index, Stop index]
	i = 0
	while (i < (len(introgression_probabilities) - 1)):
		if introgression_probabilities[i] >= posterior_probability_threshold:
			start = i
			stop = i
			variable_site_counter = 0
			while introgression_probabilities[i] >= posterior_probability_threshold and i < (len(introgression_probabilities)-1):
				variable_site_counter += 1
				stop = i
				i = i+1
			index_tract_list.append([start, stop, variable_site_counter])
		i=i+1
	
	# Append each introgression tract from index_tract_list to coordinate_tract_list in format of [start_coordinate, stop_coordinate, length in base pairs] 
	for tract in index_tract_list:
		# Format chr:coordinate. These coordinates are from the vcf file which is 1-based. Need to convert to zero-based
		start = coordinates[tract[0]]
		stop = coordinates[tract[1]]
		variable_sites = tract[2]

		chrm = start.split(":")[0]

		start_coordinate = int(start.split(":")[1])
		stop_coordinate = int(stop.split(":")[1])
		
		# Convert to zero-based
		# Only need to decrement start by 1 because in bed files, the start position is 0-based. 
		start_coordinate -= 1
		start = chrm + ":" + str(start_coordinate)

		length = stop_coordinate - start_coordinate 
		
		if length > 1:
			coordinate_tract_list.append([start, stop, length, variable_sites])
	
	# Sort coordinate_tract_list by index 2 of each list (length in bp) in order of highest to lowest 
	coordinate_tract_list.sort(reverse=True, key=itemgetter(2))

	combined_results.append(coordinate_tract_list)

# Write tract_list in format chr start stop
def write_tracts_to_bed(file_name, tract_list):
	with open(file_name,"w") as bed_file:
		for tract in tract_list:
			scaffold = str(tract).split(":")[0].split("'")[1]
			start = str(tract).split(":")[1].split("'")[0]
			stop = str(tract[1]).split(":")[1].split("'")[0]
			variable_sites = str(tract[3])
			bed_file.write(scaffold + "\t" + start + "\t" + stop + "\t" + scaffold + ":" + start + "_" + stop + "\n")

def write_sites_by_tract(tract_list):
	with open("SNV_sites_by_tract.tsv","w") as f:
		for tract in tract_list:
			scaffold = str(tract).split(":")[0].split("'")[1]
			start = str(tract).split(":")[1].split("'")[0]
			stop = str(tract[1]).split(":")[1].split("'")[0]
			variable_sites = str(tract[3])
			f.write(scaffold + "\t" + start + "\t" + stop + "\t" + scaffold + ":" + start + "_" + stop + "\t" + variable_sites + "\n")

# Intersect introgression tract file with bed file containing positions of 100kb gaps 
def intersect_tract_file_with_100kb_gaps(tract_file):
	os.system("bedtools intersect -a " + tract_file + " -b " + gaps_file + " -wo > gap_overlap.bed")

# Tracts overlapping with 100kb gaps
def filter_tracts_for_gaps(tract_file):
	gap_overlaps = open("gap_overlap.bed", "r").read().splitlines()

	gap_overlapped_tracts = dict()
	for item in gap_overlaps:
		tract = item.split("\t")[3]
		tract_start = item.split("\t")[1]
		tract_stop = item.split("\t")[2]
		gap_start = item.split("\t")[5]
		gap_stop = item.split("\t")[6]
		gap_overlapped_tracts[tract] = [tract_start, tract_stop, gap_start, gap_stop]
		
	tracts = open(tract_file,"r").read().splitlines()
	with open("tracts_gap_filter.bed","w") as f2:
		for line in tracts:
			tract = line.split("\t")[3]
			if tract not in gap_overlapped_tracts:
				f2.write(line + "\n")
			else:
				scaffold = tract.split(":")[0]
				tract_start = gap_overlapped_tracts[tract][0]
				tract_stop = gap_overlapped_tracts[tract][1]
				gap_start = gap_overlapped_tracts[tract][2]
				gap_stop = gap_overlapped_tracts[tract][3]

				new_tract = False
				new_tract_2 = False
				if gap_start > tract_start and gap_start < tract_stop:
					new_tract = [tract_start, gap_start]
					new_tract_2 = [gap_stop, tract_stop]

				elif gap_start > tract_start and gap_stop >= tract_stop:
					new_tract = [tract_start, gap_start]

				elif gap_start <= tract_start and gap_stop > tract_stop:
					new_tract = [gap_stop, tract_stop]

				elif gap_start < tract_start and gap_stop > tract_stop:
					continue

				elif gap_start == tract_start and gap_stop == tract_stop:
					continue

				if new_tract:
					f2.write(scaffold + "\t" + new_tract[0] + "\t" + new_tract[1] + "\t" + scaffold + ":" + new_tract[0] + "_" + new_tract[1] + "\n")

				if new_tract_2:
					f2.write(scaffold + "\t" + new_tract_2[0] + "\t" + new_tract_2[1] + "\t" + scaffold + ":" + new_tract_2[0] + "_" + new_tract_2[1] + "\n")

# Use mosdepth to get the mean coverage depth of each introgression tract for each species 
def run_mosdepth(tract_file, bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by " + tract_file + " --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(threads, prefix, bam_file)
	os.system(mosdepth)

# Get list of mosdepth output file paths in alphabetic order 
def get_mosdepth_output_file_list():
	find_files = "find /hb/scratch/mglasena/phylonet_hmm/ -type f -name '*.regions*' | grep -v 'csi' > mosdepth_output_files"
	os.system(find_files)
	mosdepth_output_file_list = open("mosdepth_output_files","r").read().splitlines()
	os.system("rm mosdepth_output_files")
	return sorted(mosdepth_output_file_list)

# Populate coverage_dict in the format of {tract_name: coverage_species1, coverage_species2, coverage_species3, coverage_species4}
def parse_mosdepth(mosdepth_region_file):
	with gzip.open(mosdepth_region_file,"rt") as f:
		inputs = f.read().splitlines()
	
	for record in inputs:
		tract_name = record.split("\t")[3]
		coverage = float(record.split("\t")[4])

		if tract_name not in coverage_dict:
			coverage_dict[tract_name] = [coverage]
		else:
			coverage_dict[tract_name].append(coverage)

# Filter coverage_dict for tracts where one sample has lower than 5x mean depth or higher than 100x mean depth
def filter_tracts_by_depth():
	tracts = open("tracts_gap_filter.bed","r").read().splitlines()

	for key in list(coverage_dict):
		if min(coverage_dict[key]) < 5 or max(coverage_dict[key]) >= 100:
			continue
		else:
			filtered_coverage_dict[key] = coverage_dict[key]

	with open("tracts_pf.bed","w") as f3:
		for tract in tracts:
			if tract.split("\t")[3] in filtered_coverage_dict:
				f3.write(tract + "\n")

	with open("tract_coverage.tsv","w") as f4:
		for key,value in filtered_coverage_dict.items():
			f4.write(key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\n")

	Sdro, Sfra, Spal, Hpul = ([] for i in range(4))

	for key,value in filtered_coverage_dict.items():
		Sdro.append(value[0])
		Sfra.append(value[1])
		Spal.append(value[2])
		Hpul.append(value[3])

	print("Sdro mean coverage of introgressed tracts: {}".format(mean(Sdro)))
	print("Sfra mean coverage of introgressed tracts: {}".format(mean(Sfra)))
	print("Spal mean coverage of introgressed tracts: {}".format(mean(Spal)))
	print("Hpul mean coverage of introgressed tracts: {}".format(mean(Hpul)))

def get_stats_on_tracts(tract_file):
	tracts = open(tract_file,"r").read().splitlines()

	length_dist = [int(tract.split("\t")[2]) - int(tract.split("\t")[1]) for tract in tracts]
	ten_kb_length_dist = [item for item in length_dist if item >= 10000]

	# Total number of introgression tracts 
	num_tracts = len(length_dist)
	print("Number of introgression tracts: {}".format(num_tracts))

	# Combined length of all introgression tracts 
	total_length_introgression_tracts = sum(length_dist)
	print("Combined length of the 10kb introgression tracts: {}".format(total_length_introgression_tracts))

	# Calculate mean, median and stdev of tract lengths 
	print("Mean introgression tract length: {}".format(statistics.mean(length_dist)))
	print("Median introgression tract length: {}".format(statistics.median(length_dist)))
	print("Standard Deviation for all introgression tracts: {}".format(statistics.stdev(length_dist)))

	# Number of introgression tracts 10 kb in length or greater 
	num_ten_kb_tracts = len(ten_kb_length_dist)
	print("Number of 10kb introgression tracts: {}".format(num_ten_kb_tracts))

	# Combined length of 10kb introgression tracts 
	total_length_ten_kb_introgression_tracts = sum(ten_kb_length_dist)
	print("Combined length of the 10kb introgression tracts: {}".format(total_length_ten_kb_introgression_tracts))

	# Calculate mean, median, stdev of tracts larger than 10kb
	print("Mean tract length of introgression tracts greater than 10 kb: {}".format(statistics.mean(ten_kb_length_dist)))
	print("Median tract length of introgression tracts longer than 10 kb: {}".format(statistics.median(ten_kb_length_dist)))
	print("Standard deviation of introgression tracts greater than 10kb: {}".format(statistics.stdev(ten_kb_length_dist)))

def organize_tracts_by_scaffold(tract_file):
	scaffold_list = open(scaffold_info_file,"r").read().splitlines()

	for scaffold in scaffold_list:
		scaffold, length = scaffold.split("\t")
		scaffold_dict[scaffold] = [int(length)]

	for scaffold in scaffold_list:
		scaffold = scaffold.split("\t")[0]
		first_coordinate = int(open(phylonet_hmm_alignment_dir + scaffold + "/" + "coordinates","r").read().splitlines()[0].split(":")[1]) - 1
		last_coordinate = int(open(phylonet_hmm_alignment_dir + scaffold + "/" + "coordinates","r").read().splitlines()[-1].split(":")[1])
		total_scaffold_alignment_length = last_coordinate - first_coordinate
		scaffold_dict[scaffold].append(total_scaffold_alignment_length)

	for scaffold_alignment_file in os.listdir(all_sites_dir):
		number_sites_all_sites = len(open(all_sites_dir + scaffold_alignment_file + "/" + "coordinates","r").read().splitlines())
		scaffold_dict[scaffold_alignment_file].append(number_sites_all_sites)

	for probability_file in os.listdir(probability_file_dir):
		scaffold = probability_file.split(".json")[0]
		number_sites = len(json.load(open(probability_file_dir + probability_file,"r")))
		number_sites_above_threshold = len([probability for probability in json.load(open(probability_file_dir + probability_file,"r")) if probability >= 0.9])
		scaffold_dict[scaffold].append(number_sites)
		scaffold_dict[scaffold].append(number_sites_above_threshold)
	
	# Initiate counters for number tracts and number ten_kb_tracts
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

def write_scaffold_dict_to_csv():
	csv_file = open("introgression_by_scaffold.csv","w")
	writer = csv.writer(csv_file)	
	header = ["scaffold", "length (bp)", "length spanned by first and last sites (bp)", "number of sites in all_sites alignment", "number nexus sites in phmm alignment", "number posterior probabilities > 90", "number_tracts", "number_ten_kb_tracts"]
	writer.writerow(header)
	for scaffold in scaffold_dict:
		data = [scaffold, scaffold_dict[scaffold][0], scaffold_dict[scaffold][1],scaffold_dict[scaffold][2],scaffold_dict[scaffold][3],scaffold_dict[scaffold][4],scaffold_dict[scaffold][5],scaffold_dict[scaffold][6]]
		writer.writerow(data)
	csv_file.close()

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

def write_ten_kb_tracts():
	tracts = open("tracts_pf.bed").read().splitlines()
	
	with open("ten_kb_tracts.bed", "w") as f2:
		for tract in tracts:
			if int(tract.split("\t")[2]) - int(tract.split("\t")[1]) >= 10000:
				f2.write(tract + "\n")
def main():
	# Get a list of lists of [phylonet_hmm output file, scaffold coordinate file]
	files_by_scaffold_list = get_file_paths_pairs_list()

	# Process phylonet_hmm output from each scaffold
	for scaffold_file_pair in files_by_scaffold_list:
		json_file_path = scaffold_file_pair[0]
		coordinate_file_path = scaffold_file_pair[1]
		find_tracts_scaffold(json_file_path, coordinate_file_path)

	# Flatten the combined_results list of lists into a single list in format of [[start_coordinate, stop_coordinate, length in bp],[]]
	flattened_combined_results = [item for sublist in combined_results for item in sublist]
	
	# Sort flat_combined_results list of all introgression tracts by tract length in base pairs in descending order. 
	flattened_combined_results.sort(reverse=True, key=itemgetter(2))

	write_tracts_to_bed("tracts.bed", flattened_combined_results)
	write_sites_by_tract(flattened_combined_results)

	intersect_tract_file_with_100kb_gaps("tracts.bed")

	filter_tracts_for_gaps("tracts.bed")
	
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)("tracts_gap_filter.bed", bam_file) for bam_file in bam_file_paths_list)

	mosdepth_output_file_list = get_mosdepth_output_file_list()

	for file in mosdepth_output_file_list:
		parse_mosdepth(file)

	filter_tracts_by_depth()

	get_stats_on_tracts("tracts_pf.bed")
	organize_tracts_by_scaffold("tracts_pf.bed")
	write_scaffold_dict_to_csv()

	write_tract_length_dist()
	write_ten_kb_tracts()

if __name__ == "__main__":
        main()
