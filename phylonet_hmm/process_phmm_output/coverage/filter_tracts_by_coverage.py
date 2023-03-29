import os
import gzip
from joblib import Parallel, delayed
from statistics import mean

root_dir = "/hb/scratch/mglasena/phylonet_hmm/"

original_tract_file = "/hb/scratch/mglasena/phylonet_hmm/ten_kb_tracts.bed"

threads = 4

gaps_file = "/hb/scratch/mglasena/phylonet_hmm/100kb_gaps.bed"

bam_file_paths_list = [
"/hb/groups/pogson_group/dissertation/data/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

coverage_dict = dict()

filtered_coverage_dict = dict()

# Intersect introgression tract file with bed file containing positions of 100kb gaps 
def intersect_tract_file_with_100kb_gaps():
	os.system("bedtools intersect -a " + original_tract_file + " -b " + gaps_file + " -wo > gap_overlap.bed")

# Use mosdepth to get the mean coverage depth of each introgression tract for each species 
def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by {} --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(original_tract_file, threads, prefix, bam_file)
	os.system(mosdepth)

# Get list of mosdepth output file paths in alphabetic order 
def get_mosdepth_output_file_list():
	find_files = "find /hb/scratch/mglasena/phylonet_hmm/tract_coverage/ -type f -name '*.regions*' | grep -v 'csi' > mosdepth_output_files"
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

# Filter coverage_dict for the following conditions:
# Tracts where one sample has lower than 5x depth
# Tracts where one sample has > 100x depth
# Tracts overlapping with 100kb gaps
def filter_coverage_dict():
	# Identify tracts that overlap with 100kb gaps
	with open("gap_overlap.bed", "r") as f:
		gap_overlaps = f.read().splitlines()

	gap_overlapped_tracts = set()
	for item in gap_overlaps:
		tract = item.split("\t")[3]
		gap_overlapped_tracts.add(tract)
	
	for key in list(coverage_dict):
		if min(coverage_dict[key]) < 5 or max(coverage_dict[key]) >= 100 or key in gap_overlapped_tracts:
			continue
		else:
			filtered_coverage_dict[key] = coverage_dict[key]

	with open(original_tract_file,"r") as f2:
		tracts = f2.read().splitlines()

	with open("ten_kb_tracts_pf.bed","a") as f3:
		for tract in tracts:
			if tract.split("\t")[3] in filtered_coverage_dict:
				f3.write(tract + "\n")

	with open("tract_coverage.tsv","a") as f4:
		for key,value in filtered_coverage_dict.items():
			f4.write(key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\n")

	Sdro = []
	Sfra = []
	Spal = []
	Hpul = []

	for key,value in filtered_coverage_dict.items():
		Sdro.append(value[0])
		Sfra.append(value[1])
		Spal.append(value[2])
		Hpul.append(value[3])

	print("Sdro mean coverage of introgressed tracts: {}".format(mean(Sdro)))
	print("Sfra mean coverage of introgressed tracts: {}".format(mean(Sfra)))
	print("Spal mean coverage of introgressed tracts: {}".format(mean(Spal)))
	print("Hpul mean coverage of introgressed tracts: {}".format(mean(Hpul)))

def main():
	intersect_tract_file_with_100kb_gaps()

	#Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

	mosdepth_output_file_list = get_mosdepth_output_file_list()

	for file in mosdepth_output_file_list:
		parse_mosdepth(file)

	filter_coverage_dict()

if __name__ == "__main__":
	main()
