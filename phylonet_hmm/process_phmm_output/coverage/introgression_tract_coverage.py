import os
import gzip
from joblib import Parallel, delayed

regions_file = "/hb/scratch/mglasena/phylonet_hmm/ten_kb_tracts.bed"
threads = 4

bam_file_paths_list = [
"/hb/groups/pogson_group/dissertation/data/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

coverage_dict = dict()

def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by {} --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(regions_file, threads, prefix, bam_file)
	os.system(mosdepth)

def get_mosdepth_output_file_list():
	find_files = "find /hb/scratch/mglasena/phylonet_hmm/tract_coverage/ -type f -name '*.regions*' | grep -v 'csi' > mosdepth_output_files"
	os.system(find_files)
	mosdepth_output_file_list = open("mosdepth_output_files","r").read().splitlines()
	os.system("rm mosdepth_output_files")
	return mosdepth_output_file_list

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

def write_failed_tracts():
	with open("failed_coverage_tracts","a") as f:
		for key,value in coverage_dict.items():
			if min(value) < 5:
				f.write(key + "\n")

def main():
	#Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

	mosdepth_output_file_list = get_mosdepth_output_file_list()

	for file in mosdepth_output_file_list:
		parse_mosdepth(file)

	write_failed_tracts()

if __name__ == "__main__":
	main()
