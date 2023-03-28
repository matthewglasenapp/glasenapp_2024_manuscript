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

def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by {} --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(original_tract_file, threads, prefix, bam_file)
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

def write_failed_tracts_and_report_mean_coverage_depth():
	with open("failed_coverage_tracts","a") as f:
		for key,value in coverage_dict.items():
			if min(value) < 5 or max(value) >= 100:
				f.write(key + "\n")
				del coverage_dict[key]

	Sfra = []
	Sdro = []
	Spal = []
	Hpul = []

	for key,value in coverage_dict.items():
		Sdro.append(value[0])
		Sfra.append(value[1])
		Spal.append(value[2])
		Hpul.append(value[3])

	print("Sfra mean coverage of introgressed tracts: {}".format(mean(Sfra)))
	print("Sdro mean coverage of introgressed tracts: {}".format(mean(Sdro)))
	print("Spal mean coverage of introgressed tracts: {}".format(mean(Spal)))
	print("Hpul mean coverage of introgressed tracts: {}".format(mean(Hpul)))

def rewrite_tract_file_with_failed_coverage_removed():
	with open(original_tract_file,"r") as f:
		tracts = f.read().splitlines()

	with open("failed_coverage_tracts","r") as f2:
		failed_tracts = f2.read().splitlines()

	with open("ten_kb_tracts_pf.bed","a") as f3:
		for tract in tracts:
			if tract.split("\t")[3] in failed_tracts:
					print("Failed tract!")
					continue
			else:
				f3.write(tract + "\n")

def intersect_tract_file_with_100kb_gaps():
	os.system("bedtools intersect -a " + original_tract_file + " -b " + gaps_file + " -wo > gap_overlap.bed")

def rewrite_tract_file_with_failed_gap_removed():
	with open("gap_overlap.bed", "r") as f:
		overlaps = f.read().splitlines()

	overlapped_tracts = set()
	for item in overlaps:
		tract = item.split("\t")[3]
		overlapped_tracts.add(tract)

	with open("ten_kb_tracts_pf.bed","r") as f2:
		tracts = f2.read().splitlines()

	with open("ten_kb_tracts_pfg.bed","a") as f3:
		for tract in tracts:
			if tract.split("\t")[3] in overlapped_tracts:
				continue
			else:
				f3.write(tract + "\n")

def main():
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

	mosdepth_output_file_list = sorted(get_mosdepth_output_file_list())

	print(mosdepth_output_file_list)


	for file in mosdepth_output_file_list:
		parse_mosdepth(file)

	write_failed_tracts_and_report_mean_coverage_depth()

	rewrite_tract_file_with_failed_coverage_removed()

	intersect_tract_file_with_100kb_gaps()

	print(coverage_dict)

if __name__ == "__main__":
	main()
