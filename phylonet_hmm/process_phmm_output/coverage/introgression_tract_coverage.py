import os
from joblib import Parallel, delayed

regions_file = "/hb/scratch/mglasena/phylonet_hmm/ten_kb_tracts.bed"
threads = 4

bam_file_paths_list = [
"/hb/groups/pogson_group/dissertation/data/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

def run_mosdepth(bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by {} --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(regions_file, threads, prefix, bam_file)
	os.system(mosdepth)

def main():
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)(bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()
