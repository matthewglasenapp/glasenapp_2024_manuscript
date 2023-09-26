from joblib import Parallel, delayed
import os

threads = 4

genome_metadata_dir = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/"

CDS_file = genome_metadata_dir + "nonoverlapping_unique_CDS.bed"

# Reference alignment BAM files for assessing coverage depth
bam_file_paths_list = [
"/hb/home/mglasena/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/home/mglasena/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

# Use mosdepth to get coverage depth metrics for introgressed genes by species 
def run_mosdepth(region_file, bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by " + region_file + " --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(threads, prefix, bam_file)
	os.system(mosdepth)

def main():
	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)(CDS_file, bam_file) for bam_file in bam_file_paths_list)

if __name__ == "__main__":
	main()

