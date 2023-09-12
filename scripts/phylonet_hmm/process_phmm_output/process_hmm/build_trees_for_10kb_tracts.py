import os 
from joblib import Parallel, delayed

num_cores = 24
root_dir = "/hb/scratch/mglasena/test_vcf2phylip/"
vcf2phylip_path = "/hb/home/mglasena/software/vcf2phylip/"
genotype_file = "/hb/home/mglasena/dissertation/data/genotypes/franciscanus/3bp_filtered_genotype_calls_pf.g.vcf.gz"
introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"

def get_list_of_regions():
	regions = open(introgression_tract_file,"r").read().splitlines()
	regions = [region.split("\t")[3].replace("_","-").replace("-","_",1) for region in regions]
	return regions

def subset_vcf(region):
	output_dir = root_dir + region
	os.system("mkdir " + output_dir)
	output_file = output_dir + "/" + region + ".vcf.gz"
	view = "bcftools view -r {} -o {} {}".format(region, output_file, genotype_file)
	os.system(view)

def create_fasta_alignments(region):
	input_file = root_dir + region + "/" + region + ".vcf.gz"
	output_dir = root_dir + region
	run_vcf2phylip = "python3 {}vcf2phylip.py -w --input {} -p -f --output-folder {} --output-prefix {}".format(vcf2phylip_path, input_file, output_dir, region)
	os.system(run_vcf2phylip)

def main():
	regions = get_list_of_regions()

	Parallel(n_jobs=num_cores)(delayed(subset_vcf)(region) for region in regions)
	
	Parallel(n_jobs=num_cores)(delayed(create_fasta_alignments)(region) for region in regions)

if __name__ == "__main__":
	main()
