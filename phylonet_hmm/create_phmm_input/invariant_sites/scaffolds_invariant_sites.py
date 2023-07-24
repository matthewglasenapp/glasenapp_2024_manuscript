import os 
import subprocess
import multiprocessing
from joblib import Parallel, delayed

num_jobs = 21

# Get separate multi-sample vcf file for each scaffold. Get a file with the coordinate for each SNV site applied to the matrix.

# Need to manually input the samples to exlcude from the multisample vcf file in the subset_vcf_by_scaffold() function.  

reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# File with scaffold names, one per line
scaffold_list_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/scaffolds.txt"

filtered_vcf = "/hb/scratch/mglasena/invariant_sites_vcf_phylonet_hmm/combined_vcf_files/filtered_genotype_calls_individual_genotypes.g.vcf.gz"

# Path to vcf2phylip program
vcf2phylip_path = "/hb/home/mglasena/software/vcf2phylip/"

root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/"

vcf_by_scaffold_dir = root_dir + "vcf_by_scaffold_all_sites/"
make_vcf_by_scaffold_dir = "mkdir -p {}".format(vcf_by_scaffold_dir)
os.system(make_vcf_by_scaffold_dir)

scaffold_nexus_alignment_dir = root_dir + "scaffold_nexus_alignments_all_sites/"
make_scaffold_nexus_alignment_dir = "mkdir -p {}".format(scaffold_nexus_alignment_dir)
os.system(make_scaffold_nexus_alignment_dir)

outgroup_sample_name = "QB3KMK016"
number_species = "4"

def get_scaffold_list():
	with open(scaffold_list_file, "r") as f:
		scaffold_list = f.read().splitlines()
	return scaffold_list

# Get a vcf file for each scaffold 
def subset_vcf_by_scaffold(scaffold):
	output_dir = vcf_by_scaffold_dir
	subset_vcf_by_scaffold = "gatk SelectVariants -R {} -V {} --select-type-to-exclude MNP --select-type-to-exclude INDEL --select-type-to-exclude SYMBOLIC --select-type-to-exclude MIXED -O {}{}.vcf.gz -L {} --exclude-filtered true".format(reference_genome, filtered_vcf, output_dir, scaffold, scaffold)
	os.system(subset_vcf_by_scaffold)

# Run vcf2phylip on each individual scaffold
def convert_vcf_to_nexus(scaffold):
	input_dir = vcf_by_scaffold_dir
	output_dir = scaffold_nexus_alignment_dir + scaffold + "/"
	run_vcf2phylip = "python3 {}vcf2phylip.py -w --input {}{}.vcf.gz --min-samples-locus {} -p -n -r --outgroup {} --output-folder {} --output-prefix {}".format(vcf2phylip_path, input_dir, scaffold, number_species, outgroup_sample_name, output_dir, scaffold)
	
	os.system(run_vcf2phylip)

# Reformat vcf2phylip used sites output file for compatibility with phylonet_hmm output. 
def reformat_coordinate_files(scaffold):
	input_dir = scaffold_nexus_alignment_dir + scaffold
	os.chdir(input_dir)
	
	os.environ['coordinate_file'] = "{}.min{}.used_sites.tsv".format(scaffold, number_species)


	create_local_coordinate_file = 'cat {}.min{}.used_sites.tsv | grep -v POS | awk "{{print $2}}" > {}_coordinates'.format(scaffold, number_species, scaffold)

	create_global_coordinate_file = '''cat $coordinate_file | grep -v POS | awk '{ print $1 ":" $2 }' > coordinates'''
	
	os.system(create_local_coordinate_file)
	subprocess.call(create_global_coordinate_file, shell=True);

	remove_local_coordinate_file = "rm {}_coordinates".format(scaffold)
	os.system(remove_local_coordinate_file)

def main():
	scaffold_list = get_scaffold_list()
	
	Parallel(n_jobs=num_jobs)(delayed(subset_vcf_by_scaffold)(scaffold) for scaffold in scaffold_list)

	Parallel(n_jobs=num_jobs)(delayed(convert_vcf_to_nexus)(scaffold) for scaffold in scaffold_list)
	
	Parallel(n_jobs=num_jobs)(delayed(reformat_coordinate_files)(scaffold) for scaffold in scaffold_list)

if __name__ == "__main__":
        main()
