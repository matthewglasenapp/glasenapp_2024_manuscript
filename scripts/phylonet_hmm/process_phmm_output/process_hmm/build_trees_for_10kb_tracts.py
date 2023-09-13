import os 
from joblib import Parallel, delayed

sample_names = {
'QB3KMK013': 'fragilis',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
'QB3KMK016': 'pulcherrimus',
}

num_cores = 24
output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/ten_kb_tracts_gene_trees/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

vcf2phylip_path = "/hb/home/mglasena/software/vcf2phylip/"
genotype_file = "/hb/home/mglasena/dissertation/data/genotypes/franciscanus/3bp_filtered_genotype_calls_pf.g.vcf.gz"
introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"

def get_list_of_regions():
	regions = open(introgression_tract_file,"r").read().splitlines()
	regions = [region.split("\t")[3].replace("_","-").replace("-","_",1) for region in regions]
	return regions

def subset_vcf(region):
	output_dir = output_dir + region
	os.system("mkdir " + output_dir)
	output_file = output_dir + "/" + region + ".vcf.gz"
	view = "bcftools view -r {} -o {} {}".format(region, output_file, genotype_file)
	os.system(view)

def create_fasta_alignments(region):
	input_file = output_dir + region + "/" + region + ".vcf.gz"
	output_dir = output_dir + region
	run_vcf2phylip = "python3 {}vcf2phylip.py -w --input {} -p -f --output-folder {} --output-prefix {}".format(vcf2phylip_path, input_file, output_dir, region)
	os.system(run_vcf2phylip)

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

def get_fasta_alignment_paths():
	fasta_files = os.listdir("vcf2fasta_CDS/")
	fasta_file_paths_lst = [root_dir + "vcf2fasta_CDS/" + item for item in fasta_files]
	return fasta_file_paths_lst

# Reconstruct gene trees for each fasta alignment in vcf2fasta_gene/ directory
def run_iqtree(fasta_file):
	run_iqtree = "iqtree2 -s vcf2fasta_gene/{} -m MFP -B 1000".format(fasta_file)
	os.system(run_iqtree)

# Delete unwanted files 
def clean_up_iqtree_files():
	delete_gz = 'find ./ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./-type f -name "*.bionj" -delete'
	delete_contree = 'find ./ -type f -name "*.contree" -delete'
	delete_iqtree = 'find ./ -type f -name "*.iqtree" -delete'
	delete_log = 'find ./ -type f -name "*.log" -delete'
	delete_mldist = 'find ./ -type f -name "*.mldist" -delete'
	delete_nex = 'find ./ -type f -name "*.nex" -delete'
	delete_phy = 'find ./ -type f -name "*.phy" -delete'
	os.system(delete_gz)
	os.system(delete_bionj)
	os.system(delete_contree)
	os.system(delete_iqtree)
	os.system(delete_log)
	os.system(delete_mldist)
	os.system(delete_nex)
	os.system(delete_phy)

# Replace sample id numbers with species scientific names 
def edit_tree_files(input_file):
	with open(input_file, "r") as f:
		tree_list = f.read().splitlines()
	
	output_file = input_file.split(".min4")[0] + ".nwk"
	with open(output_file,"a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

def main():
	os.chdir(output_dir)

	regions = get_list_of_regions()

	Parallel(n_jobs=num_cores)(delayed(subset_vcf)(region) for region in regions)
	
	Parallel(n_jobs=num_cores)(delayed(create_fasta_alignments)(region) for region in regions)

	fasta_files = 

	Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_files)

	clean_up_iqtree_files()

	tree_files = 

if __name__ == "__main__":
	main()
