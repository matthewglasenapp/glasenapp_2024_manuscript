import os 
from joblib import Parallel, delayed
import csv

num_cores = 24

sample_names = {
'QB3KMK013': 'fragilis',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
'QB3KMK016': 'pulcherrimus',
}

outgroup = "pulcherrimus"

output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/ten_kb_tracts_gene_trees/"
make_output_dir = "mkdir -p {}".format(output_dir)
os.system(make_output_dir)

vcf2phylip_path = "/hb/home/mglasena/software/vcf2phylip/"

genotype_file = "/hb/scratch/mglasena/phylonet_hmm/phylonet_hmm_variant_sites_vcf/3bp_filtered_genotype_calls_pf.g.vcf.gz"

introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"

tract_info_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts/tract_info.tsv"

# nw_utils directory
nw_utils = "/hb/home/mglasena/software/newick_utils/src/"

def get_list_of_regions():
	regions = open(introgression_tract_file,"r").read().splitlines()
	regions = [region.split("\t")[3].replace("_","-").replace("-","_",1) for region in regions]
	return regions

def subset_vcf(region):
	# vcf files are 1 indexed
	one_indexed_region = region.split(":")[0] + ":" + str(int(region.split("-")[0].split(":")[1])+1) + "-" + region.split("-")[1]
	print(region)
	print(one_indexed_region)

	output_directory = output_dir + region
	os.system("mkdir " + output_directory)
	output_file = output_directory + "/" + region + ".vcf.gz"
	view = "bcftools view -r {} -o {} {}".format(one_indexed_region, output_file, genotype_file)
	print(view)
	os.system(view)

def create_fasta_alignments(region):
	input_file = output_dir + region + "/" + region + ".vcf.gz"
	output_directory = output_dir + region
	run_vcf2phylip = "python3 {}vcf2phylip.py -w --input {} -p -f --output-folder {} --output-prefix {}".format(vcf2phylip_path, input_file, output_directory, region)
	os.system(run_vcf2phylip)

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

def get_fasta_alignment_paths():
	find_fasta_files = 'find $PWD -type f -name *.fasta* > fasta_files'
	os.system(find_fasta_files)
	fasta_files = open("fasta_files","r").read().splitlines()
	os.system("rm fasta_files")
	return fasta_files

# Reconstruct gene trees for each fasta alignment in vcf2fasta_gene/ directory
def run_iqtree(fasta_file):
	run_iqtree = "iqtree2 -s {} -m MFP -B 1000".format(fasta_file)
	os.system(run_iqtree)

# Delete unwanted files 
def clean_up_iqtree_files():
	delete_gz = 'find ./ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./ -type f -name "*.bionj" -delete'
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

def get_tree_file_path_list():
	get_tree_file_paths = 'find $PWD -type f -name *.treefile* > tree_file_paths'
	os.system(get_tree_file_paths)
	tree_file_paths = open("tree_file_paths", "r").read().splitlines()
	os.system("rm tree_file_paths")
	return tree_file_paths

# Replace sample id numbers with species scientific names 
def edit_tree_files(input_file):
	with open(input_file, "r") as f:
		tree_list = f.read().splitlines()
	
	edited_tree_file = input_file.split(".min4")[0] + ".nwk"
	with open(edited_tree_file,"a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

	rooted_file = edited_tree_file.split(".nwk")[0] + "_rooted.nwk"
	root_and_prune = "nw_reroot {} {} -s | nw_prune - {} > {}".format(edited_tree_file, outgroup, outgroup, rooted_file)
	os.system(root_and_prune)

	ordered_rooted_file = rooted_file.split(".nwk")[0] + "_bl_removed.nwk"
	order_rooted_file = "nw_topology {} | nw_order - > {}".format(rooted_file, ordered_rooted_file)
	os.system(order_rooted_file)

def update_tract_info_csv_file():
	get_rooted_tree_file_lst = 'find $PWD -type f -name *rooted.nwk* > rooted_tree_file_lst.txt'
	get_rooted_bl_removed_tree_file_lst = 'find $PWD -type f -name *rooted_bl_removed.nwk* > rooted_bl_removed_tree_file_lst'
	os.system(get_rooted_tree_file_lst)
	os.system(get_rooted_bl_removed_tree_file_lst)
	rooted_tree_file_lst = open("rooted_tree_file_lst.txt","r").read().splitlines()
	rooted_bl_removed_tree_file_lst = open("rooted_bl_removed_tree_file_lst","r").read().splitlines()
	os.system("rm rooted_tree_file_lst.txt")
	os.system("rm rooted_bl_removed_tree_file_lst")

	tract_info_lines = open(tract_info_file,"r").read().splitlines()[1:]
	
	updated_csv = tract_info_file.split("tract_info.tsv")[0] + "updated_tract_info.tsv"
	output_csv = open(updated_csv,"w")
	writer = csv.writer(output_csv, delimiter = "\t")		

	header = ["tract_name", "Scaffold", "Start", "Stop", "Length", "SNV Sites", "SNV/bp", "Overlapping Coding Bases", "Percent Coding", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Genes", "Sdro 1x", "Sdro 10x", "Sfra 1x", "Sfra 10x", "Spal 1x", "Spal 10x", "Hpul 1x", "Hpul 10x", "Gene Tree", "Sdro-Spla sister?", "Bootstrap Support"]
	writer.writerow(header)

	counter = 0 
	for line in tract_info_lines:
		counter +=1 
		print(counter)
		length = int(line.split("\t")[4])
		if length >= 10000:
			tract = line.split("\t")[0].replace("-", "_", 1)
			gene_tree_file = output_dir + tract + "/" + tract + "_rooted.nwk"
			
			if gene_tree_file in rooted_tree_file_lst:
				gene_tree = open(gene_tree_file,"r").read().split("\n")[0]
				topology_file = output_dir + tract + "/" + tract + "_rooted_bl_removed.nwk"
				topology = open(topology_file,"r").read().split("\n")[0]
				substring = "(droebachiensis,pallidus)"
				if substring in topology:
					dro_pal_sister = "Yes"
					bootstrap = topology.split("(droebachiensis,pallidus)")[1].split(",")[0]
				else:
					dro_pal_sister = "No"
					bootstrap = "n/a"
			else:
				gene_tree = ""
				dro_pal_sister = ""
				bootstrap = ""

		else:
			gene_tree = ""
			dro_pal_sister = ""
			bootstrap = ""
		
		data = line.split("\t")
		data.append(gene_tree)
		data.append(dro_pal_sister)
		data.append(bootstrap)
		writer.writerow(data)
	
	output_csv.close()

def main():
	os.chdir(output_dir)

	regions = get_list_of_regions()

	Parallel(n_jobs=num_cores)(delayed(subset_vcf)(region) for region in regions)
	
	Parallel(n_jobs=num_cores)(delayed(create_fasta_alignments)(region) for region in regions)

	fasta_files = get_fasta_alignment_paths()

	Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_files)

	clean_up_iqtree_files()

	tree_files = get_tree_file_path_list()

	Parallel(n_jobs=num_cores)(delayed(edit_tree_files)(tree_file) for tree_file in tree_files)

	update_tract_info_csv_file()

if __name__ == "__main__":
	main()
