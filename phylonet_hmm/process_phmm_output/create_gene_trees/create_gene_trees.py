import os
from joblib import Parallel, delayed
num_cores = 24

# Mapping of DNA sample names to species names 
sample_names = {
'QB3KMK013': 'fragilis',
'QB3KMK010': 'franciscanus',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
'QB3KMK016': 'pulcherrimus',
'QB3KMK012': 'intermedius',
'SPUR.00': 'purpuratus',
}

# nw_utils directory
nw_utils = "/hb/groups/pogson_group/dissertation/software/newick_utils/src/"

# Path to vcf2fasta.py
vcf2fasta = "/hb/groups/pogson_group/dissertation/software/vcf2fasta/vcf2fasta.py"

# Feature of gff file that vcf2fasta.py will build alignments for
feature = "gene"

# Path to S. purpuratus reference genome
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# S. purpuratus gff3 file
gff_file = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Path to filtered multisample vcf file
vcf_file = "/hb/scratch/mglasena/data/genotypes/franciscanus/insertions_removed.g.vcf.gz"

#with open("intersect.csv","r") as f:
	#inputs = f.read().splitlines()
#gene_list = [item.split(",")[5] for item in inputs]

gene_lst = ["LOC115921720","LOC591845","LOC586866","LOC115921376","LOC586606","LOC115917953"]

# Using S. purpuratus gff3 file, make gff file for each gene that passed previous filters
def make_sco_gff(gene):
	command = "grep {} {} > single_gene_gff_records/{}.record".format(gene, gff_file, gene)
	os.system(command)

def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff sco_gff.gff --feat {}".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

def run_iqtree(fasta_file):
	run_iqtree = "iqtree2 -s vcf2fasta_gene/{} -m MFP -B 1000".format(fasta_file)
	os.system(run_iqtree)

def clean_up_iqtree_files():
	delete_gz = 'find ./vcf2fasta_gene/ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./vcf2fasta_gene/ -type f -name "*.bionj" -delete'
	delete_contree = 'find ./vcf2fasta_gene/ -type f -name "*.contree" -delete'
	delete_iqtree = 'find ./vcf2fasta_gene/ -type f -name "*.iqtree" -delete'
	delete_log = 'find ./vcf2fasta_gene/ -type f -name "*.log" -delete'
	delete_mldist = 'find ./vcf2fasta_gene/ -type f -name "*.mldist" -delete'
	os.system(delete_gz)
	os.system(delete_bionj)
	os.system(delete_contree)
	os.system(delete_iqtree)
	os.system(delete_log)
	os.system(delete_mldist)

def edit_tree_files(input_file):
	with open("vcf2fasta_gene/" + input_file, "r") as f:
		tree_list = f.read().splitlines()
	
	output_file = input_file.split(".fas")[0] + ".nwk"
	with open(output_file,"a") as f2:
		for tree in tree_list:
			for sample_name in sample_names.keys():
				if sample_name in tree:
					new_tree = tree.replace(sample_name, sample_names[sample_name])
					tree = new_tree
			f2.write(tree + "\n")

#For topology parsing, remove branch lengths and bootstraps. Arrange identical topologies to have the same texual representation. 
# nw_topology creates cladogram. Option -I gets rid of branc lengths. 
# nw_order orders the tree so that trees with identical topologies will have identical newick strings. Ooption -c d reorders the tree in such a way as to remove the ladder. 
def clean_gene_trees(input_file, output_file):
	clean = "{}nw_topology -I {} | {}nw_order -c d - > {}".format(nw_utils, input_file, nw_utils, output_file)
	os.system(clean)

def main():
	os.system("mkdir single_gene_gff_records/")
	Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_lst)

	# Concatenate all single gene gff records into "sco_gff.gff" file
	os.system('find ./single_gene_gff_records/ -type f -name "*.record" -exec cat {} \\; > sco_gff.gff')
	
	# Delete the single gene records
	os.system('find ./single_gene_gff_records/ -type f -name "*.record" -delete')
	os.system('rmdir single_gene_gff_records/')

	run_vcf2fasta()
	replace_missing_genotype_char()
	
	fasta_file_list = os.listdir("vcf2fasta_gene/")
	Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_file_list)
	
	tree_file_lst = [item for item in os.listdir("vcf2fasta_gene/") if "treefile" in item]
	Parallel(n_jobs=num_cores)(delayed(edit_tree_files)(input_file) for input_file in tree_file_lst)

if __name__ == "__main__":
	main()