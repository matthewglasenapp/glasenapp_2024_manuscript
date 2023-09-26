import os
import csv
from joblib import Parallel, delayed
num_cores = 24

# Mapping of DNA sample names to species names 
sample_names = {
'QB3KMK013': 'fragilis',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
'QB3KMK016': 'pulcherrimus',
}

output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/create_gene_trees/"

# nw_utils directory
nw_utils = "/hb/home/mglasena/software/newick_utils/src/"

# Path to vcf2fasta.py
vcf2fasta = "/hb/home/mglasena/software/vcf2fasta/vcf2fasta.py"

# Feature of gff file that vcf2fasta.py will build alignments for
feature = "gene"

# Path to S. purpuratus reference genome
reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# S. purpuratus gff3 file
gff_file = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Path to filtered multisample vcf file
vcf_file = "/hb/scratch/mglasena/phylonet_hmm/phylonet_hmm_variant_sites_vcf/filtered_genotype_calls.g.vcf.gz"

# Intersection file containing introgressed genes
intersect_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts/intersect.tsv"

# PSG intersection file
psg_intersect_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts/updated_psg_intersect.tsv"

# Get list of introgressed gene NCBI LOC IDs
def get_gene_list():
	introgressed_genes = [gene.split("\t")[0] for gene in open(intersect_file,"r").read().splitlines()[1:]]
	return introgressed_genes 

# Using S. purpuratus gff3 file, make gff file for each gene that passed previous filters
def make_sco_gff(gene):
	if "/" in gene:
		gene = gene.replace("/","_")
	command = "grep {} {} > single_gene_gff_records/{}.record".format(gene, gff_file, gene)
	os.system(command)

# Create fasta matrix for each gene in intersect_gff.gff file  
def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff intersect_gff.gff --feat {}".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_gene/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_gene/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

# Reconstruct gene trees for each fasta alignment in vcf2fasta_gene/ directory
def run_iqtree(fasta_file):
	run_iqtree = "iqtree2 -s vcf2fasta_gene/{} -m MFP -B 1000".format(fasta_file)
	os.system(run_iqtree)

# Delete unwanted files 
def clean_up_iqtree_files():
	delete_gz = 'find ./vcf2fasta_gene/ -type f -name "*.gz" -delete'
	delete_bionj = 'find ./vcf2fasta_gene/ -type f -name "*.bionj" -delete'
	delete_contree = 'find ./vcf2fasta_gene/ -type f -name "*.contree" -delete'
	delete_iqtree = 'find ./vcf2fasta_gene/ -type f -name "*.iqtree" -delete'
	delete_log = 'find ./vcf2fasta_gene/ -type f -name "*.log" -delete'
	delete_mldist = 'find ./vcf2fasta_gene/ -type f -name "*.mldist" -delete'
	delete_nex = 'find ./vcf2fasta_gene/ -type f -name "*.nex" -delete'
	delete_phy = 'find ./vcf2fasta_gene/ -type f -name "*.phy" -delete'
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

# Write new "gene_intersect_file.tsv" file that includes the gene tree topology for each introgressed gene along with all the existing info in intersect_file
def update_gene_intersection_file():
	csv_file = open("gene_intersect_file.tsv","w")
	writer = csv.writer(csv_file, delimiter="\t")
	header = ["NCBI Gene ID", "Name", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Introgression Tract(s)", "ECB Gene ID", "Gene Synonyms", "Curation Status", "Info", "Echinobase GO", "Echinobase KO", "Tu GO", "Gene Tree"]

	writer.writerow(header)

	for record in open(intersect_file).read().splitlines()[1:]:
		data = record.split("\t")
		if data[0] + ".nwk" in os.listdir():
			gene_tree = open(str(data[0]) + ".nwk", "r").readline()
			data.append(gene_tree)
		else:
			data.append("missing")

		writer.writerow(data)
	
	csv_file.close()

# Write new "psg_intersection_file.tsv" file that includes the gene tree topology for each introgressed, positively selected gene along with all the existing info in psg_intersect_file
def update_psg_intersection_file():
	csv_file = open("psg_intersection_file.tsv","w")
	writer = csv.writer(csv_file, delimiter="\t")
	header = ["NCBI Gene ID", "Name", "Synonyms", "Kober and Pogson Gene ID", "Kober and Pogson Name", "Kober and Pogson Synonyms", "PSG #", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Overlapping introgression tract(s)", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Coding Bases", "Percent Coding", "Gene Tree"]
	writer.writerow(header)

	for record in open(psg_intersect_file).read().splitlines()[1:]:
		data = record.split("\t")
		if data[0] + ".nwk" in os.listdir():
			gene_tree = open(str(data[0]) + ".nwk", "r").readline()
			data.append(gene_tree)
		else:
			data.append("missing")

		writer.writerow(data)
	
	csv_file.close()

# Concatenate all gene tree topologies into one file called "concatenated_trees.nwk"
def concatenate_trees():
	concatenate_trees = "cat *.nwk > concatenated_trees.nwk"
	os.system(concatenate_trees)

# For topology parsing, remove branch lengths and bootstraps. Arrange identical topologies to have the same texual representation. 
# nw_topology creates cladogram. Option -I gets rid of branch lengths. 
# nw_order orders the tree so that trees with identical topologies will have identical newick strings. Ooption -c d reorders the tree in such a way as to remove the ladder. 
def clean_gene_trees(input_file, output_file):
	outgroup = "franciscanus"
	clean = "{}nw_reroot {} {} | {}nw_topology -I - | {}nw_order -c d - > {}".format(nw_utils, input_file, outgroup, nw_utils, nw_utils, output_file)
	os.system(clean)

def main():
	os.chdir(output_dir)

	gene_lst = get_gene_list()

	os.system("mkdir single_gene_gff_records/")
	Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_lst)

	# Concatenate all single gene gff records into "intersect_gff.gff" file
	os.system('find ./single_gene_gff_records/ -type f -name "*.record" -exec cat {} \\; > intersect_gff.gff')
	
	# Delete the single gene records
	os.system('find ./single_gene_gff_records/ -type f -name "*.record" -delete')
	os.system('rmdir single_gene_gff_records/')

	run_vcf2fasta()
	# replace_missing_genotype_char()
	
	# fasta_file_list = os.listdir("vcf2fasta_gene/")
	# Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_file_list)
	
	# tree_file_lst = [item for item in os.listdir("vcf2fasta_gene/") if "treefile" in item]
	# Parallel(n_jobs=num_cores)(delayed(edit_tree_files)(input_file) for input_file in tree_file_lst)

	# update_gene_intersection_file()
	# update_psg_intersection_file()

	# clean_up_iqtree_files()
	# concatenate_trees()
	# clean_gene_trees("concatenated_trees.nwk", "clean_single_locus_trees.nwk")

if __name__ == "__main__":
	main()