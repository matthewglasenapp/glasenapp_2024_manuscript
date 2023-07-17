import os 
from joblib import Parallel, delayed
import multiprocessing

# Specify number of cores accessible
num_cores = 24

# Change this to parent directory containing the gene directories with alignments and trees 
parent_directory = "/hb/home/mglasena/local_ancestry_project/paml_sco/single_copy_ortholog_fasta_alignments/"

def get_gene_paths():
	subdirectory_lst = os.listdir(parent_directory)
	gene_paths_lst = [parent_directory + subdir for subdir in subdirectory_lst]
	return gene_paths_lst

def run_paml(gene_path):
	os.chdir(gene_path)
	print("Running paml")
	os.system("codeml")

def main():
	gene_paths_list = get_gene_paths()
	print("Gene paths length: {}".format(len(gene_paths_list)))
	Parallel(n_jobs=num_cores)(delayed(run_paml)(gene) for gene in gene_paths_list)

if __name__ == "__main__":
	main()
