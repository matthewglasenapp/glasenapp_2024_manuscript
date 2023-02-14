import os
import multiprocessing
from joblib import Parallel, delayed

memory = 5

phylonet_hmm_path = "/hb/groups/pogson_group/dissertation/software/phylonet/PHiMM.jar"

input_file_dir = "/hb/scratch/mglasena/phylonet_hmm/run_3/hmm_input/hmm_nexus_files/"

root_dir = "/hb/scratch/mglasena/phylonet_hmm/run_3/"
hmm_dir = root_dir + "hmm/"

def get_scaffold_input_nexus_file_path_list():
	with open("scaffold_input_nexus_file_paths_file", "r") as f:
		scaffold_input_nexus_file_path_list = f.read().splitlines()
	
	return scaffold_input_nexus_file_path_list

def run_hmm(scaffold):
	run_hmm = "java -Xmx{}g -jar {} {}".format(memory, phylonet_hmm_path, scaffold)
	os.system(run_hmm)

def main():
	array_id = os.environ["array_id"]
	print("Array ID: {}".format(array_id))
	scaffold_input_nexus_file_path_list = get_scaffold_input_nexus_file_path_list()
	scaffold = scaffold_input_nexus_file_path_list[int(array_id)]
	print("Scaffold: {}".format(scaffold))

	os.chdir(hmm_dir)
	run_hmm(scaffold)

if __name__ == "__main__":
        main()
