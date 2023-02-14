#!/usr/bin/env python

import os
import multiprocessing
from joblib import Parallel, delayed

phylonet_hmm_path = "/hb/groups/pogson_group/dissertation/software/phylonet/PHiMM.jar"
input_nexus_file_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/hmm_nexus_files/"

hmm_output_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm/"
os.system("mkdir " + hmm_output_dir)

memory = 5
number_runs = 100
num_scaffolds = 21

def get_scaffold_file_lst():
	get_nexus_file_lst = "find {} -type f > scaffold_input_nexus_file_paths_file".format(input_nexus_file_dir)
	os.system(get_nexus_file_lst)

def create_run_lst(number_runs):
	run_lst = []
	for i in range(1,number_runs + 1):
		run_lst.append("run_" + str(i))
	return run_lst

def run_hmm(scaffold):
	os.system("java -Xmx{}g -jar {} {}".format(memory,phylonet_hmm_path,scaffold))

def main():
	get_scaffold_file_lst()
	inputs = open("scaffold_input_nexus_file_paths_file", "r").read().splitlines()
	run_lst = create_run_lst(number_runs)
	for item in run_lst:
		os.mkdir(hmm_output_dir + item)
		os.chdir(hmm_output_dir + item)
		Parallel(n_jobs=num_scaffolds)(delayed(run_hmm)(scaffold) for scaffold in inputs)

if __name__ == "__main__":
        main()