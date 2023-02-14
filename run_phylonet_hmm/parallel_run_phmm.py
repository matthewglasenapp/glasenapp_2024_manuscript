#!/usr/bin/env python

import os
import multiprocessing
from joblib import Parallel, delayed

phylonet_hmm_path = "/hb/groups/pogson_group/dissertation/software/phylonet/PHiMM.jar"
output_dir = "/hb/scratch/mglasena/phylonet_hmm/run_3/hmm/"
memory = 5

num_jobs = 21

scaffold_alignments = "scaffold_input_nexus_file_paths_file"

def run_hmm(scaffold):
	os.system("java -Xmx{}g -jar {} {}".format(memory,phylonet_hmm_path,scaffold))

def main():
	inputs = open(scaffold_alignments, "r").read().splitlines()
	os.chdir(output_dir)
	Parallel(n_jobs=num_jobs)(delayed(run_hmm)(scaffold) for scaffold in inputs)

if __name__ == "__main__":
        main()