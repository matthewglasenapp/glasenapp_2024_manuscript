import numpy as np
import os
import json

hmm_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm/"

json_output_dir = "/hb/scratch/mglasena/phylonet_hmm/probability_files/"

scaffold_output_file_dict = {}

scaffold_introgression_probabilities_dict = {}

final_probability_lst = {}

def get_all_bestrun_files():
	get_all_bestrun_files = 'find {} -type f -name "rawOutput.json" | grep "bestrun" > out_files'.format(hmm_dir)
	os.system(get_all_bestrun_files)

def get_bestrun_files_by_scaffold():
	with open("out_files","r") as f:
		files = f.read().splitlines()

	for line in files:
		if line.split("/")[7] in scaffold_output_file_dict.keys():
			scaffold_output_file_dict[line.split("/")[7]].append(line)
		else:
			scaffold_output_file_dict[line.split("/")[7]] = [line]

def get_introgression_probabilities_by_scaffold(scaffold):
		scaffold_introgression_probabilities_dict[scaffold] = []
		for file in scaffold_output_file_dict[scaffold]:
			with open(file,"r") as probability_file:
				scaffold_introgression_probabilities_dict[scaffold].append(json.load(probability_file).get("posteriorProbabilityOfSpeciesTrees")[0])

def average_probabilities_across_runs_by_scaffold(scaffold):
	arrays = [np.array(probability_lst) for probability_lst in scaffold_introgression_probabilities_dict[scaffold]]
	final_probability_lst[scaffold] = [np.mean(k) for k in zip(*arrays)]

def write_intorgression_probabilities_to_json(scaffold):
	introgression_probabilities_lst = final_probability_lst[scaffold]
	json_object = json.dumps(introgression_probabilities_lst)
	file_name = "{}{}.json".format(json_output_dir,scaffold)
	with open(file_name,"w") as outfile:
		outfile.write(json_object)

def main():
	get_all_bestrun_files()
	
	print("get_bestrun_files_by_scaffold")
	get_bestrun_files_by_scaffold()
	
	for scaffold in scaffold_output_file_dict:
		get_introgression_probabilities_by_scaffold(scaffold)

	print("Average probabilities across runs")
	for scaffold in scaffold_introgression_probabilities_dict:
		average_probabilities_across_runs_by_scaffold(scaffold)

	for scaffold in final_probability_lst:
		write_intorgression_probabilities_to_json(scaffold)

if __name__ == "__main__":
	main()