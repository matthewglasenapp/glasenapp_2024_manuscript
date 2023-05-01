import os 
import json
phylonet_hmm_alignment_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments"
probability_file_dir = "/hb/scratch/mglasena/phylonet_hmm/probability_files/"

coordinate_by_scaffold_dict = dict()
probability_by_scaffold_dict = dict()

def get_coordinate_file_paths():
	print("Indexing scaffold coordinate files from {}".format(phylonet_hmm_alignment_dir))
	find_files = "find {} -type f -name '*coordinates*' > coordinate_files".format(phylonet_hmm_alignment_dir)
	os.system(find_files)
	coordinate_file_list = open("coordinate_files","r").read().splitlines()
	os.system("rm coordinate_files")
	return coordinate_file_list

def create_coordinate_dict(coordinate_file_list):
	for coordinate_file in coordinate_file_list:
		scaffold = open(coordinate_file,"r").readline().split(":")[0]
		coordinate_by_scaffold_dict[scaffold] = [item.split(":")[1] for item in open(coordinate_file,"r").read().splitlines()]

def get_probability_file_paths():
	find_files = "find {} -type f -name '*.json' > probability_files".format(probability_file_dir)
	os.system(find_files)
	probability_file_list = open("probability_files","r").read().splitlines()
	os.system("rm probability_files")
	return probability_file_list

def create_probability_dict(probability_file_list):
	for probability_file in probability_file_list:
		scaffold = probability_file.split(".json")[0].split("/")[-1]
		probability_by_scaffold_dict[scaffold] = json.load(open(probability_file,"r"))

def get_tract_list():
	tract_lst = [[item.split("\t")[0] + ":" + str(int(item.split("\t")[1]) + 1), item.split("\t")[0] + ":" + item.split("\t")[2]] for item in open("tracts_pf.bed","r").read().splitlines()]
	return tract_lst

def check_tract_probabilities(tract):
	scaffold = tract[0].split(":")[0]
	start = tract[0].split(":")[1]
	stop = tract[1].split(":")[1]

	if start in coordinate_by_scaffold_dict[scaffold] and stop in coordinate_by_scaffold_dict[scaffold]:
		start_index = coordinate_by_scaffold_dict[scaffold].index(start)
		stop_index = coordinate_by_scaffold_dict[scaffold].index(stop)

		for i in range(start_index, stop_index + 1):
			print(probability_by_scaffold_dict[scaffold][i])
			if probability_by_scaffold_dict[scaffold][i] < 0.9:
				print("Tract {} has an error!".format(tract))
	else:
		print(tract)

def main():
	coordinate_files = get_coordinate_file_paths()
	create_coordinate_dict(coordinate_files)

	probability_files = get_probability_file_paths()
	create_probability_dict(probability_files)

	tract_lst = get_tract_list()
	for tract in tract_lst:
		check_tract_probabilities(tract)

if __name__ == "__main__":
	main()