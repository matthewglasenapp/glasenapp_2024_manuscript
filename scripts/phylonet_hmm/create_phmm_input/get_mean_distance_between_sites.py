import os
from statistics import mean

root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments/"

gap_lst = []

def append_to_gap_lst(input_coordinate_file):
	coordinates = [int(item.split(":")[1]) for item in open(input_coordinate_file,"r").read().splitlines()]
	for i in range(1, len(coordinates)):
		gap_lst.append(coordinates[i] - coordinates[i-1])

def get_coordinate_files_list():
	find = "find $PWD -type f -name '*coordinates*' > coordinate_files_list"
	os.system(find)
	coordinate_files_list = open("coordinate_files_list","r").read().splitlines()
	os.system("rm coordinate_files_list")
	return coordinate_files_list

def main():
	os.chdir(root_dir)
	coordinate_files_list = get_coordinate_files_list()
	for coordinate_file in coordinate_files_list:
		append_to_gap_lst(coordinate_file)
	print("The average gap between variable SNV sites in the PhyloNet-HMM input alignments is {} base pairs".format(mean(gap_lst)))
	print("The largest gap between variable SNV sites in the PhyloNet-HMM input alignments is {} base pairs".format(max(gap_lst)))

if __name__ == "__main__":
	main()
