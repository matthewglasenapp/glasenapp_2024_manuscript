import os

root_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/nexus_alignments_all_sites"

gap_dict = dict()

gap_threshold = 100000

def get_coordinate_files():
	find_coordinate_files = "find {} -type f -name *coordinates* > coordinate_files".format(root_dir)
	os.system(find_coordinate_files)
	with open("coordinate_files","r") as f:
		scaffold_coordinate_file_lst = f.read().splitlines()
	os.system("rm coordinate_files")
	return scaffold_coordinate_file_lst

def find_gaps(scaffold_coordinate_file):
	with open(scaffold_coordinate_file,"r") as f:
		inputs = f.read().splitlines()

	scaffold_name = inputs[0].split(":")[0]
	gap_dict[scaffold_name] = []
	coordinates = [item.split(":")[1] for item in inputs]

	gap_lst = []
	gap_len_lst = []
	
	i = 0
	first_base = 1
	first_coordinate = coordinates[0]
	gap_size = int(first_coordinate) - first_base
	if gap_size >= gap_threshold:
		gap = "{}:{}-{}".format(scaffold_name, first_base, first_coordinate)
		gap_lst.append(gap)

	while i < len(coordinates) - 1:
		position_1 = str(int(coordinates[i]) - 1)
		position_2 = coordinates[i+1]
		gap_size = int(position_2) - int(position_1)
		if gap_size >= gap_threshold:
			gap = "{}:{}-{}".format(scaffold_name, position_1, position_2)
			gap_lst.append(gap)
		i = i+1

	for item in gap_lst:
		gap_dict[scaffold_name].append(item)

def get_gap_stats():
	gap_lst = []
	for key,value in gap_dict.items():
		for item in value:
			gap_length = int(item.split(":")[1].split("-")[1]) - int(item.split(":")[1].split("-")[0])
			gap_lst.append(gap_length)

	Mb_gap_lst = [item for item in gap_lst if item >= 1000000]

	print("There were {} gaps greater than 100kb in length".format(len(gap_lst)))
	print("The largest gap was {}.".format(max(gap_lst)))
	print("There were {} gaps greater than 1Mb in length".format(len(Mb_gap_lst)))

	with open("gaps_by_scaffold.tsv","a") as f:
		for key,value in gap_dict.items():
			scaffold = key
			number_gaps = len(value)
			f.write(scaffold + "\t" + str(number_gaps) + "\n") 

def write_bed_file():
	with open("100kb_gaps.bed","w") as f:
		for key,value in gap_dict.items():
			for gap in value:
				scaffold = gap.split(":")[0]
				start = gap.split(":")[1].split("-")[0]
				stop = gap.split(":")[1].split("-")[1]
				f.write(scaffold + "\t" + start + "\t" + stop + "\t" + gap + "\n")

def main():
	scaffold_coordinate_file_lst = get_coordinate_files()
	for coordinate_file in scaffold_coordinate_file_lst:
		find_gaps(coordinate_file)
	number_gaps = 0
	for key,value in gap_dict.items():
		for item in value:
			number_gaps += 1
	get_gap_stats()
	write_bed_file()

if __name__ == "__main__":
	main()