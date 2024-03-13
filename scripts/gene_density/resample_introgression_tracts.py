import os
import random

cores = 24

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/gene_density/introgression_tracts/"

number_replicates = 1000

introgression_tract_file = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm/ten_kb_tracts.bed"

introgression_tract_list = open(introgression_tract_file,"r").read().splitlines()

# Directory for replicate interval files 
replicate_dir = root_dir + "replicate_interval_files/"
make_replicate_dir = "mkdir -p {}".format(replicate_dir)
os.system(make_replicate_dir)

def get_replicate_file_lst():
	replicate_file_lst = ["replicate_" + str(i) for i in range(1,number_replicates + 1)]
	return replicate_file_lst

def get_random_tract():
	random_number = random.randint(0, len(introgression_tract_list) - 1)
	random_tract = introgression_tract_list[random_number]
	return random_tract

def create_pseudoreplicate(output_file):
	resampled_tract_lst = []
	for i in range(1, len(introgression_tract_list) + 1):
		resampled_tract_lst.append(get_random_tract())
	
	with open(output_file + ".bed", "w") as f:
		for tract in resampled_tract_lst:
			f.write(tract + "\n")
	
def main():
	os.chdir(replicate_dir)
	replicate_file_lst = get_replicate_file_lst()
	for replicate in replicate_file_lst:
		create_pseudoreplicate(replicate)

if __name__ == "__main__":
	main()