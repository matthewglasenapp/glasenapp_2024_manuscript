import os
output_dir = "/hb/scratch/mglasena/phylonet_hmm/bestrun_files/"

bestrun_file_list = open("out_files","r").read().splitlines()

scaffold_list = ['NW_022145605.1', 'NW_022145615.1', 'NW_022145595.1', 'NW_022145594.1', 'NW_022145614.1', 'NW_022145604.1', 'NW_022145606.1', 'NW_022145596.1', 'NW_022145597.1', 'NW_022145607.1', 'NW_022145599.1', 'NW_022145603.1', 'NW_022145613.1', 'NW_022145609.1', 'NW_022145612.1', 'NW_022145602.1', 'NW_022145598.1', 'NW_022145600.1', 'NW_022145610.1', 'NW_022145611.1', 'NW_022145601.1']

def copy(file):
	scaffold = file.split("/")[-3]
	run = file.split("/")[-4]
	copy = "cp {} {}{}/{}_rawOutput.json".format(file, output_dir, scaffold, run)
	os.system(copy)

def main():
	os.system("mkdir " + output_dir)
	
	for scaffold in scaffold_list:
		mkdir = "mkdir {}{}".format(output_dir, scaffold)
		os.system(mkdir)

	for file in bestrun_file_list:
		copy(file)

if __name__ == "__main__":
	main()

	