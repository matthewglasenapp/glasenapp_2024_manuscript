import os
import csv

parent_directory = "/hb/home/mglasena/local_ancestry_analysis/paml_sco/single_copy_ortholog_alignments/"

dN_dS_dict = dict()

def get_gene_paths():
	subdirectory_lst = os.listdir(parent_directory)
	gene_paths_lst = [parent_directory + subdir for subdir in subdirectory_lst]
	return gene_paths_lst

def extract_paml():
	if("." in gene):
		return()
  	
  	with open(gene+'/paml_out.txt') as f:
		lines = f.readlines() # list containing lines of file
		columns = [] # To store column names
		i = 1
		count=0
		
		for line in lines:
	  		line = line.strip()
	  		
	  		if "omega (dN/dS) =  " in line:
				omega=float(line.split()[-1])
				count+=1
	 		
	 		if "tree length for dN:" in line:
				dn=float(line.split()[-1])
				count+=1
	  	
	  		if "tree length for dS:" in line:
				ds=float(line.split()[-1])
				count+=1
	 
	 		if count == 3:
				dN_dS_dict[gene.split('/')[-1]]=[ds,dn,omega]
				count = 0
				break
	f.close()

def write_data_to_csv():
	output_tsv_file = open("dNdS.tsv","w")
	writer = csv.writer(output_tsv_file, delimiter="\t")    
	header = ["gene", "dS", "dN", "dN/dS"]
	writer.writerow(header)

	for key,value in dN_dS_dict.items():
		data_to_write = [key, value[0], value[1], value[2]]
		writer.writerow(data_to_write)

	output_tsv_file.close()

def main():
	gene_paths_list = get_gene_paths()

	for gene in gene_paths_list:
		extract_paml()

	write_data_to_csv()

if __name__ == "__main__":
	main()
