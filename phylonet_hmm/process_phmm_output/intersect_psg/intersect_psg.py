import csv 

# CSV file of 6,520 single copy orthologs from Kober and Pogson. Rows 5 - 1,012 contain the 1,008 genes with significant tests for positive selection
psg_file = "psg.csv"

# CSV file of the genes overlapping with introgression tracts 
introgressed_genes = "intersect.csv"

def find_psg_overlap():
	with open(psg_file,"r") as file1, open(introgressed_genes,"r") as file2:
		csv1 = csv.reader(file1)
		csv2 = csv.reader(file2)
		psg_list = list(csv1)
		introgressed_genes_list = list(csv2)

	# Get list of SPU identifiers from the set of positively selected genes 
	psg_list = [n for n in psg_list[4:1012]]

	for gene in psg_list:
		for record in introgressed_genes_list:
			if gene[0] in record or gene[1] in record or gene[2] in record:
				print(gene)
				print(record[5])
				print(record[6])
				print(record)
				print("Positively Selected Gene #{}".format(psg_list.index(gene)))
				print(" ")

def main():
	find_psg_overlap()

if __name__ == "__main__":
	main()
