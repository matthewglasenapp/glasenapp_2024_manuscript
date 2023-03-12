import csv 

overlap_lst = []

csv_file = "psg.csv"
csv_file_2 = "intersect.csv"

with open(csv_file,"r") as file1, open(csv_file_2,"r") as file2:
	csv1 = csv.reader(file1)
	csv2 = csv.reader(file2)
	inputs = list(csv1)
	inputs2 = list(csv2)

gene_id_lst = [n[0] for n in inputs[4:1012]]

for item in gene_id_lst:
	for record in inputs2:
		if item in record[7]:
			if float(record[5]) == 100.0:
				print(item)
				print(record[1])
				print(record[2])
				print(gene_id_lst.index(item))
				print(" ")

