import os

#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
inputs = open("GenePageGeneralInfo_AllGenes.txt","r").read().splitlines()

record_lst = [record.split("\t") for record in inputs]
counter = 0 

for lst in record_lst[0:20]:
	if len(lst) == 6:

		print(lst)
		lst = list(filter(None, lst))
		print(lst)

		echinobase_gene_id = lst[0]
		print("Echinobase ID: {}".format(echinobase_gene_id))
		symbol = lst[1]
		print("Symbol: {}".format(symbol))
		name = lst[2]
		print("name: {}".format(name))

		del lst[0:3]

		for item in lst:
			if "curated" in item:
				curation_status = item
				print("Curation Status: {}".format(curation_status))
				lst.pop(lst.index(curation_status))

		print(lst)

