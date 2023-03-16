import os
import csv 

# Introgression tract file in bed format
tract_file = "ten_kb_tracts.bed"

# Bed file of all S. purpuratus protein coding genes 
gene_list_file = "protein_coding_genes.bed"

exon_list_file = "unique_exons.bed"

# Gene metadata from Echinobase 
#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
gene_info_file = "GenePageGeneralInfo_AllGenes.txt"

# Gene GO Terms from Echinobase 
#os.system("wget https://download.echinobase.org/pub/GenePageReports/GeneGoTerms.txt")
gene_go_terms = "GeneGoTerms.txt"

# Gene KO Terms from Echinobase
#os.system("wget https://download.echinobase.org/pub/GenePageReports/GeneKoTerms.txt")
gene_ko_terms = "GeneKoTerms.txt"

# Tu et al. 2012 S. purpuratus GO terms
# https://genome.cshlp.org/content/22/10/2079/suppl/DC1
# Supp table S2
tu_go_terms = "tu_2012_go.csv"

gene_dictionary = dict()

ECB_SPU_dict = dict()

gene_intersection_dict = dict()

# Create gene dictionary from GenePageGeneralInfo_AllGenes.txt and write to csv
def create_gene_dictionary():
	#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
	inputs = open(gene_info_file,"r").read().splitlines()

	record_lst = [record.split("\t") for record in inputs]

	for lst in record_lst:
		lst = list(filter(None, lst))
		echinobase_gene_id = lst[0]
		symbol = lst[1]
		name = lst[2]

		del lst[0:3]

		curation_status = "n/a"
		synonyms = "n/a"
		synonyms_2 = "n/a"
		synonyms_3 = "n/a"
		info = "n/a"

		for item in lst:
			if "curated" in item:
				curation_status = item
				lst.pop(lst.index(curation_status))

		if len(lst) == 1:
			synonyms = lst[0]
			lst.pop(0)

		elif len(lst) >= 1:
			for item in lst:			
				# I checked to make sure that "|" only appears in the list once. Won't have problems with overwriting synonyms_2
				if "|" in item:
					synonyms_2 = item
					lst.pop(lst.index(item))
			
		if len(lst) > 1:
			for item in lst:
				if "SPU_" in item:
					synonyms_3 = item
					lst.pop(lst.index(synonyms_3))

		if len(lst) == 1:
			info = lst[0].strip()
			lst.pop(0)

		if "|" in name and synonyms == "n/a":
			synonyms = name
			name = "missing"

		if len(name) > 200 and info == "n/a":
			info = name
			name = "missing"
	
		# There is only one case of synonyms_3 and its record does not have synonyms or synonyms_1

		if synonyms != "n/a" and synonyms_2 == "n/a" and synonyms_3 == "n/a":
			combined_synonyms = synonyms
	
		elif synonyms == "n/a" and synonyms_2 != "n/a" and synonyms_3 == "n/a":
			combined_synonyms = synonyms_2
	
		elif synonyms == "n/a" and synonyms_2 == "n/a" and synonyms_3 != "n/a":
			combined_synonyms = synonyms_3

		elif synonyms != "n/a" and synonyms_2 != "n/a":
			combined_synonyms = synonyms.strip() + "|" + synonyms_2.strip()

		else:
			combined_synonyms = "n/a"

		# Get rid of any "||". I will splitting the synonyms into a list by the "|" delimiter 
		if "||" in combined_synonyms:
			combined_synonyms = "|".join(combined_synonyms.split("||"))

		combined_synonyms = combined_synonyms.split("|")
		
		# Population Gene Dictionary with info
		if not symbol in gene_dictionary.keys():
			gene_dictionary[echinobase_gene_id] = [symbol, name, combined_synonyms, curation_status, info]
		
		# Deal with duplicate records 
		else:
			if gene_dictionary[echinobase_gene_id][3] == "curated" or gene_dictionary[echinobase_gene_id][3] == "manually curated":
				if curation_status != "curated" or curation_status != "manually curated":
					continue
			elif gene_dictionary[echinobase_gene_id][3] != "curated" or gene_dictionary[echinobase_gene_id][3] != "manually curated":
				if curation_status == "curated" or curation_status == "manually curated":
					gene_dictionary[echinobase_gene_id] = [symbol, name, combined_synonyms, curation_status, info]

				# Just stick with the first entry. Was more reliable in the only two edge cases
				else:
					continue

def create_ECB_SPU_mapping_dict():
	for key,value in gene_dictionary.items():
		if "SPU_" in gene_dictionary[key][1]:
			spu_id = gene_dictionary[key][1]
			ECB_SPU_dict[key] = [spu_id]
		
		else:
			spu_lst = []
			for item in gene_dictionary[key][2]:
				if "SPU_" in item:
					spu_lst.append(item)
				ECB_SPU_dict[key] = spu_lst

def add_GO_KO_termns_to_gene_dictionary():
	gene_go_terms_list = [gene.split("\t") for gene in open(gene_go_terms,"r").read().splitlines()]

	gene_ko_terms_list = [gene.split("\t") for gene in open(gene_ko_terms,"r").read().splitlines()]

	# File has one line header
	tu_go_terms_list = [gene.split(",",1) for gene in open("tu_2012_go.csv","r").read().splitlines()[1:]]

	# Add Echinobase GO Terms to gene dictionary
	for gene in gene_go_terms_list:
		# Account for duplicate records in GeneGoTerms.txt
		if len(gene_dictionary[gene[0]]) == 5:
			gene_dictionary[gene[0]].append([gene[3]])
		elif len(gene_dictionary[gene[0]]) == 6:
			gene_dictionary[gene[0]][5].append(gene[3])

	# Append "n\a" to gene dictionary records if didn't append GO Terms
	for key,value in gene_dictionary.items():
		if len(value) == 5:
			gene_dictionary[key].append("n/a")

	# Add Echinobase KO Terms to Gene Dictionary 
	for gene in gene_ko_terms_list:
		# Account for duplicate records in GeneKoTerms.txt
		if len(gene_dictionary[gene[0]]) == 6:
			gene_dictionary[gene[0]].append([gene[3]])
		elif len(gene_dictionary[gene[0]]) == 7:
			gene_dictionary[gene[0]][5].append(gene[3])

	# Append "n\a" to gene dictionary records if didn't append GO Terms
	for key,value in gene_dictionary.items():
		if len(value) == 6:
			gene_dictionary[key].append("n/a")

	# Add Tu et al. GO Terms 
	# First, create dictionary of format {SPU_number : GO_Terms}
	tu_go_dict = dict()
	for item in tu_go_terms_list:
		tu_go_dict[item[0]] = item[1]

	for gene in gene_dictionary:
		if len(ECB_SPU_dict[gene]) > 1:
			spu_lst = ECB_SPU_dict[gene]
			go_lst =[]
			for spu_number in spu_lst:
				if spu_number in tu_go_dict.keys():
					go_lst.append(tu_go_dict[spu_number])
			if len(go_lst) >= 1:
				go = ",".join(go_lst)
				gene_dictionary[gene].append(go)
			else:
				gene_dictionary[gene].append("n/a")
			
		elif len(ECB_SPU_dict[gene]) == 0:
			gene_dictionary[gene].append("n/a")
		
		else:
			spu = ECB_SPU_dict[gene][0]
			if spu in tu_go_dict.keys():
				gene_dictionary[gene].append(tu_go_dict[spu])
			else:
				gene_dictionary[gene].append("no_tu")

	# Eliminate redundant items and change format of GO/KO records from list to string
	for value in gene_dictionary.values():
		if type(value[5]) is list:
			value[5] = ",".join(list(set(value[5])))

		if type(value[6]) is list:
			value[6] = ",".join(list(set(value[6])))

def write_gene_dictionary_to_csv():
	csv_file = open("test.csv","w")
	writer = csv.writer(csv_file)	
	header = ["symbol", "echinobase_gene_id", "name", "synonyms", "curation status", "info", "ECB GO Terms", "ECB KO Terms", "Tu et al. GO Terms"]
	writer.writerow(header)

	# Write gene info to csv file
	for key,value in gene_dictionary.items():
		data = [value[0], key, value[1], "|".join(value[2]), value[3], value[4], value[5], value[6], value[7]]
		writer.writerow(data)

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file, gene_file, outfile):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_file + " -wo > " + outfile)

def create_gene_intersection_dict(overlap_file):
	overlaps = open(overlap_file,"r").read().splitlines()

	for record in overlaps[0:1]:
		record = record.split("\t")
		tract_id = record[3].replace("_","-")
		#LOC ID
		gene_id = record[7].split("-",1)[1].split(";",1)[0]

		#try:
			#gene_name = gene_dict[gene_id][0]
		#except KeyError:
			#try:
				#for key,value in gene_dict.items():
					#if gene_id in value[1]:
						#gene_name = key
			#except KeyError:
				#gene_name = "missing"

		#try:
			#gene_info = gene_dict[gene_id][1]
		#except KeyError:
			#try: 
				#for key,value in gene_dict.items():
					#if gene_id in value[1]:
						#gene_info = value
			
			#except:
				#gene_info = "missing record"

		#try:
			#go_terms = gene_dict[gene_id][2]
		#except KeyError:
			#try:
				#for key,value in gene_dict.items():
					#if gene_id in value[1]:
						#go_terms = value[2]
			#except KeyError:
				#go_terms = "missing"

		#try:
			#ko_terms = gene_dict[gene_id][3]
		#except KeyError:
			#try:
				#for key,value in gene_dict.items():
					#if gene_id in value[1]:
						#ko_terms = value[3]
			#except KeyError:
				#ko_terms = "missing"

		#try:
			#tu_go_terms = gene_dict[gene_id][4]
		#except KeyError:
			#try:
				#for key,value in gene_dict.items():
					#if gene_id in value[1]:
						#tu_go_terms = value[4]
			#except KeyError:
				#ko_terms = "missing"

		gene_coordinates = str(record[4]) + ":" + str(record[5]) + "-" + str(record[6])
		gene_length = int(record[6]) - int(record[5])
		bp_overlap = record[14]
		percent_introgressed = (int(bp_overlap) / gene_length)*100

def create_gene_intersection_file():
	csv_file = open("intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "overlapping_gene", "gene_name", "gene_length", "overlapping_bases", "percent_introgressed", "gene_coordinates", "Echinobase GO", "Echinobase KO", "gene_info", "Tu GO"]
	writer.writerow(header)
	
	for key,value in gene_intersection_dict.values():
		data = []
		writer.writerow(data)
	
	csv_file.close()

def main():
	create_gene_dictionary()
	create_ECB_SPU_mapping_dict()
	add_GO_KO_termns_to_gene_dictionary()
	write_gene_dictionary_to_csv()
	
	intersect_genes(tract_file, gene_list_file, "gene_overlap.bed")
	intersect_genes(tract_file, exon_list_file, "exon_overlap.bed")
	
	create_gene_intersection_dict("gene_overlap.bed")
	#create_gene_intersection_dict("exon_overlap.bed")

	#create_gene_intersection_file()
	#create_exon_intersection_file()

	
	# Identify the problematic duplicates:
	lst = []
	for key,value in gene_dictionary.items():
		lst.append(value[0])
	
	unique = set()
	duplicates = []
	for item in lst:
		if item in unique:
			duplicates.append(item)
		else:
			unique.add(item)

	print(duplicates)

	for key,value in gene_dictionary.items():
		for item in duplicates:
			if value[0] == item:
				print(key)



if __name__ == "__main__":
	main()

