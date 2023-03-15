#!/usr/bin/env python

import os 
import subprocess
import csv

working_directory = "/Users/matt/desktop/investigate_tracts/"
tract_file = "/Users/matt/desktop/investigate_tracts/ten_kb_tracts.bed"
gene_list_file = "/Users/matt/desktop/investigate_tracts/protein_coding_genes.bed"

#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
gene_info_file = "GenePageGeneralInfo_AllGenes.txt"

#os.system("wget https://download.echinobase.org/pub/GenePageReports/GeneGoTerms.txt")
gene_go_terms = "/Users/matt/desktop/investigate_tracts/GeneGoTerms.txt"

#os.system("wget https://download.echinobase.org/pub/GenePageReports/GeneKoTerms.txt")
gene_ko_terms = "/Users/matt/desktop/investigate_tracts/GeneKoTerms.txt"

# Tu et al. 2012 S purpuratus GO terms 
tu_go_terms = "/Users/matt/desktop/investigate_tracts/tu_2012_go.csv"

def intersect_genes(tract_file,gene_list_file):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_list_file + " -wo > overlap")

# Build dictionary of gene identifiers, names, and curation status. {ID:[name,curation_stats]}
def get_gene_identifier_dict():
	# Initialize dictionary with the format {LOCID:name,info,go_terms,ko_terms)
	gene_dict = {}

	loc_spu_dict = dict()

	gene_list = open(gene_info_file,"r").read().splitlines()

	gene_go_terms_list = [gene.split("\t") for gene in open(gene_go_terms,"r").read().splitlines()]

	gene_ko_terms_list = [gene.split("\t") for gene in open(gene_ko_terms,"r").read().splitlines()]

	tu_go_terms_list = [gene.split(",",1) for gene in open("tu_2012_go.csv","r").read().splitlines()[1:]]
	
	# Exclude records that have incorrect formatting
	passed_records = [gene for gene in gene_list if len(gene.split("\t")) >= 3]
	
	# Record excluded records 
	excluded_records = [gene for gene in gene_list if len(gene.split("\t")) <3]
	#print("Excluded Records:")
	#[print(gene + "\n") for gene in excluded_records]

	for gene in passed_records:
		gene = gene.split("\t")
		LOC_ID = gene[1]
		if "SPU_" in gene[4]:
			spu_lst = []
			for item in gene[4].split("SPU_")[1:]:
				if len(item) >= 1:
					spu_number = item.split("|")[0]
					spu_lst.append("SPU_" + str(spu_number))
			loc_spu_dict[LOC_ID] = spu_lst

	for gene in passed_records:
		echinobase_gene_id = gene.split("\t",3)[0]
		symbol = gene.split("\t",3)[1]
		name = gene.split("\t",3)[2]
		info = gene.split("\t",3)[3].replace("\t"," ")

		gene_dict[symbol] = [name,info]

	GO_dict = {}

	for item in gene_go_terms_list:
		GO_dict[item[2]] = item[3]

	KO_dict = {}

	for item in gene_ko_terms_list:
		KO_dict[item[2]] = item[3]

	tu_go_dict = dict()

	for item in tu_go_terms_list:
		tu_go_dict[item[0]] = item[1]

	for gene in gene_dict:
		if gene in GO_dict.keys():
			gene_dict[gene].append(GO_dict[gene])
		else:
			gene_dict[gene].append("n/a")
		
		if gene in KO_dict.keys():
			gene_dict[gene].append(KO_dict[gene])
		else:
			gene_dict[gene].append("n/a")

		if gene in loc_spu_dict.keys():
			if len(loc_spu_dict[gene]) > 1:
				spu_lst = loc_spu_dict[gene]
				go_lst =[]
				for spu_number in spu_lst:
					if spu_number in tu_go_dict.keys():
						go_lst.append(tu_go_dict[spu_number])
				if len(go_lst) >= 1:
					go = ",".join(go_lst)
					gene_dict[gene].append(go)
				else:
					gene_dict[gene].append("n/a")
			
			else:
				spu = loc_spu_dict[gene][0]
				if spu in tu_go_dict.keys():
					gene_dict[gene].append(tu_go_dict[spu])
				else:
					gene_dict[gene].append("n/a")
		else:
			gene_dict[gene].append("n/a")

	return gene_dict

# Get file of introgression_tract\tgene_id\tgene_name\tbase_pair_overlap
def create_intersection_file():
	csv_file = open("/Users/matt/desktop/investigate_tracts/intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "overlapping_gene", "gene_name", "gene_length", "overlapping_bases", "percent_introgressed", "gene_coordinates", "Echinobase GO", "Echinobase KO", "gene_info", "Tu GO"]
	writer.writerow(header)
	overlaps = open("/Users/matt/desktop/investigate_tracts/overlap").read().splitlines()

	gene_dict = get_gene_identifier_dict()

	for record in overlaps:
		record = record.split("\t")
		tract_id = record[3].replace("_","-")
		gene_id = record[7].split("-",1)[1].split(";",1)[0]

		try:
			gene_name = gene_dict[gene_id][0]
		except KeyError:
			try:
				for key,value in gene_dict.items():
					if gene_id in value[1]:
						gene_name = key
			except KeyError:
				gene_name = "missing"

		try:
			gene_info = gene_dict[gene_id][1]
		except KeyError:
			try: 
				for key,value in gene_dict.items():
					if gene_id in value[1]:
						gene_info = value
			
			except:
				gene_info = "missing record"

		try:
			go_terms = gene_dict[gene_id][2]
		except KeyError:
			try:
				for key,value in gene_dict.items():
					if gene_id in value[1]:
						go_terms = value[2]
			except KeyError:
				go_terms = "missing"

		try:
			ko_terms = gene_dict[gene_id][3]
		except KeyError:
			try:
				for key,value in gene_dict.items():
					if gene_id in value[1]:
						ko_terms = value[3]
			except KeyError:
				ko_terms = "missing"

		try:
			tu_go_terms = gene_dict[gene_id][4]
		except KeyError:
			try:
				for key,value in gene_dict.items():
					if gene_id in value[1]:
						tu_go_terms = value[4]
			except KeyError:
				ko_terms = "missing"

		gene_coordinates = str(record[4]) + ":" + str(record[5]) + "-" + str(record[6])
		gene_length = int(record[6]) - int(record[5])
		bp_overlap = record[14]
		percent_introgressed = (int(bp_overlap) / gene_length)*100
		data = [tract_id,gene_id,gene_name,gene_length,bp_overlap,percent_introgressed,gene_coordinates,go_terms,ko_terms,gene_info,tu_go_terms]
		writer.writerow(data)
	
	csv_file.close()

def main():
	os.chdir(working_directory)

	get_gene_identifier_dict()

	intersect_genes(tract_file,gene_list_file)

	create_intersection_file()

if __name__ == "__main__":
	main()
