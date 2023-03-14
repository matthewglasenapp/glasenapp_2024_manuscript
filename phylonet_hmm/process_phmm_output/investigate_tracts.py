#!/usr/bin/env python

import os 
import subprocess
import csv

working_directory = "/Users/matt/desktop/investigate_tracts/"
tract_file = "/Users/matt/desktop/investigate_tracts/ten_kb_tracts.bed"
gene_list = "/Users/matt/desktop/investigate_tracts/protein_coding_genes.bed"

#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
gene_info_file = "GenePageGeneralInfo_AllGenes.txt"

#os.system("wget https://download.echinobase.org/pub/GenePageReports/GeneGoTerms.txt")
gene_go_terms = "/Users/matt/desktop/investigate_tracts/GeneGoTerms.txt"

def intersect_genes(tract_file,gene_list):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_list + " -wo > overlap")

# Build dictionary of gene identifiers, names, and curation status. {ID:[name,curation_stats]}
def get_gene_identifier_dict():
	# Initialize dictionary with the format {LOCID:name,info,go_terms)
	gene_dict = {}

	gene_list = open(gene_info_file,"r").read().splitlines()

	gene_go_terms_list = [gene.split("\t") for gene in open("GeneGoTerms.txt","r").read().splitlines()]
	
	# Exclude records that have incorrect formatting
	passed_records = [gene for gene in gene_list if len(gene.split("\t")) >= 3]
	
	# Record excluded records 
	excluded_records = [gene for gene in gene_list if len(gene.split("\t")) <3]
	print("Excluded Records:")
	[print(gene + "\n") for gene in excluded_records]

	counter = 0
	for gene in passed_records:
		echinobase_gene_id = gene.split("\t",3)[0]
		symbol = gene.split("\t",3)[1]
		name = gene.split("\t",3)[2]
		info = gene.split("\t",3)[3].replace("\t"," ")

		gene_dict[symbol] = [name,info]

	
	GO_dict = {}

	for item in gene_go_terms_list:
		GO_dict[item[2]] = item[3]

	for gene in gene_dict:
		if gene in GO_dict.keys():
			gene_dict[gene].append(GO_dict[gene])
			counter += 1
		else:
			gene_dict[gene].append("n/a")

	print(gene_dict)
	return gene_dict

# Get file of introgression_tract\tgene_id\tgene_name\tbase_pair_overlap
def create_intersection_file():
	csv_file = open("/Users/matt/desktop/investigate_tracts/intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "overlapping_gene", "gene_name", "gene_length", "overlapping_bases", "percent_introgressed", "gene_coordinates", "Echinobase GO Terms", "gene_info"]
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
						print(gene_id)
						print(gene_name)
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

		gene_coordinates = str(record[4]) + ":" + str(record[5]) + "-" + str(record[6])
		gene_length = int(record[6]) - int(record[5])
		bp_overlap = record[14]
		percent_introgressed = (int(bp_overlap) / gene_length)*100
		data = [tract_id,gene_id,gene_name,gene_length,bp_overlap,percent_introgressed,gene_coordinates,go_terms,gene_info]
		writer.writerow(data)
	
	csv_file.close()

def main():
	os.chdir(working_directory)

	get_gene_identifier_dict()

	intersect_genes(tract_file,gene_list)

	create_intersection_file()

if __name__ == "__main__":
	main()
