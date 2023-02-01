#!/usr/bin/env python

# Add function to build phylogenetic tree from the vcf2phylip nexus matrix that was the original input to phylonet-hmm.
# Use index_tract_list to index into the nexus matrices and pull out alignments to be run through iqtree. How do these compare to the other alignments?

import os 
import subprocess
import csv

working_directory = "/Users/matt/desktop/investigate_tracts/"
tract_file = "/Users/matt/desktop/investigate_tracts/ten_kb_tracts.bed"
gene_list = "/Users/matt/desktop/investigate_tracts/protein_coding_genes.bed"

# Organization!!! Figure out where I want these files to be located. What is the file structure I want??

def intersect_genes(tract_file,gene_list):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_list + " -wo > overlap")

# Build dictionary of gene identifiers, names, and curation status. {ID:[name,curation_stats]}
def get_gene_identifier_dict():
	#os.system("wget http://ftp.echinobase.org/pub/GenePageReports/GenePageGeneralInfo_AllGenes.txt")
	gene_page_download_dir = "/Users/matt/desktop/process_hmm/"
	os.chdir(gene_page_download_dir)
	gene_list = open("GenePageGeneralInfo_AllGenes.txt","r").read().splitlines()
	gene_dict = {}
	# Exclude records that have incorrect formatting
	new_list = [n for n in gene_list if len(n.split("\t")) >= 3]
	# Record excluded records 
	excluded_records = [n for n in gene_list if len(n.split("\t")) <3]
	#[print(n + "\n") for n in excluded_records]

	for n in new_list:
		split_list = n.split("\t",3)
		echinobase_gene_id=split_list[0]
		symbol = split_list[1]
		name = split_list[2]
		info = split_list[3].replace("\t"," ")

		gene_dict[symbol] = [name,info]

	return gene_dict

# Get file of introgression_tract\tgene_id\tgene_name\tbase_pair_overlap
def create_intersection_file():
	csv_file = open("/Users/matt/desktop/investigate_tracts/intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "overlapping_gene", "gene_name", "overlapping_bases", "gene_length", "percent_introgressed", "gene_coordinates", "gene_info"]
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

		gene_coordinates = str(record[4]) + ":" + str(record[5]) + "-" + str(record[6])
		gene_length = int(record[6]) - int(record[5])
		bp_overlap = record[14]
		percent_introgressed = (int(bp_overlap) / gene_length)*100
		data = [tract_id,gene_id,gene_name,gene_length,bp_overlap,percent_introgressed,gene_coordinates,gene_info]
		writer.writerow(data)
	csv_file.close()

def main():
	os.chdir(working_directory)

	intersect_genes(tract_file,gene_list)

	create_intersection_file()

if __name__ == "__main__":
	main()
