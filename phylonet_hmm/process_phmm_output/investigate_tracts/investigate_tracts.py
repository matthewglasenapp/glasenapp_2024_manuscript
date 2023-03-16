#!/usr/bin/env python

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

# Build dictionary of gene identifiers, names, and curation status. {ID:[name,curation_stats]}
def get_gene_identifier_dict():
	# Initialize dictionary with the format {LOC_ID : name, info, go_terms, ko_terms)
	gene_dict = dict()

	# Initialize dictionary with the format {LOC_ID : SPU_IDs}
	loc_spu_dict = dict()

	# Initialize GO term dictionary with the format {LOC_ID: GO Terms}
	GO_dict = dict()

	# Initialize KO term dictionary with the format {LOC_ID: GO Terms }
	KO_dict = {}

	# Initialize Tu et al. (2012) GO term dictionary with the format{SPU_ID : GO Terms}
	tu_go_dict = dict()

	gene_page_general_info_list = open(gene_info_file,"r").read().splitlines()

	gene_go_terms_list = [gene.split("\t") for gene in open(gene_go_terms,"r").read().splitlines()]

	gene_ko_terms_list = [gene.split("\t") for gene in open(gene_ko_terms,"r").read().splitlines()]

	tu_go_terms_list = [gene.split(",",1) for gene in open("tu_2012_go.csv","r").read().splitlines()[1:]]
	
	# Exclude records that have incorrect formatting
	passed_records = [gene for gene in gene_page_general_info_list if len(gene.split("\t")) >= 3]
	
	# Record excluded records 
	excluded_records = [gene for gene in gene_page_general_info_list if len(gene.split("\t")) <3]
	#print("Excluded Records:")
	#[print(gene + "\n") for gene in excluded_records]

	# Populate gene_dict
	for gene in passed_records:
		LOC_ID = gene.split("\t",3)[1]
		echinobase_gene_id = gene.split("\t",3)[0]
		name = gene.split("\t",3)[2]
		info = gene.split("\t",3)[3].replace("\t"," ")

		gene_dict[LOC_ID] = [name,info]

	# Populate loc_spu_dict 
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

	# Population GO_dict
	for item in gene_go_terms_list:
		GO_dict[item[2]] = item[3]

	# Populate KO_dict 
	for item in gene_ko_terms_list:
		KO_dict[item[2]] = item[3]

	# Population tu_go_dict
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

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file, gene_file, outfile):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_file + " -wo > " + outfile)

# Get file of introgression_tract\tgene_id\tgene_name\tbase_pair_overlap
def create_gene_intersection_file(overlap_file):
	csv_file = open("/Users/matt/desktop/investigate_tracts/intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "overlapping_gene", "gene_name", "gene_length", "overlapping_bases", "percent_introgressed", "gene_coordinates", "Echinobase GO", "Echinobase KO", "gene_info", "Tu GO"]
	writer.writerow(header)
	overlaps = open(overlap_file,"r").read().splitlines()

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
	get_gene_identifier_dict()

	intersect_genes(tract_file, gene_list_file, "gene_overlap.bed")
	intersect_genes(tract_file, exon_list_file, "exon_overlap.bed")

	create_gene_intersection_file("gene_overlap.bed")
	#create_exon_intersection_file("exon_overlap.bed")

if __name__ == "__main__":
	main()
