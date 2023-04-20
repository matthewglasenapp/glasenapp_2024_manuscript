import os
import csv 

# Total bases introgressed
#prob_90
total_bases_introgressed = 7903066
#prob_80
#total_bases_introgressed = 35708711
#prob_75
#total_bases_introgressed = 50712387

# Introgression tract file in bed format
tract_file = "tracts.bed"

# Coverage depth of each tract for each species 
tract_coverage_file = "tract_coverage.tsv"

# Bed file of all S. purpuratus protein coding genes 
gene_list_file = "protein_coding_genes.bed"

# File of unique protein-coding exons 
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

# Gene dictionary in the format of {ECB-GENEPAGE_ID: symbol, name, synonyms, curation_status, info, ECB GO Terms, ECB KO Terms, Tu GO Terms}
gene_dictionary = dict()

# Dictionary of {ECCB-GENEPAGE_ID: [SPU_IDs]}
ECB_SPU_dict = dict()

# Same dictionary as gene_dictionary, but formatted as {LOC_ID: ECB-GENEPAGE-ID, name, synonyms, curation_stats, info, ECB GO Terms, ECB KO Terms, Tu GO Terms}
LOC_gene_dictionary = dict()

# Dictionary for genes overlapping with an introgression tract {LOC_ID: data}
gene_intersection_dict = dict()

# Dictionary of the average coverage depths of each introgression tract for each species
tract_coverage_depth_dict = dict()

# Create gene dictionary from GenePageGeneralInfo_AllGenes.txt
def create_gene_dictionary():
	inputs = open(gene_info_file,"r").read().splitlines()

	record_lst = [record.split("\t") for record in inputs]

	# Parse each record for info
	for lst in record_lst:
		
		# Delete any empty lists that occurred due to consecutive tab characters. Each record was split into a list on each "\t"
		lst = list(filter(None, lst))
		
		# Each record begins with an ECB-GENEPAGE- numer 
		echinobase_gene_id = lst[0]
		
		# Symbol is generally the corresponding LOC number, but sometimes it is a name
		symbol = lst[1]
		name = lst[2]

		# Delete echinobase_gene_id, symbol, and name form list after saving them to variables. This makes parsing the rest of the list easier. 
		del lst[0:3]

		# Define variables to be filled later.
		curation_status = "n/a"
		synonyms = "n/a"
		info = "n/a"

		# Populate curation_status variable and remove item from list. Each record should have a status of either "partially curated" or "curated"
		for item in lst:
			if "curated" in item:
				curation_status = item
				lst.pop(lst.index(curation_status))

		# If there is only one item left in the list, that item is the gene synonyms. Assign this item to synonyms variable
		if len(lst) == 1:
			synonyms = lst[0]
			lst.pop(0)

		# If there are multiple items left in the list, one of them is synonyms and the other is a general info column. 
		# Find the synonyms record and assign it to synonyms. Remove that item from the list 
		elif len(lst) > 1:
			for item in lst:			
				# I checked to make sure that "|" only appears in one item per list. If it appeared twice, synonyms would be overwritten
				if "|" in item:
					synonyms = item
					lst.pop(lst.index(item))
			
		# There is one weird edge case where the list still has two records: ['Metallothioneins have a high content of cysteine residues that bind various heavy metals.', 'SPU_017134']
		# I checked to ensure that this records synonyms variable is "n/a"
		# Assigned 'SPU_017134' to synonyms variable 
		if len(lst) > 1:
			for item in lst:
				if "SPU_" in item:
					synonyms = item
					lst.pop(lst.index(synonyms))

		# Now all records are a list of 1 or empty lists. The records with a list of one are the info column. Populate info variable.
		# After this chunk, all records will be an empty list 
		if len(lst) == 1:
			info = lst[0].strip()
			lst.pop(0)

		# Control for weird edge cases where what we stored as the name is actually a list of synonyms. Example: name = SPU_010482|Sp-PolypL_68|polyprotein-like-68
		if "|" in name and synonyms == "n/a":
			synonyms = name
			name = "missing"

		# Control for three weird edge case where what we stored as the name is actually the info
		if len(name) > 200 and info == "n/a":
			info = name
			name = "missing"

		# Get rid of any "||". I will splitting the synonyms into a list by the "|" delimiter 
		if "||" in synonyms:
			synonyms = "|".join(synonyms.split("||"))

		synonyms = synonyms.split("|")
		
		# Populate Gene Dictionary with info
		gene_dictionary[echinobase_gene_id] = [symbol, name, synonyms, curation_status, info]

# Dictionary of {ECCB-GENEPAGE_ID: [SPU_ID(s)]}. This will be needed when trying to assign Tu et al. 2012 GO Terms to gene_dictionary
def create_ECB_SPU_mapping_dict():
	for key,value in gene_dictionary.items():
		
		# Nine cases where the name (value[1]) is an SPU_ number. 
		if "SPU_" in gene_dictionary[key][1]:
			spu_id = gene_dictionary[key][1]
			ECB_SPU_dict[key] = [spu_id]
		
		# For all other cases, parse the synonyms list for "SPU_" numbers
		else:
			spu_lst = []
			for item in gene_dictionary[key][2]:
				if "SPU_" in item:
					spu_lst.append(item)
				ECB_SPU_dict[key] = spu_lst

def add_GO_KO_termns_to_gene_dictionary():
	# List of GO Terms by ECB-GENEPAGE record 
	gene_go_terms_list = [gene.split("\t") for gene in open(gene_go_terms,"r").read().splitlines()]

	# List of KO Terms by ECB-GENEPAGE record
	gene_ko_terms_list = [gene.split("\t") for gene in open(gene_ko_terms,"r").read().splitlines()]

	# List of Tu et al. 2012 GO Terms by SPU_ number 
	# File has one line header
	tu_go_terms_list = [gene.split(",",1) for gene in open("tu_2012_go.csv","r").read().splitlines()[1:]]

	# Add Echinobase GO Terms to gene dictionary
	# Add GO Terms as a list because there are duplicate ECB-GENEPAGE numbers in the GO Terms file
	for gene in gene_go_terms_list:
		# If the length of the value is 5, GO terms have not been added for this key (ECB-GENEPAGE ID)
		if len(gene_dictionary[gene[0]]) == 5:
			gene_dictionary[gene[0]].append([gene[3]])
		
		# If the length of the value is 6, GO terms have already been added for this key (ECB-GENEPAGE ID)
		# Append the additional GO terms to the pre-existing list 
		elif len(gene_dictionary[gene[0]]) == 6:
			gene_dictionary[gene[0]][5].append(gene[3])

	# Append "n/a" to gene dictionary records if GO Terms weren't appended
	for key,value in gene_dictionary.items():
		if len(value) == 5:
			gene_dictionary[key].append("n/a")

	# Add Echinobase KO Terms to Gene Dictionary 
	# Add KO Terms as a list because there are duplicate ECB-GENEPAGE numbers in the KO Terms file
	for gene in gene_ko_terms_list:
		# If the length of the value is 6, KO terms have not been added for this key (ECB-GENEPAGE ID)
		if len(gene_dictionary[gene[0]]) == 6:
			gene_dictionary[gene[0]].append([gene[3]])
		
		# If the length of the value is 7, KO terms have already been added for this key (ECB-GENEPAGE ID)
		# Append the additional KO terms to the pre-existing list 
		elif len(gene_dictionary[gene[0]]) == 7:
			gene_dictionary[gene[0]][5].append(gene[3])

	# Append "n\a" to gene dictionary records if didn't append KO Terms
	for key,value in gene_dictionary.items():
		if len(value) == 6:
			gene_dictionary[key].append("n/a")

	# Add Tu et al. GO Terms 
	# First, create dictionary of format {SPU_number : GO_Terms}
	tu_go_dict = dict()
	for item in tu_go_terms_list:
		tu_go_dict[item[0]] = item[1]

	# Is this appropriate? When multiple SPU IDs map to one ECB-GENEPAGE ID? Is it appropriate to append all the SPU Tu et al. GO Terms to the one ECB-GENEPAGE ID?
	for gene in gene_dictionary:
		# If there are multiple SPU_ IDs for the given ECB-GENEPAGE number, create go_lst and add the GO terms for all SPU numbers to this list
		if len(ECB_SPU_dict[gene]) > 1:
			spu_lst = ECB_SPU_dict[gene]
			go_lst =[]
			for spu_number in spu_lst:
				if spu_number in tu_go_dict.keys():
					go_lst.append(tu_go_dict[spu_number])
			
			# If there are multiple items (independent strings) in the go_lst, combine these into a single item separated by a "|"
			if len(go_lst) >= 1:
				go = "|".join(go_lst)
				gene_dictionary[gene].append(go)
			
			# If there were no Tu et al. 2012 GO terms, add "n/a" to gene_dictionary
			else:
				gene_dictionary[gene].append("n/a")
			
		
		# If there isn't a corresponding SPU number of gene_dictionary key, add "n/a" for Tu et al. GO Terms
		elif len(ECB_SPU_dict[gene]) == 0:
			gene_dictionary[gene].append("n/a")
		
		# If there is exactly one SPU_ ID for the given ECB-GENEPAGE number 
		elif len(ECB_SPU_dict[gene]) == 1:
			spu = ECB_SPU_dict[gene][0]
			
			# If SPU_ ID has Tu et al GO Terms, add these to gene_dictionary
			if spu in tu_go_dict.keys():
				gene_dictionary[gene].append(tu_go_dict[spu])
			
			# If SPU_ID has no Tu et al. GO Terms, add "n/a" to gene_dictionary
			else:
				gene_dictionary[gene].append("n/a")

	# Because there were duplicated ECB-GENEPAGE records in GeneGoTerms.txt and GeneKoTerms.txt, the same GO term may have been added twice to a given ECB-GENEPAGE record
	# Eliminate redundant items using list(set(GO/KO term list))
	# Change format of GO/KO records from list to string separating independent terms by comma
	for value in gene_dictionary.values():
		
		# When terms were added, if a ECB-GENEPAGE record had no terms, the string "n/a" was added, so not all GO/KO records are a list 
		if type(value[5]) is list:
			value[5] = ",".join(list(set(value[5])))

		if type(value[6]) is list:
			value[6] = ",".join(list(set(value[6])))

	# Deal with redundant Tu et al. terms here!

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file, gene_file, outfile):
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_file + " -wo > " + outfile)

def handle_duplicates_and_reverse_dictionary():
	# Write code to remove records with duplicate symbols or merge them together into a single record 

	# Population LOC_gene_dictionary
	for key,value in gene_dictionary.items():
		if value[0] not in LOC_gene_dictionary.keys():
			LOC_gene_dictionary[value[0]] = [key, value[1], value[2], value[3], value[4], value[5], value[6], value[7]]
		
		# If duplicate entry, search for LOC ID
		else:
			for item in value[2]:
				if "LOC" in item:
					id = item
					break
			
			# Add record to LOC_gene_dictionary, add the previous id/name to the synonyms column
			if id not in LOC_gene_dictionary.keys():
				synonyms = value[2]
				gene_name = value[0]
				synonyms.insert(0, gene_name)
				LOC_gene_dictionary[id] = [key, value[1], synonyms, value[3], value[4], value[5], value[6], value[7]]

			# There are four duplicate records that contain no mapping to a LOC ID.
			#else:
				#print(key,value)

def write_gene_dictionary_to_csv():
	csv_file = open("all_genes.csv","w")
	writer = csv.writer(csv_file)	
	header = ["symbol", "echinobase_gene_id", "name", "synonyms", "curation status", "info", "ECB GO Terms", "ECB KO Terms", "Tu et al. GO Terms"]
	writer.writerow(header)

	# Write gene info to csv file
	for key,value in LOC_gene_dictionary.items():
		data = [key, value[0], value[1], "|".join(value[2]), value[3], value[4], value[5], value[6], value[7]]
		writer.writerow(data)

def make_coverage_depths_dict():
	coverage_depths = open(tract_coverage_file,"r").read().splitlines()
	for record in coverage_depths:
		tract_name = record.split("\t")[0].replace("_","-")
		dro = record.split("\t")[1]
		fra = record.split("\t")[2]
		pal = record.split("\t")[3]
		pul = record.split("\t")[4]
		tract_coverage_depth_dict[tract_name] = [dro,fra,pal,pul]

def create_gene_intersection_dict(overlap_file):
	overlaps = open(overlap_file,"r").read().splitlines()

	for record in overlaps:
		record = record.split("\t")
		tract_id = record[3].replace("_","-")
		#LOC ID
		gene_id = record[7].split("-",1)[1]

		dro_cov = tract_coverage_depth_dict[tract_id][0]
		fra_cov = tract_coverage_depth_dict[tract_id][1]
		pal_cov = tract_coverage_depth_dict[tract_id][2]
		pul_cov = tract_coverage_depth_dict[tract_id][3]

		gene_name = "missing"
		
		# Edge case for REJ8 and Coel1, which are rej8 and coel1 in the dictionary
		if gene_id.lower() in LOC_gene_dictionary.keys():
			gene_id = gene_id.lower()

		if gene_id in LOC_gene_dictionary.keys():
			gene_ecb_id = LOC_gene_dictionary[gene_id][0]
			gene_name = LOC_gene_dictionary[gene_id][1]
			gene_synonyms = LOC_gene_dictionary[gene_id][2]
			gene_curation_status = LOC_gene_dictionary[gene_id][3]
			gene_info = LOC_gene_dictionary[gene_id][4]
			gene_ECB_GO = LOC_gene_dictionary[gene_id][5]
			gene_ECB_KO = LOC_gene_dictionary[gene_id][6]
			gene_tu_GO = LOC_gene_dictionary[gene_id][7]
		
		# If gene_id isn't in LOC_gene_dictionary.keys(), check to see if it is the synonyms of another record
		else:
			# I verified that this doesn't lead to finding multiple records and overwriting. 
			for key,value in LOC_gene_dictionary.items():
				if gene_id in value[2]:
					gene_ecb_id = value[0]
					gene_name = value[1]
					gene_synonyms = value[2]
					gene_curation_status = value[3]
					gene_info = value[4]
					gene_ECB_GO = value[5]
					gene_ECB_KO = value[6]
					gene_tu_GO = value[7]
		
		gene_coordinates = str(record[4]) + ":" + str(record[5]) + "-" + str(record[6])
		gene_length = int(record[6]) - int(record[5])
		bp_overlap = record[14]
		percent_introgressed = (int(bp_overlap) / gene_length)*100

		# Populate gene_intersection_dict with all of the assigned variables
		gene_intersection_dict[gene_id] = [tract_id, dro_cov, fra_cov, pal_cov, pul_cov, gene_name, gene_length, bp_overlap, percent_introgressed, gene_coordinates, gene_ecb_id, gene_synonyms, gene_curation_status, gene_info, gene_ECB_GO, gene_ECB_KO, gene_tu_GO]

def create_gene_intersection_file():
	csv_file = open("intersect.csv","w")
	writer = csv.writer(csv_file)	
	header = ["introgression_tract", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Gene", "Gene Name", "Gene Length", "Overlapping Bases", "Percent Introgressed", "Gene Coordinates", "ECB Gene ID", "Gene Synonyms", "Curation Status", "Gene Info", "Echinobase GO", "Echinobase KO", "Tu GO"]
	writer.writerow(header)
	
	for key,value in gene_intersection_dict.items():
		data = [value[0], value[1], value[2], value[3], value[4], key, value[5], value[6], value[7], value[8], value[9], value[10], "|".join(value[11]), value[12], value[13], value[14], value[15], value[16]]
		writer.writerow(data)
	
	csv_file.close()

def write_introgressed_genes_to_bed():
	with open("introgressed_genes.bed","w") as f:
		for key,value in gene_intersection_dict.items():
			if value[8] == 100:
				scaffold = value[9].split(":")[0]
				start = value[9].split(":")[1].split("-")[0]
				stop = value[9].split(":")[1].split("-")[1]
				f.write(scaffold + "\t" + start + "\t" + stop + "\t" + key + "\t" + value[5] + "\n")

def get_exon_overlap():
	overlapping_bases = []

	with open("exon_overlap.bed","r") as f:
		inputs = f.read().splitlines()

	for overlap in inputs:
		overlapping_bases.append(int(overlap.split("\t")[14]))

	return sum(overlapping_bases)

def get_gene_overlap():
	overlapping_bases = []

	with open("gene_overlap.bed","r") as f:
		inputs = f.read().splitlines()

	for overlap in inputs:
		overlapping_bases.append(int(overlap.split("\t")[14]))

	return sum(overlapping_bases)


def main():
	create_gene_dictionary()
	create_ECB_SPU_mapping_dict()
	add_GO_KO_termns_to_gene_dictionary()
	
	intersect_genes(tract_file, gene_list_file, "gene_overlap.bed")
	intersect_genes(tract_file, exon_list_file, "exon_overlap.bed")
	
	handle_duplicates_and_reverse_dictionary()

	write_gene_dictionary_to_csv()

	make_coverage_depths_dict()

	create_gene_intersection_dict("gene_overlap.bed")
	#create_gene_intersection_dict("exon_overlap.bed")

	create_gene_intersection_file()
	#create_exon_intersection_file()

	overlapping_coding_bases = get_exon_overlap()
	overlapping_genic_bases = get_gene_overlap()
	overlapping_intronic_bases = overlapping_genic_bases - overlapping_coding_bases
	overlapping_intergenic_bases = total_bases_introgressed - overlapping_genic_bases

	print("Overlapping genic bases: {}".format(overlapping_genic_bases))
	print("Overlapping coding bases: {}".format(overlapping_coding_bases))
	print("Overlapping intronic bases: {}".format(overlapping_intronic_bases))
	print("Overlapping intergenic bases: {}".format(overlapping_intergenic_bases))

	print("Total introgressed bases: {}".format(total_bases_introgressed))
	print(overlapping_coding_bases + overlapping_intronic_bases + overlapping_intergenic_bases)

	write_introgressed_genes_to_bed()

if __name__ == "__main__":
	main()
