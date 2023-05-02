import os
import csv 
import gzip
from joblib import Parallel, delayed

process_hmm_output_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/process_hmm_species_tree/"
working_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts_species_tree/"
phylonet_hmm_alignment_dir = "/hb/scratch/mglasena/phylonet_hmm/hmm_input/scaffold_nexus_alignments/"

# Reference alignment BAM files for assessing coverage dpeth
bam_file_paths_list = [
"/hb/groups/pogson_group/dissertation/data/bam_files/droebachiensis_SRR5767286_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/fragilis_SRR5767279_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pallidus_SRR5767285_dedup_aligned_reads.bam", 
"/hb/groups/pogson_group/dissertation/data/bam_files/pulcherrimus_SRR5767283_dedup_aligned_reads.bam"
]

threads = 4

introgressed_gene_coverage_dict = dict()

# Total bases introgressed
tract_length_dist = list(csv.reader(open(process_hmm_output_dir + "tract_length_dist.csv","r")))
total_bases_introgressed = sum([int(item) for item in tract_length_dist[0]])
print("Total Bases Introgressed: {}".format(total_bases_introgressed))

# Introgression tract file in bed format
tract_file = process_hmm_output_dir + "tracts_pf.bed"

# Coverage depth of each tract for each species 
tract_coverage_file = process_hmm_output_dir + "tract_coverage.tsv"

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

# CSV file of 6,520 single copy orthologs from Kober and Pogson. Rows 5 - 1,012 contain the 1,008 genes with significant tests for positive selection
psg_file = "psg.csv"

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

coordinate_by_scaffold_dict = dict()

tract_intersection_dict = dict()

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
	print("Intersecting {} with {}. Writing results to {}".format(tract_file, gene_file, outfile))
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
	print("Adding {} records to tract_coverage_depth_dict".format(len(coverage_depths)))
	for record in coverage_depths:
		tract_name = record.split("\t")[0].replace("_","-")
		dro = record.split("\t")[1]
		fra = record.split("\t")[2]
		pal = record.split("\t")[3]
		pul = record.split("\t")[4]
		tract_coverage_depth_dict[tract_name] = [dro,fra,pal,pul]

def get_coordinate_file_paths():
	print("Indexing scaffold coordinate files from {}".format(phylonet_hmm_alignment_dir))
	find_files = "find {} -type f -name '*coordinates*' > coordinate_files".format(phylonet_hmm_alignment_dir)
	os.system(find_files)
	coordinate_file_list = open("coordinate_files","r").read().splitlines()
	os.system("rm coordinate_files")
	return coordinate_file_list

def create_coordinate_dict(coordinate_file_list):
	for coordinate_file in coordinate_file_list:
		scaffold = open(coordinate_file,"r").readline().split(":")[0]
		coordinate_by_scaffold_dict[scaffold] = [item.split(":")[1] for item in open(coordinate_file,"r").read().splitlines()]

def create_gene_intersection_dict(overlap_file):
	overlaps = open(overlap_file,"r").read().splitlines()
	print("Number of records in overlap file: {}".format(len(overlaps)))

	for record in overlaps:
		record = record.split("\t")
		#LOC ID
		gene_id = record[7].split("-",1)[1]

		if gene_id in gene_intersection_dict:
			additional_tract = record[3].replace("_","-")
			bp_overlap = int(record[14])
			gene_intersection_dict[gene_id][5].append(additional_tract)
			gene_intersection_dict[gene_id][2] += bp_overlap
			gene_intersection_dict[gene_id][3] = gene_intersection_dict[gene_id][2] / gene_intersection_dict[gene_id][1] * 100

		else:
			tract_id = [record[3].replace("_","-")]

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
			bp_overlap = int(record[14])
			percent_introgressed = bp_overlap / gene_length * 100

			# Populate gene_intersection_dict with all of the assigned variables
			gene_intersection_dict[gene_id] = [gene_name, gene_length, bp_overlap, percent_introgressed, gene_coordinates, tract_id, gene_ecb_id, gene_synonyms, gene_curation_status, gene_info, gene_ECB_GO, gene_ECB_KO, gene_tu_GO]

	print("Length of gene_intersection_dict: {}".format(len(gene_intersection_dict)))

def create_tract_info_file(overlap_file):
	overlaps = open(overlap_file,"r").read().splitlines()
	value_error_counter = 0 
	
	for key,value in tract_coverage_depth_dict.items():
		scaffold = key.split(":")[0].replace("-","_")
		start = str(int(key.split(":")[1].split("-")[0]) + 1)
		stop = str(key.split(":")[1].split("-")[1])
		
		if start in coordinate_by_scaffold_dict[scaffold] and stop in coordinate_by_scaffold_dict[scaffold]:
			SNV_sites = int(coordinate_by_scaffold_dict[scaffold].index(stop)) - int(coordinate_by_scaffold_dict[scaffold].index(start)) + 1
		
		else:
			print(key)
			print(start)
			print(stop)
			
			if start in coordinate_by_scaffold_dict[scaffold]:
				start_pos = int(start)
				print("Start in coordinate_by_scaffold_dict: {}".format(start_pos))
			else:
				print("Finding new Start")
				start_pos = 0
				i = 0
				while start_pos < int(start):
					start_pos = int(coordinate_by_scaffold_dict[scaffold][i])
					i += 1

				start_pos = int(coordinate_by_scaffold_dict[scaffold][i+1])
			
				print("New Start: {}".format(start_pos))

			if stop in coordinate_by_scaffold_dict[scaffold]:
				print("Stop in coordinate_by_scaffold_dict: {}".format(stop_pos))
				stop_pos = int(stop)
			else:
				print("Finding New Stop")
				stop_pos = int(coordinate_by_scaffold_dict[scaffold][-1])
				i = len(coordinate_by_scaffold_dict[scaffold]) - 1
				while stop_pos > int(stop):
					i -= 1
					stop_pos = int(coordinate_by_scaffold_dict[scaffold][i])
			
				stop_pos = int(coordinate_by_scaffold_dict[scaffold][i-1])

				print("New Stop: {}".format(stop_pos))

				SNV_sites = int(coordinate_by_scaffold_dict[scaffold].index(str(stop_pos))) - int(coordinate_by_scaffold_dict[scaffold].index(str(start_pos))) + 1

		tract_coverage_depth_dict[key].append(SNV_sites)
		
	for overlap in overlaps:
		tract_id = overlap.split("\t")[3].replace("_","-")
		gene_id = overlap.split("\t")[7].split("gene-")[1]
		if tract_id in tract_intersection_dict:
			tract_intersection_dict[tract_id].append(gene_id)
		else:
			tract_intersection_dict[tract_id] = [gene_id]

	for key,value in tract_coverage_depth_dict.items():
		if key in tract_intersection_dict:
			tract_coverage_depth_dict[key].append(tract_intersection_dict[key])
		else:
			tract_coverage_depth_dict[key].append(["n/a"])

	csv_file = open("tract_info.csv","w")
	writer = csv.writer(csv_file)	

	header = ["tract_name", "Scaffold", "Start", "Stop", "Length", "SNV Sites", "SNV/bp", "Sdro", "Sfra", "Spal", "Hpul",  "Overlapping Genes"]
	writer.writerow(header)
	
	for key,value in tract_coverage_depth_dict.items():
		scaffold = key.split(":")[0]
		start = key.split(":")[1].split("-")[0]
		stop = key.split(":")[1].split("-")[1]
		length = int(stop) - int(start)
		snv_per_bp = int(value[4])/length
		data = [key, scaffold, start, stop, length, value[4], snv_per_bp, value[0], value[1], value[2], value[3], ",".join(value[5])]
		writer.writerow(data)
	
	csv_file.close()

def write_introgressed_genes_to_bed():
	with open("introgressed_genes.bed","w") as f:
		for key,value in gene_intersection_dict.items():
				scaffold = value[4].split(":")[0]
				start = value[4].split(":")[1].split("-")[0]
				stop = value[4].split(":")[1].split("-")[1]
				f.write(scaffold + "\t" + str(start) + "\t" + str(stop) + "\t" + key + "\n")

def run_mosdepth(tract_file, bam_file):
	prefix = bam_file.split("/")[-1].split("_dedup")[0]
	mosdepth = "mosdepth --by " + tract_file + " --no-per-base --thresholds 1,10,20 -t {} --fast-mode {} {}".format(threads, prefix, bam_file)
	os.system(mosdepth)

# Get list of mosdepth output file paths in alphabetic order 
def get_mosdepth_output_file_list():
	find_files = "find {} -type f -name '*.regions*' | grep -v 'csi' > mosdepth_output_files".format(working_dir)
	os.system(find_files)
	mosdepth_output_file_list = open("mosdepth_output_files","r").read().splitlines()
	os.system("rm mosdepth_output_files")
	return sorted(mosdepth_output_file_list)

# Populate coverage_dict in the format of {tract_name: coverage_species1, coverage_species2, coverage_species3, coverage_species4}
def parse_mosdepth(mosdepth_region_file):
	with gzip.open(mosdepth_region_file,"rt") as f:
		inputs = f.read().splitlines()
	
	for record in inputs:
		gene = record.split("\t")[3]
		coverage = float(record.split("\t")[4])

		if gene not in introgressed_gene_coverage_dict:
			introgressed_gene_coverage_dict[gene] = [coverage]
		else:
			introgressed_gene_coverage_dict[gene].append(coverage)

def update_gene_intersection_dict():
	for key,value in gene_intersection_dict.items():
		values = value
		gene_intersection_dict[key] = [values[0], values[1], values[2], values[3], values[4], introgressed_gene_coverage_dict[key][0], introgressed_gene_coverage_dict[key][1], introgressed_gene_coverage_dict[key][2], introgressed_gene_coverage_dict[key][3], values[5], values[6], values[7], values[8], values[9], values[10], values[11], values[12]]

	with open("fully_introgressed_genes_coverage.tsv","w") as f:
		for key,value in introgressed_gene_coverage_dict.items():
			if gene_intersection_dict[key][3] == 100:
				f.write(key + "\t" + str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" + str(value[3]) + "\n")

def create_gene_intersection_file():
	csv_file = open("intersect.tsv","w")
	writer = csv.writer(csv_file, delimiter="\t")	
	header = ["NCBI Gene ID", "Name", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Introgression Tract(s)", "ECB Gene ID", "Gene Synonyms", "Curation Status", "Info", "Echinobase GO", "Echinobase KO", "Tu GO"]
	writer.writerow(header)
	
	for key,value in gene_intersection_dict.items():
		data = [key, value[0], value[1], value[2], value[3], value[4], value[5], value[6], value[7], value[8], ",".join(value[9]), value[10], "|".join(value[11]), value[12], value[13], value[14], value[15], value[16]]
		writer.writerow(data)
	
	csv_file.close()

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

def find_psg_overlap():
	psg_list = list(csv.reader(open(psg_file,"r"), delimiter = ","))
	introgressed_genes_list = list(csv.reader(open("intersect.tsv","r"), delimiter = "\t"))

	# Get list of SPU identifiers from the set of positively selected genes 
	psg_list = [n for n in psg_list[4:1012]]

	csv_file = open("psg_intersect.tsv","w")
	writer = csv.writer(csv_file, delimiter = "\t")	
	header = ["NCBI Gene ID", "Name", "Synonyms", "Kober and Pogson Gene ID", "Kober and Pogson Name", "Kober and Pogson Synonyms", "PSG #", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Overlapping introgression tract(s)", "Sdro", "Sfra", "Spal", "Hpul"]
	writer.writerow(header)
	
	for gene in psg_list:
			for record in introgressed_genes_list:
				if gene[0] in record or gene[1] in record or gene[2] in record:
					data = [record[0], record[1], record[12], gene[0], gene[1], gene[2], psg_list.index(gene), record[2], record[3], record[4], record[5], record[10], record[6], record[7], record[8], record[9]]
					writer.writerow(data)
	
	csv_file.close()

def main():
	os.chdir(working_dir)
	create_gene_dictionary()
	create_ECB_SPU_mapping_dict()
	add_GO_KO_termns_to_gene_dictionary()
	
	intersect_genes(tract_file, gene_list_file, "gene_overlap.bed")
	intersect_genes(tract_file, exon_list_file, "exon_overlap.bed")
	
	handle_duplicates_and_reverse_dictionary()

	write_gene_dictionary_to_csv()

	make_coverage_depths_dict()

	coordinate_files = get_coordinate_file_paths()
	create_coordinate_dict(coordinate_files)

	create_gene_intersection_dict("gene_overlap.bed")

	#create_gene_intersection_dict("exon_overlap.bed")

	create_tract_info_file("gene_overlap.bed")

	write_introgressed_genes_to_bed()

	Parallel(n_jobs=len(bam_file_paths_list))(delayed(run_mosdepth)("introgressed_genes.bed", bam_file) for bam_file in bam_file_paths_list)

	mosdepth_output_file_list = get_mosdepth_output_file_list()

	for file in mosdepth_output_file_list:
		parse_mosdepth(file)

	update_gene_intersection_dict()

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

	find_psg_overlap()

if __name__ == "__main__":
	main()
