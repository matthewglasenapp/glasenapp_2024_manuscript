import os
import gzip
from itertools import islice
import csv
from joblib import Parallel, delayed
import multiprocessing
import subprocess
import statistics
import multiprocessing
from Bio import SeqIO

num_cores = 24

# nw_utils directory
nw_utils = "/hb/home/mglasena/software/newick_utils/src/"

# Path to vcf2fasta.py
vcf2fasta = "/hb/home/mglasena/software/vcf2fasta/vcf2fasta.py"

# Feature of gff file that vcf2fasta.py will build alignments for
feature = "CDS"

# Path to S. purpuratus reference genome
reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# S. purpuratus gff3 file
gff_file = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.gff"

# Path to filtered multisample vcf file
vcf_file = "/hb/scratch/mglasena/invariant_sites_vcf_phylonet_hmm/combined_vcf_files/filtered_genotype_calls_individual_genotypes.g.vcf.gz"

# Bed file containing a record for each protein coding gene in the S. purpuratus assembly. See the ncbi/ directory for scripts to generate this file
protein_coding_genes_bed_file = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/protein_coding_genes.bed"

# Directory containing the output files from mrna.py
bed_file_dir = "/hb/home/mglasena/dissertation/data/mosdepth/mrna_cov/"

# Sample bed file from the output of mrna.py
bed_file = "/hb/home/mglasena/dissertation/data/mosdepth/mrna_cov/pallidus_SRR5767285.regions.bed.gz"

# File containing genes with the majority of ther bases introgressed 
introgressed_genes = "majority_bases_introgressed_genes_coverage.tsv"

# Specify species to include for ortholog finder. MUST BE ALPHABETICAL!
# Strongylocentrotidae Subset
#subset_sample_list = ['depressus_SRR5767284', 'droebachiensis_SRR5767286', 'fragilis_SRR5767279', 'franciscanus_SRR5767282', 'intermedius_SRR5767280', 'nudus_SRR5767281', 'pallidus_SRR5767285', 'pulcherrimus_SRR5767283', 'purpuratus_SRR7211988']
subset_sample_list = ['fragilis_SRR5767279', 'pulcherrimus_SRR5767283']

# Specify thresholds for filtering. 
min_cov_threshold = 5

prop_1x_threshold = 0.75

prop_10x_threshold = 0.0

# Required gap between two adjacent loci in base pairs 
required_gap = 0

# Average coverage of S. purpuratus exons for each sample
mean_coverage_spur5_exons = {
"depressus_SRR5767284": 47.5,
"droebachiensis_SRR5767286" : 41.5,
"fragilis_SRR5767279": 46.8,
"franciscanus_SRR5767282": 33.8,
"intermedius_SRR5767280": 44.2,
"nudus_SRR5767281": 40.5,
"pallidus_SRR5767285": 15,
"pulcherrimus_SRR5767283": 44.3,
"purpuratus_SRR7211988": 100.3,
}

# Mapping of DNA sample names to species names 
sample_names = {
#'(4': "(lividus",
'QB3KMK013': 'fragilis',
'QB3KMK011': 'nudus',
'QB3KMK010': 'franciscanus',
'QB3KMK015': 'depressus',
'QB3KMK002': 'pallidus',
'QB3KMK014': 'droebachiensis',
'QB3KMK016': 'pulcherrimus',
'QB3KMK012': 'intermedius',
'SPUR.00': 'purpuratus',
}

# Initialize new dictionary containing the average coverage of S. purpuratus exons for each sample in subset_sample_list variable
subset_mean_coverage_spur5_exons = dict()

# Initialize dictionary in format of {"rna": [average coverage depth for each species], [prop_1x], [prop_10x], [prop_20x])}
rna_dict = dict()

# Initialize dictionary for mRNAs passing initial filter by coverage depth 
passed_rna_dict = dict()

# Initialize dictionary mapping mRNAs passing intial filter by coverage depth to the name of their parent gene
mrna_gene_dict = dict()

# Initialize dictionay in format of "Scaffold":[[scaffold,rna,start,stop]] for mRNAs passing initial filter. 
# This dictionary will be used to filter out mRNAs that are too close to one another
scaffold_dict = {}

# Initialize final list for genes passing both the coverage depth and gap filters
passed_rnas = []

# Initialize dictionary mapping mRNA names to parent gene name for mRNAs that passed all filters
filtered_mrna_gene_dict = dict()

# Initialize dictionary mapping CDS names to parent mRNA names for CDS of genes that had an mRNA pass previous filters. 
cds_parent_rna_dict = dict()

# Populate subset_mean_coverage_spur5_exons variable with the samples included in the subset_sample_list variable
def subset_coverage_dict():
	for sample in subset_sample_list:
		subset_mean_coverage_spur5_exons[sample] = mean_coverage_spur5_exons[sample]

# Get zipped list of regions and thresholds files for each species 
def get_zipped_bed_file_list():
	get_regions_file_paths = "find {} -type f -name *.regions.bed.gz* > regions_files".format(bed_file_dir)
	os.system(get_regions_file_paths)

	get_thresholds_file_paths = "find {} -type f -name *.thresholds.bed.gz* > thresholds_files".format(bed_file_dir)
	os.system(get_thresholds_file_paths)

	with open("regions_files", "r") as f1, open("thresholds_files","r") as f2:
		file_list = zip(sorted(f1.read().splitlines()),sorted(f2.read().splitlines()))

	os.system("rm regions_files")
	os.system("rm thresholds_files")

	return list(file_list)

# Populate keys and create structure for rna_dict in format of {"rna": [average coverage depth for each species], [prop_1x], [prop_10x], [prop_20x])}
def initialize_rna_dict():
	with gzip.open(bed_file,"rt") as f:
		for line in f:
			gene = line.split("\t")[0]
			rna_dict[gene] = [],[],[],[]

# Populate rna dict with coverage metrics for each species. No filtering involved in this step 
def fill_rna_dict(regions_file, thresholds_file):
	with gzip.open(regions_file, "rt") as file1, gzip.open(thresholds_file, "rt") as file2:
		for line_file1, line_file2 in zip(file1, islice(file2, 0, None)):
			gene_id = line_file1.split("\t")[-2]
			mean_depth = float(line_file1.split("\t")[-1].strip())

			total_base_count = int(line_file2.split("\t")[-4])
			
			one_x_count = int(line_file2.split("\t")[-3])
			ten_x_count = int(line_file2.split("\t")[-2])
			twenty_x_count = int(line_file2.split("\t")[-1])
			
			try:
				prop_1x = float(one_x_count / total_base_count)
				prop_10x = float(ten_x_count / total_base_count)
				prop_20x = float(twenty_x_count / total_base_count)
			except ZeroDivisionError:
				prop_1x = 0.0
				prop_10x = 0.0
				prop_20x = 0.0

			rna_dict[gene_id][0].append(mean_depth)
			rna_dict[gene_id][1].append(prop_1x)
			rna_dict[gene_id][2].append(prop_10x)
			rna_dict[gene_id][3].append(prop_20x)

# Write rna_dict to csv file 
def write_all_rna_dict_csv():
	csv_file = open("all_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	
	# Create header row for csv file
	for i in range(4):
		for sample in subset_sample_list:
			header.append(sample)
	
	writer.writerow(header)

	# Write row for each mRNA containing mean depth, prop1x, prop10x, and prop20x for each species/sample
	counter = 0
	for key,value in rna_dict.items():
		row = [key] + value[0] + value[1] + value[2] + value[3]
		writer.writerow(row)
		counter += 1

	csv_file.close()
	
	print("There were {} mRNA records pre-filter".format(len(rna_dict)))
	print("{} mRNA records written to all_rna.csv".format(counter))

# Filter rna_dict by coverage metrics. Add mRNA's passing filter to passed_rna_dict
def filter_rna_dict():
	record_counter = 0
	failed_counter = 0
	passed_counter = 0

	for key, value in rna_dict.items():
		record_counter += 1
		mean_depth_lst = [item for item in value[0]]
		one_x_lst = [item for item in value[1]]
		ten_x_lst = [item for item in value[2]]
		twenty_x_lst = [item for item in value[3]]
		
		if min(mean_depth_lst) >= min_cov_threshold and min(one_x_lst) >= prop_1x_threshold and min(ten_x_lst) >= prop_10x_threshold:
			test_var = True
			
			# Filter any species who's coverage for that mRNA is double or greater the average coverage for exons
			for counter, sample in enumerate(subset_mean_coverage_spur5_exons):
				if mean_depth_lst[counter] >= (subset_mean_coverage_spur5_exons[sample] * 2):
					test_var = False
				
			if test_var == True:
				passed_rna_dict[key] = value
				passed_counter += 1 
			else:
				failed_counter += 1
		else:
			failed_counter += 1

	print("{} mRNA records processed".format(record_counter))
	print("{} mRNA records passed intial coverage depth filters".format(passed_counter))
	print("{} mRNA records failed intial coverage depth filters".format(failed_counter))
	print("{} mRNA records in passed_rna_dict".format(len(passed_rna_dict)))

# Get file of mRNA records from gff3 file
def get_mrna_gff():
	get_mrna_from_gff = '''awk '$3 == "mRNA"' {} > mrna_records.txt'''.format(gff_file)
	os.system(get_mrna_from_gff)
	
# Create dictionary for mRNAs passing initial filter in the format of "Scaffold":[[scaffold,rna,start,stop]]. This dictionary will be used for filtering by distance on chromosome
def create_scaffold_dict():
	
	# Creates mapping of rna: parent_gene for mRNAs in passed_rna_dict 
	mt_rna_counter = 0
	record_counter = 0

	with open("mrna_records.txt", "r") as f:
		for line in f:
			rna = line.split("\t")[8].split(";")[0].split("rna-")[1]
			if rna in passed_rna_dict.keys():
				scaffold = line.split("\t")[0]
				start = line.split("\t")[3]
				stop = line.split("\t")[4]
				parent_gene = line.split("\t")[8].split(";")[1].split("gene-")[1]
				
				if scaffold == "NC_001453.1":
					mt_rna_counter += 1

				else:
					mrna_gene_dict[rna] = parent_gene

					if scaffold_dict.get(scaffold):
						scaffold_dict[scaffold].append([scaffold,rna,start,stop])
					else:
						scaffold_dict[scaffold] = [[scaffold,rna,start,stop]]

					record_counter += 1

	os.system("rm mrna_records.txt")
	print("{} mitochondrial mRNAs were removed".format(mt_rna_counter))
	print("{} mRNA records in scaffold_dict".format(record_counter))
	print("{} mRNA records in mrna_gene_dict".format(len(mrna_gene_dict)))

# Get list of passed mRNAs and create new filtered_mrna_gene_dict that only includes mRNAs that passed the gap filter. 
def get_passed_rnas():
	for rna_list in scaffold_dict.values():
		for rna in rna_list:
			passed_rnas.append(rna[1])

	for key,value in mrna_gene_dict.items():
		if key in passed_rnas:
			filtered_mrna_gene_dict[key] = value

	print("{} mRNA records in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))

# Write new bed file ("unlinked_loci.bed") of parent genes with an mRNA transcript that passed all filters. This will be used to build the initial vcf2fasta alignments 
def write_new_bed_file():
	records_written = 0
	records_written_lst = []
	
	with open(protein_coding_genes_bed_file,"r") as f:
		with open("unlinked_loci.bed","a") as f2:
			for line in f:
				if line.split("\t")[3].split("gene-")[1] in filtered_mrna_gene_dict.values():
					records_written_lst.append(line.split("\t")[3].split("gene-")[1])
					records_written += 1
					f2.write(line)

	print("{} records written to unlinked_loci.bed".format(records_written))

def cross_reference():
	gene_lst = open(introgressed_genes,"r").read().splitlines()
	gene_lst = [item.split("\t")[0] for item in gene_lst]
	print("{} genes with the majority of their bases introgressed".format(len(gene_lst)))

	passed_rnas = open("unlinked_loci.bed","r").read().splitlines()
	passed_rnas = [item.split("\t")[3].split("gene-")[1] for item in passed_rnas]

	# Initialize list for genes that have a single-copy mRNA that passes filter
	sco_introgressed_genes = []
	
	for gene in gene_lst:
		if gene in passed_rnas:
			sco_introgressed_genes.append(gene)

	# Write new bed file for passing gene records
	original_bed_file = open("unlinked_loci.bed", "r").read().splitlines()
	with open("passed_records_majority_introgressed","w") as f:
		for line in original_bed_file:
			if line.split("\t")[3].split("gene-")[1] in sco_introgressed_genes:
				f.write(line + "\n")

	print("{} genes with the majority of their bases introgressed had single-copy mRNAs that passed filter".format(len(sco_introgressed_genes)))

# Get list of parent gene identifiers for those genes that passed all filters. Example: Dbxref=GeneID:582406
def get_gene_ids():
	get_info_column = "awk '{ print $10 }' passed_records_majority_introgressed > gene_list"
	os.system(get_info_column)

	with open("gene_list","r") as f, open("gene_ids","a") as f2:
		gene_list = f.read().splitlines()
		for gene in gene_list:
			identifier = gene.split(";")[1]
			f2.write(identifier + "\n")

	os.system("rm gene_list")
	gene_ids = open("gene_ids", "r").read().splitlines()
	os.system("rm gene_ids")
	return gene_ids

# Using S. purpuratus gff3 file, make gff file for each gene that passed previous filters
def make_sco_gff(gene):
	command = "grep {} {} > single_gene_gff_records/{}.record".format(gene, gff_file, gene)
	os.system(command)

# Use vcf2fasta to create fasta alignments for all mRNAs passing filter
def run_vcf2fasta():
	run_vcf2fasta = "{} --fasta {} --vcf {} --gff sco_gff.gff --feat {} --blend".format(vcf2fasta, reference_genome, vcf_file, feature)
	os.system(run_vcf2fasta)

# Vcf2fasta makes alignments for all isoforms of each gene. This function identifies the unwanted isoform alignment files and deletes them.
def remove_redundant_isoforms():
	# List of redundant isoforms to delete following vcf2fasta
	records_to_delete = []

	get_cds_gff = '''awk '$3 == "CDS"' sco_gff.gff > cds.gff'''
	os.system(get_cds_gff)

	records = open("cds.gff","r").read().splitlines()

	for record in records:
		cds_name = record.split("\t")[8].split(";")[0].split("cds-")[1]
		parent_rna_name = record.split("\t")[8].split(";")[1].split("rna-")[1]
		cds_parent_rna_dict[cds_name] = parent_rna_name
	
	passed_rnas_lst = list(filtered_mrna_gene_dict.keys())

	for key,value in cds_parent_rna_dict.items():
		if not value in passed_rnas_lst:
			records_to_delete.append(key)

	for record in records_to_delete:
		cds_parent_rna_dict.pop(record)
		delete = "rm vcf2fasta_CDS/{}.fas".format(record)
		os.system(delete)

	#os.system("rm sco_gff.gff")
	#os.system("rm cds.gff")

def remove_not_multiple_of_three():
	get_file_lst = 'find ./vcf2fasta_CDS/ -type f -name "*.fas" > fasta_file_list'
	os.system(get_file_lst)

	cds_length_dict = dict()
	not_multiple_of_three_counter = 0
	not_multiple_of_three_lst = []
	passed_cds_length_dict = dict()

	files = open("fasta_file_list", "r").read().splitlines()

	os.system("rm fasta_file_list")

	for file in files:
		cds = file.split("/")[-1].split(".fas")[0]
		length = len(open(file,"r").read().splitlines()[1])
		cds_length_dict[cds] = length

	for key,value in cds_length_dict.items():
		if value % 3 != 0:
			not_multiple_of_three_counter += 1
			not_multiple_of_three_lst.append(key)
		else:
			passed_cds_length_dict[key] = value

	os.system("mkdir not_multiple_of_three")
	
	for cds in not_multiple_of_three_lst:
		rna = cds_parent_rna_dict[cds]
		passed_rnas.remove(rna)
		filtered_mrna_gene_dict.pop(rna)
		cds_parent_rna_dict.pop(cds)
		move = "mv vcf2fasta_CDS/{}.fas not_multiple_of_three/".format(cds)
		os.system(move)

	print("Number of concatenated CDS records: {}".format(len(cds_length_dict)))
	print("Number of concatenated CDS records that are not a multiple of 3: {}".format(not_multiple_of_three_counter))

def write_passed_rna_dict_csv():
	csv_file = open("passed_rna.csv","w")
	writer = csv.writer(csv_file)	
	header = ["mRNA"]
	records_written = 0 
	
	for i in range(4):
		for sample in subset_sample_list:
			header.append(sample)
	
	writer.writerow(header)

	for key,value in passed_rna_dict.items():
		if key in filtered_mrna_gene_dict.keys():
			row = [key] + value[0] + value[1] + value[2] + value[3]
			records_written += 1
			writer.writerow(row)

	csv_file.close()

	print("{} records written to passed_rna.csv".format(records_written))

	print("{} mRNAs in passed_rnas list".format(len(passed_rnas)))
	print("{} mRNAs in filtered_mrna_gene_dict".format(len(filtered_mrna_gene_dict)))

	with open("passed_rnas.txt","a") as f:
		for record in passed_rnas:
			f.write(record + "\n")

# Iqtree does not tolerate the '*'' symbol. Replace '*' with 'N'
def replace_missing_genotype_char():
	replace_missing_genotypes = r'find ./vcf2fasta_CDS/ -type f -exec sed -i.bak "s/\*/N/g" {} \;'
	os.system(replace_missing_genotypes)
	
	delete_bak_files = 'find ./vcf2fasta_CDS/ -type f -name "*.bak" -delete'
	os.system(delete_bak_files)

	replace_question_character = r'find ./vcf2fasta_CDS/ -type f -exec sed -i.bak "s/\?/N/g" {} \;'
	os.system(replace_question_character)
	
	os.system(delete_bak_files)

def get_fasta_alignment_paths():
	fasta_files = os.listdir("vcf2fasta_CDS/")
	fasta_file_paths_lst = ["/hb/scratch/mglasena/test_rna_metrics/vcf2fasta_CDS/" + item for item in fasta_files]
	return fasta_file_paths_lst

# def convert_fasta_to_phylip(input_fasta_file):
# 	output_parent_dir = "single_copy_ortholog_fasta_alignments/"
# 	make_output_dir = "mkdir -p {}".format(output_parent_dir)
# 	os.system(make_output_dir)

# 	output_dir = "/hb/scratch/mglasena/test_rna_metrics/" + output_parent_dir + input_fasta_file.split("/")[-1].split(".fas")[0] + "/"
# 	os.mkdir(output_dir)

# 	record_dict = {key: value.translate(table='Standard', stop_symbol='*') for key,value in SeqIO.to_dict(SeqIO.parse(input_fasta_file, "fasta")).items()}

# 	for key,value in record_dict.items():
# 		value.id = key

# 	records = list(record_dict.values())

# 	SeqIO.write(records, output_dir + "consensAlignAA.ordered.phylip", "phylip-sequential")

def convert_fasta_to_phylip(input_fasta_file):
	output_parent_dir = "single_copy_ortholog_fasta_alignments/"
	make_output_parent_dir = "mkdir -p {}".format(output_parent_dir)
	os.system(make_output_parent_dir)

	output_dir = "/hb/scratch/mglasena/test_rna_metrics/" + output_parent_dir + input_fasta_file.split("/")[-1].split(".fas")[0] + "/"
	make_output_dir = "mkdir -p {}".format(output_dir)
	os.system(make_output_dir)

	records = SeqIO.parse(input_fasta_file, "fasta")
	#SeqIO.write(records, output_dir + "consensAlign.ordered.phylip", "phylip-sequential")
	SeqIO.write(records, output_dir + "consensAlign.ordered.phylip", "phylip-sequential")

def get_fasta_alignment_paths_2():
	subdirectory_lst = os.listdir("single_copy_ortholog_fasta_alignments/")
	fasta_file_paths_lst = ["/hb/scratch/mglasena/test_rna_metrics/single_copy_ortholog_fasta_alignments/" + subdir + "/" + "consensAlign.ordered.phylip" for subdir in subdirectory_lst]
	return fasta_file_paths_lst

def reformat_phylip(input_phylip_file):
	input_file = open(input_phylip_file,"r").read().splitlines()
	output_dir = input_phylip_file.split("consens")[0]

	line_1 = input_file[0]
	line_2 = "  ".join(input_file[1].split(" "))
	line_3 = "  ".join(input_file[2].split(" "))
	line_4 = "  ".join(input_file[3].split(" "))
	line_5 = "  ".join(input_file[4].split(" "))

	with open(output_dir + "consensAlign.ordered.phylip","w") as output_file:
		output_file.write(line_1 + "\n")
		output_file.write(line_2 + "\n")
		output_file.write(line_3 + "\n")
		output_file.write(line_4 + "\n")
		output_file.write(line_5 + "\n")

def run_iqtree(fasta_alignment_path):
	run_iqtree = "iqtree -s {} -m MFP".format(fasta_alignment_path)
	os.system(run_iqtree)

def get_gene_paths():
	parent_directory = "/hb/scratch/mglasena/test_rna_metrics/single_copy_ortholog_fasta_alignments/"
	subdirectory_lst = os.listdir(parent_directory)
	gene_paths_lst = [parent_directory + subdir for subdir in subdirectory_lst]
	return gene_paths_lst

def run_paml(gene_path):
	os.chdir(gene_path)
	print("Running paml")
	os.system("codeml")

def main():
	# subset_coverage_dict()

	# bed_file_list = get_zipped_bed_file_list()
	
	# initialize_rna_dict()

	# for regions_file, thresholds_file in bed_file_list:
	# 	for sample in subset_sample_list:
	# 		if sample in regions_file and sample in thresholds_file:
	# 			fill_rna_dict(regions_file, thresholds_file)

	# write_all_rna_dict_csv()
	# filter_rna_dict()
	# get_mrna_gff()
	# create_scaffold_dict()
	# get_passed_rnas()
	# write_new_bed_file()
	# cross_reference()
	# get_gene_ids()

	# gene_ids = get_gene_ids()

	# os.system("mkdir single_gene_gff_records/")
	# Parallel(n_jobs=num_cores)(delayed(make_sco_gff)(gene) for gene in gene_ids)
	
	# # Concatenate all single gene gff records into "sco_gff.gff" file
	# os.system('find ./single_gene_gff_records/ -type f -name "*.record" -exec cat {} \\; > sco_gff.gff')
	
	# # Delete the single gene records
	# os.system('find ./single_gene_gff_records/ -type f -name "*.record" -delete')
	# os.system('rmdir single_gene_gff_records/')

	# run_vcf2fasta()

	# remove_redundant_isoforms()
	# remove_not_multiple_of_three()
	# write_passed_rna_dict_csv()
	# replace_missing_genotype_char()

	#fasta_alignment_paths_list = get_fasta_alignment_paths()
	#Parallel(n_jobs=num_cores)(delayed(convert_fasta_to_phylip)(fasta_file) for fasta_file in fasta_alignment_paths_list)

	fasta_alignment_paths_list_2 = get_fasta_alignment_paths_2()
	#Parallel(n_jobs=num_cores)(delayed(reformat_phylip)(fasta_file) for fasta_file in fasta_alignment_paths_list_2)
	# Parallel(n_jobs=num_cores)(delayed(run_iqtree)(fasta_file) for fasta_file in fasta_alignment_paths_list_2)

	for file in fasta_alignment_paths_list_2:
		file_dir = file.split("consensAlign.ordered.phylip")[0]
		copy = "cp codeml.ctl " + file_dir
		os.system(copy)

	gene_paths_list = get_gene_paths()
	print("Gene paths length: {}".format(len(gene_paths_list)))
	Parallel(n_jobs=num_cores)(delayed(run_paml)(gene) for gene in gene_paths_list)

if __name__ == "__main__":
	main()
