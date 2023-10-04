import os
import csv
import subprocess

# Specify the directory that contains protein_coding_genes.bed, unique_CDSs.bed, GenePageGeneralInfo_AllGenes.txt, GeneGoTerms.txt, GeneKoTerms.txt, tu_2012_go.csv, and psg.csv
genome_metadata_dir = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/"

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_80/investigate_tracts/"

# File of unique protein-coding CDSs 
merged_CDS_file = genome_metadata_dir + "nonoverlapping_unique_CDS.bed"

# Unmerged CDS file
full_CDS_file = genome_metadata_dir + "unique_CDS.bed"

# Overlapped genes with significant branch sites tests
Sdro_bsites_overlap = root_dir + "Sdro_bsites_intersect.tsv"
Spal_bsites_overlap = root_dir + "Spal_bsites_intersect.tsv"

Sdro_bsites_tract_LOC_dict = dict()
Spal_bsites_tract_LOC_dict = dict()

Sdro_LOC_coding_base_overlap_dict = dict()
Spal_LOC_coding_base_overlap_dict = dict()

def get_psg_overlapped_tracts(overlap_file):
	sample = overlap_file.split("/")[-1].split("_")[0]
	psg_list = open(overlap_file,"r").read().splitlines()[1:]
	output_file_name = overlap_file.split("/")[-1].split("_intersect.tsv")[0]
	
	with open(root_dir + output_file_name + "_overlapped_tracts.bed","w") as f:
		for psg in psg_list:
			overlapping_tracts = psg.split("\t")[12].split(",")
			
			for tract in overlapping_tracts:
				LOC_ID = psg.split("\t")[0]
				gene_coordinates = psg.split("\t")[11]
				
				if sample == "Sdro":
					Sdro_bsites_tract_LOC_dict[tract] = LOC_ID
				elif sample == "Spal":
					Spal_bsites_tract_LOC_dict[tract] = LOC_ID
				
				scaffold = tract.split(":")[0].replace("-","_")
				start = tract.split(":")[1].split("-")[0]
				stop = tract.split(":")[1].split("-")[1]
				f.write(scaffold + "\t" + start + "\t" + stop + "\t" + tract + "\t" + LOC_ID + "_" + gene_coordinates + "\n")

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file, gene_file, outfile):
	print("Intersecting {} with {}. Writing results to {}".format(root_dir + tract_file, gene_file, root_dir + outfile))
	os.system("bedtools intersect -a " + root_dir + tract_file + " -b " + gene_file + " -wo > " + root_dir + outfile)

def get_coding_base_counts(CDS_overlap_file):
	sample = CDS_overlap_file.split("/")[-1].split("_")[0]
	overlaps = open(CDS_overlap_file,"r").read().splitlines()

	for record in overlaps:
		tract = record.split("\t")[3]
		base_overlap = int(record.split("\t")[-1])

		if sample == "Sdro":
			LOC_ID = Sdro_bsites_tract_LOC_dict[tract]
			# Bedtools merge gets rid of identifer info columns. Make sure that the overlapped CDS is actually within the LOC coordinates
			# Sometimes an introgression tract spans CDS from more than one gene, and that CDS was getting double-counted
			CDS_start = int(record.split("\t")[6])
			CDS_end = int(record.split("\t")[7])
			LOC_start = int(record.split("\t")[4].split(":")[1].split("-")[0])
			LOC_end = int(record.split("\t")[4].split(":")[1].split("-")[1])
			
			if LOC_ID in Sdro_LOC_coding_base_overlap_dict and CDS_start >= LOC_start and CDS_end <= LOC_end:
				Sdro_LOC_coding_base_overlap_dict[LOC_ID] += base_overlap
			elif CDS_start >= LOC_start and CDS_end <= LOC_end:
				Sdro_LOC_coding_base_overlap_dict[LOC_ID] = base_overlap

		elif sample == "Spal":
			LOC_ID = Spal_bsites_tract_LOC_dict[tract]

			# Bedtools merge gets rid of identifer info columns. Make sure that the overlapped CDS is actually within the LOC coordinates
			# Sometimes an introgression tract spans CDS from more than one gene, and that CDS was getting double-counted
			CDS_start = int(record.split("\t")[6])
			CDS_end = int(record.split("\t")[7])
			LOC_start = int(record.split("\t")[4].split(":")[1].split("-")[0])
			LOC_end = int(record.split("\t")[4].split(":")[1].split("-")[1])


			if LOC_ID in Spal_LOC_coding_base_overlap_dict and CDS_start >= LOC_start and CDS_end <= LOC_end:
				Spal_LOC_coding_base_overlap_dict[LOC_ID] += base_overlap
			elif CDS_start >= LOC_start and CDS_end <= LOC_end:
				Spal_LOC_coding_base_overlap_dict[LOC_ID] = base_overlap

def update_psg_intersect_csv(original_intersect_file):
	sample = original_intersect_file.split("/")[-1].split("_")[0]
	psg_intersect_file_lines = open(original_intersect_file,"r").read().splitlines()[1:]

	csv_output = root_dir + "updated_" + original_intersect_file.split("/")[-1]
	csv_file = open(csv_output,"w")
	writer = csv.writer(csv_file, delimiter = "\t")	
	header = ["NCBI Gene ID", "Name", "Synonyms", "Kober and Pogson Gene ID", "Kober and Pogson Name", "Kober and Pogson Synonyms", "PSG ?", "PSG Rank", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Overlapping introgression tract(s)", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Coding Bases", "Total Coding Bases", "Percent Coding Bases Introgressed"]
	writer.writerow(header)

	for line in psg_intersect_file_lines:
		LOC_ID = line.split("\t")[0]
		
		if sample == "Sdro":
			if LOC_ID in Sdro_LOC_coding_base_overlap_dict:
				base_overlap = Sdro_LOC_coding_base_overlap_dict[LOC_ID]

				# Deal with weird edge case of different case across different file types for eef1G, six1
				if LOC_ID == "eef1g":
					LOC_ID = "eef1G"
				if LOC_ID == "six1":
					LOC_ID = "Six1"
					
				# Get total coding base count for LOC
				count_LOC_coding_bases = "cat {} | grep '{}' | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum}}'".format(full_CDS_file, LOC_ID)
				total_coding_bases_LOC = int(subprocess.check_output(count_LOC_coding_bases, shell=True, universal_newlines=True))

				percent_coding = float((base_overlap / total_coding_bases_LOC) * 100) 
				data = line.split("\t")
				data.append(base_overlap)
				data.append(total_coding_bases_LOC)
				data.append(percent_coding)
				writer.writerow(data)
			
			else:
				data = line.split("\t")
				data.append("0")
				data.append("0.0")
				data.append("0.0")
				writer.writerow(data)
		
		elif sample == "Spal":		
			if LOC_ID in Spal_LOC_coding_base_overlap_dict:
				base_overlap = Spal_LOC_coding_base_overlap_dict[LOC_ID]

				# Deal with weird edge case of different case across different file types for eef1G, six1
				if LOC_ID == "eef1g":
					LOC_ID = "eef1G"
				if LOC_ID == "six1":
					LOC_ID = "Six1"
					
				# Get total coding base count for LOC
				count_LOC_coding_bases = "cat {} | grep '{}' | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum}}'".format(full_CDS_file, LOC_ID)
				total_coding_bases_LOC = int(subprocess.check_output(count_LOC_coding_bases, shell=True, universal_newlines=True))

				percent_coding = float((base_overlap / total_coding_bases_LOC) * 100) 
				data = line.split("\t")
				data.append(base_overlap)
				data.append(total_coding_bases_LOC)
				data.append(percent_coding)
				writer.writerow(data)
			
			else:
				data = line.split("\t")
				data.append("0")
				data.append("0.0")
				data.append("0.0")
				writer.writerow(data)
		
	csv_file.close()

def main():
	os.chdir(root_dir)
	
	get_psg_overlapped_tracts(Sdro_bsites_overlap)
	get_psg_overlapped_tracts(Spal_bsites_overlap)
	
	#intersect_genes("Sdro_bsites_overlapped_tracts.bed", merged_CDS_file, "Sdro_bsites_CDS_overlap.bed")
	#intersect_genes("Spal_bsites_overlapped_tracts.bed", merged_CDS_file, "Spal_bsites_CDS_overlap.bed")

	get_coding_base_counts(root_dir + "Spal_bsites_CDS_overlap.bed")
	get_coding_base_counts(root_dir + "Sdro_bsites_CDS_overlap.bed")
	
	update_psg_intersect_csv(Sdro_bsites_overlap)
	update_psg_intersect_csv(Spal_bsites_overlap)

if __name__ == "__main__":
	main()
