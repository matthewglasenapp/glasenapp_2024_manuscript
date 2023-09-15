import os
import csv

# Specify the directory that contains protein_coding_genes.bed, unique_CDSs.bed, GenePageGeneralInfo_AllGenes.txt, GeneGoTerms.txt, GeneKoTerms.txt, tu_2012_go.csv, and psg.csv
genome_metadata_dir = "/hb/home/mglasena/dissertation/scripts/phylonet_hmm/genome_metadata/"

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/investigate_tracts/"

# File of unique protein-coding CDSs 
CDS_list_file = genome_metadata_dir + "nonoverlapping_unique_CDS.bed"

psg_tract_LOC_dict = dict()

LOC_coding_base_overlap_dict = dict()

def get_psg_overlapped_tracts():
	psg_list = open(root_dir + "psg_intersect.tsv","r").read().splitlines()[1:]
	
	with open("psg_overlapped_tracts.bed","w") as f:
		for psg in psg_list:
			overlapping_tracts = psg.split("\t")[11].split(",")
			for tract in overlapping_tracts:
				psg_tract_LOC_dict[tract] = psg.split("\t")[0]
				scaffold = tract.split(":")[0].replace("-","_")
				start = tract.split(":")[1].split("-")[0]
				stop = tract.split(":")[1].split("-")[1]
				f.write(scaffold + "\t" + start + "\t" + stop + "\t" + tract + "\n")

# Use bedtools intersect to create file bed file containing the overlap between introgression tracts and protein coding genes 
def intersect_genes(tract_file, gene_file, outfile):
	print("Intersecting {} with {}. Writing results to {}".format(tract_file, gene_file, outfile))
	os.system("bedtools intersect -a " + tract_file + " -b " + gene_file + " -wo > " + outfile)

def get_coding_base_counts():
	overlaps = open("psg_tracts_CDS_overlap.bed","r").read().splitlines()

	for record in overlaps:
		tract = record.split("\t")[3]
		base_overlap = int(record.split("\t")[-1])
		LOC_ID = psg_tract_LOC_dict[tract]
		if LOC_ID in LOC_coding_base_overlap_dict:
			LOC_coding_base_overlap_dict[LOC_ID] += base_overlap
		else:
			LOC_coding_base_overlap_dict[LOC_ID] = base_overlap

	print(LOC_coding_base_overlap_dict)

def update_psg_intersect_csv():
	psg_intersect_file_lines = open(root_dir + "psg_intersect.tsv","r").read().splitlines()[1:]

	csv_output = root_dir + "updated_psg_intersect.tsv"
	csv_file = open(csv_output,"w")
	writer = csv.writer(csv_file, delimiter = "\t")	
	header = ["NCBI Gene ID", "Name", "Synonyms", "Kober and Pogson Gene ID", "Kober and Pogson Name", "Kober and Pogson Synonyms", "PSG #", "Length", "Introgressed Bases", "Percent Bases Introgressed", "Coordinates", "Overlapping introgression tract(s)", "Sdro", "Sfra", "Spal", "Hpul", "Overlapping Coding Bases", "Percent Coding"]
	writer.writerow(header)

	for line in psg_intersect_file_lines:
		LOC_ID = line.split("\t")[0]
		if LOC_ID in LOC_coding_base_overlap_dict:
			base_overlap = LOC_coding_base_overlap_dict[LOC_ID]
			percent_coding = float(base_overlap / int(line.split("\t")[7]) * 100) 
			data = line.split("\t")
			data.append(base_overlap)
			data.append(percent_coding)
			writer.writerow(data)
		else:
			data = line.split("\t")
			data.append("0")
			data.append("0.0")
			writer.writerow(data)

	csv_file.close()

def main():
	get_psg_overlapped_tracts()
	intersect_genes("psg_overlapped_tracts.bed", CDS_list_file, "psg_tracts_CDS_overlap.bed")
	get_coding_base_counts()
	update_psg_intersect_csv()

if __name__ == "__main__":
	main()
