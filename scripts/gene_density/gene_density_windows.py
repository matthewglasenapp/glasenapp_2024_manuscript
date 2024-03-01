import csv

output_file = "/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/gene_density/chromosome_window_density.csv"

protein_coding_genes_file = "/Users/matt/Documents/Github/glasenapp_2024_manuscript/data/spur_genome_metadata/protein_coding_genes.bed"
scaffold_lengths_file = "/Users/matt/Documents/Github/glasenapp_2024_manuscript/scripts/pixy_dxy/scaffolds.tsv"

scaffold_lengths_dict = dict()
scaffold_gene_dict = dict()
scaffold_windows_dict = dict()

window_size = 1000000

def create_scaffold_lengths_dict():
	scaffold_lengths = open(scaffold_lengths_file,"r").read().splitlines()
	for line in scaffold_lengths:
		scaffold_lengths_dict[line.split("\t")[0]] = int(line.split("\t")[1])

def create_scaffold_gene_dict():
	protein_coding_genes = open(protein_coding_genes_file,"r").read().splitlines()
	scaffolds = [scaffold.split("\t")[0] for scaffold in open(scaffold_lengths_file).read().splitlines()]
	
	for gene in protein_coding_genes:
		scaffold = gene.split("\t")[0]
		stop = int(gene.split("\t")[2])
		start = int(gene.split("\t")[1])
		if scaffold in scaffold_gene_dict:
			scaffold_gene_dict[scaffold].append([start, stop])
		else:
			scaffold_gene_dict[scaffold] = [[start, stop]]

def create_scaffold_windows_dict():
	for key,value in scaffold_lengths_dict.items():
		scaffold = key
		scaffold_windows_dict[scaffold] = []
		length = value 

		# Account for the zero indexing of BED files
		# Example: Scaffold NW_022145612.1 is 62154208 bases
		# The BED coordinates would be 0, 62154207
		i = 0
		while i < length - 1:
			window_start = i
			window_stop = i + window_size - 1
			if window_stop >= length - 1:
				window_stop = length - 1

			window = [window_start, window_stop]
			scaffold_windows_dict[scaffold].append(window)
			i = window_stop + 1

def calculate_gene_density_by_window():
	for key,value in scaffold_windows_dict.items():
		scaffold = key
		window_lst = value
		gene_lst = scaffold_gene_dict[scaffold]

		for window in window_lst:
			window_start = window[0]
			window_stop = window[1]
			genic_bases = 0
			
			for gene in gene_lst:
				gene_start = gene[0]
				gene_stop = gene[1]
				
				if gene_start >= window_start and gene_start <= window_stop:
					
					if gene_stop <= window_stop:
						genic_bases += gene_stop - gene_start + 1
					
					if gene_stop > window_stop:
						genic_bases += window_stop - gene_start + 1

				elif gene_start < window_start and gene_stop >= window_start:
					if gene_stop == window_start:
						genic_bases += 1

					if gene_stop <= window_stop:
						genic_bases += gene_stop - window_start + 1

					if gene_stop > window_stop:
						genic_bases += window_stop - window_start 


			window.append(genic_bases)
			gene_density = genic_bases / (window_stop - window_start + 1)
			window.append(gene_density)

def write_to_output_file():
	output_csv = open(output_file, "w")
	writer = csv.writer(output_csv)
	header = ["chromosome", "window_start", "window_stop", "gene_density"]
	for key,value in scaffold_windows_dict.items():
		scaffold = key
		for value in value:
			window_start = value[0]
			window_stop = value[1]
			gene_density = value[3]
			data = [scaffold, window_start, window_stop, gene_density]
			writer.writerow(data)
	output_csv.close()

def main():
	create_scaffold_lengths_dict()
	create_scaffold_gene_dict()
	create_scaffold_windows_dict()
	calculate_gene_density_by_window()
	write_to_output_file()

if __name__ == "__main__":
	main()