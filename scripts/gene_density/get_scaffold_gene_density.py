protein_coding_genes_file = "/Users/matt/Documents/Github/dissertation_chapter_2/data/spur_genome_metadata/protein_coding_genes.bed"
scaffold_lengths_file = "/Users/matt/Documents/Github/dissertation_chapter_2/scripts/pixy_dxy/scaffolds.tsv"

scaffold_lengths_dict = dict()
scaffold_gene_dict = dict()

def create_scaffold_lengths_dict():
	scaffold_lengths = open(scaffold_lengths_file,"r").read().splitlines()
	for line in scaffold_lengths:
		scaffold_lengths_dict[line.split("\t")[0]] = int(line.split("\t")[1])

def create_scaffold_gene_dict():
	protein_coding_genes = open(protein_coding_genes_file,"r").read().splitlines()
	scaffolds = [scaffold.split("\t")[0] for scaffold in open(scaffold_lengths_file).read().splitlines()]

	for scaffold in scaffolds:
		scaffold_gene_dict[scaffold] = 0
	
	for gene in protein_coding_genes:
		scaffold = gene.split("\t")[0]
		if scaffold in scaffold_gene_dict:
			gene_length = int(gene.split("\t")[2]) - int(gene.split("\t")[1])
			scaffold_gene_dict[scaffold] += gene_length

def write_results():
	with open("scaffold_gene_density.tsv","w") as f:
		for key,value in scaffold_gene_dict.items():
			f.write(key + "\t" + str(value / scaffold_lengths_dict[key]) + "\n")

def main():
	create_scaffold_lengths_dict()
	create_scaffold_gene_dict()
	write_results()

if __name__ == "__main__":
	main()