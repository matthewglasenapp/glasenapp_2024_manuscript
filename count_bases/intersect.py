import os

scaffold_lst = ['NW_022145612.1', 'NW_022145594.1', 'NW_022145606.1', 'NW_022145607.1', 'NW_022145611.1', 'NW_022145605.1', 'NW_022145609.1', 'NW_022145613.1', 'NW_022145599.1', 'NW_022145598.1', 'NW_022145600.1', 'NW_022145610.1', 'NW_022145615.1', 'NW_022145601.1', 'NW_022145603.1', 'NW_022145597.1', 'NW_022145604.1', 'NW_022145614.1', 'NW_022145596.1', 'NW_022145602.1', 'NW_022145595.1']

exon_file = "unique_exons.bed"
protein_coding_genes_file = "protein_coding_genes.bed"

coding_bases_all_scaffolds = 82267625

coding_bases_introgressed_90 = 843339
coding_bases_introgressed_80 = 4557598

intron_bases_introgressed_90 = 3461115
intron_bases_introgressed_80 = 18918426

intergenic_bases_introgressed_90 = 4773182
intergenic_bases_introgressed_80 = 124087808

total_bases_21_largest_scaffolds = 828464848

def subset_unique_exons_file():
	unique_exons = gzip.open(exon_file,"r").read().splitlines()
	subset_unique_exons = [item for item in unique_exons if item.split("\t")[0] in scaffold_lst]
	
	with open("subset_unqiue_exons.bed","w") as f:
		for item in subset_unique_exons:
			f.write(item)
			f.write("\n")

def count_exon_bases():
	exon_lst = gzip.open("subset_unqiue_exons.bed","r").read().splitlines()

	base_counter = 0 

	for line in exon_lst:
		exon_length = int(line.split("\t")[2]) - int(line.split("\t")[1])
		base_counter += exon_length

	return base_counter

def subset_protein_coding_genes_files():
	protein_coding_genes = gzip.open(protein_coding_genes_file,"r").read().splitlines()
	subset_protein_coding_genes = [item for item in protein_coding_genes if item.split("\t")[0] in scaffold_lst]
	
	with open("subset_protein_coding_genes","w") as f:
		for item in subset_protein_coding_genes:
			f.write(item)
			f.write("\n")

def count_genic_bases():
	gene_lst = gzip.open("subset_protein_coding_genes","r").read().splitlines()

	base_counter = 0 

	for line in gene_lst:
		gene_length = int(line.split("\t")[2]) - int(line.split("\t")[1])
		base_counter += gene_length

	return base_counter

def get_percent_coding_bases_introgressed():
	subset_unique_exons_file()
	subset_protein_coding_genes_files()

	coding_bases_21_largest_scaffolds = count_exon_bases()
	genic_bases_21_largest_scaffolds = count_genic_bases()
	intron_bases_21_largest_scaffolds = genic_bases_21_largest_scaffolds - coding_bases_21_largest_scaffolds
	intergenic_bases_21_largest_scaffolds = total_bases_21_largest_scaffolds - genic_bases_21_largest_scaffolds

	#print("There are {} total coding bases in the S. purpuratus genome".format(coding_bases_all_scaffolds))
	#print("There are {} coding bases on the 21 largest scaffolds.".format(coding_bases_21_largest_scaffolds))
	#print("{} percent of the total coding bases are on the 21 largest scaffolds.".format((coding_bases_21_largest_scaffolds/coding_bases_all_scaffolds)*100))
	#print("There were {} coding bases declared introgressed at the 90% threshold on the 21 largest scaffolds by PhyloNet_HMM.".format(coding_bases_introgressed_90))
	#print("{} percent of the coding bases on the 21 largest scaffolds were introgressed at the 90% threshold!".format((coding_bases_introgressed_90 / coding_bases_21_largest_scaffolds)*100))

	#print("There were {} coding bases declared introgressed at the 80% threshold on the 21 largest scaffolds by PhyloNet_HMM.".format(coding_bases_introgressed_80))
	print("{} percent of the coding bases on the 21 largest scaffolds were introgressed at the 80% threshold!".format((coding_bases_introgressed_80 / coding_bases_21_largest_scaffolds)*100))

	#print("There are {} genic bases on the 21 largest scaffolds.".format(genic_bases_21_largest_scaffolds))

	#print("There are {} intron bases on the 21 largest scaffolds.".format(intron_bases_21_largest_scaffolds))
	#print("{} percent of the intron bases on the 21 largest scaffolds were introgressed at the 90% threshold!".format((intron_bases_introgressed_90 / intron_bases_21_largest_scaffolds)*100))
	print("{} percent of the intron bases on the 21 largest scaffolds were introgressed at the 80% threshold!".format((intron_bases_introgressed_80 / intron_bases_21_largest_scaffolds)*100))


	#print("There are {} intergenic bases on the 21 largest scaffolds.".format(intergenic_bases_21_largest_scaffolds))
	#print("{} percent of the intergenic bases on the 21 largest scaffolds were introgressed at the 90% threshold!".format((intergenic_bases_introgressed_90 / intergenic_bases_21_largest_scaffolds)*100))
	print("{} percent of the intergenic bases on the 21 largest scaffolds were introgressed at the 80% threshold!".format((intergenic_bases_introgressed_80 / intergenic_bases_21_largest_scaffolds)*100))


def main():
	subset_unique_exons_file()
	get_percent_coding_bases_introgressed()
	
if __name__ == "__main__":
	main()