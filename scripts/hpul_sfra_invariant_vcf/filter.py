import os 

reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Raw genotype calls file. Leave this commented out. 
genotype_calls = "/hb/scratch/mglasena/phylonet_hmm/hpul_sfra_invariant_sites_vcf/combined_vcf_files/genotype_calls.g.vcf.gz"
genotype_calls_split_multiallelics = "/hb/scratch/mglasena/phylonet_hmm/hpul_sfra_invariant_sites_vcf/combined_vcf_files/genotype_calls_split_multiallelics.g.vcf.gz"

output_directory = "/hb/scratch/mglasena/dxy/combined_vcf_files/"

samples_to_include = {
"fragilis_SRR5767279" : "QB3KMK013",
"pulcherrimus_SRR5767283" : "QB3KMK016",
}

# Only run once to get split multiallelic file!
def split_multiallelics():
	input_file = genotype_calls
	output_file = output_directory + "genotype_calls_split_multiallelics.g.vcf.gz"
	norm = "bcftools norm -m- -f {} -Oz -o {} {}".format(reference_genome, output_file, input_file)
	os.system(norm)

def separate_variant_nonvariant():
	sample_string = ""
	for sample in samples_to_include.values():
		sample_string += "-sn " + sample + " "
	sample_string = sample_string.strip()
	
	input_file = genotype_calls_split_multiallelics
	output_variant = output_directory + "genotype_calls_variant.g.vcf.gz"
	output_nonvariant = output_directory + "genotype_calls_nonvariant.g.vcf.gz"
	get_variant = "gatk SelectVariants -V {} {} --select-type-to-include SNP --exclude-non-variants --output {}".format(input_file, sample_string, output_variant)
	get_nonvariant = "gatk SelectVariants -V {} {} --select-type-to-include NO_VARIATION --output {}".format(input_file, sample_string, output_nonvariant)
	os.system(get_variant)
	os.system(get_nonvariant)

def filter_variants():
	input_variant = output_directory + "genotype_calls_variant.g.vcf.gz"
	output_variant = output_directory + "filtered_variant.g.vcf.gz"

	filter_variant = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'.format(output_variant, input_variant)

	os.system(filter_variant)
	#os.system("rm " + input_variant)

def select_passed_variants():
	filtered_vcf = output_directory + "filtered_variant.g.vcf.gz"
	output_file = output_directory + "filtered_variant_pf.g.vcf.gz"
	select_variants = "gatk SelectVariants -R {} -V {} -O {} --exclude-filtered true".format(reference_genome, filtered_vcf, output_file)
	os.system(select_variants)
	
def merge_vcfs():
	input_variant = output_directory + "filtered_variant_pf.g.vcf.gz"
	input_nonvariant = output_directory + "genotype_calls_nonvariant.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	
	merge_vcfs = "gatk MergeVcfs -I {} -I {} -O {}".format(input_variant, input_nonvariant, output_file)
	os.system(merge_vcfs)
	os.system("rm " + input_variant)
	os.system("rm " + input_nonvariant)

# Set individual genotypes with low quality or read depth to missing: -S . -e 'FMT/DP<3 | FMT/GQ<20'
# Join multiallelics
def bcftools_filter():
	input_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls_individual_genotypes.g.vcf.gz"
	filter = "bcftools filter -S . -e 'FMT/DP<8 | FMT/GQ<30 | FMT/RGQ<30' -Ou {} | bcftools filter -e 'F_MISSING > 0.5' -Ou | bcftools norm -m +any -f {} -Oz -o {}".format(input_file, reference_genome, output_file)
	os.system(filter)
	os.system("rm " + input_file)

def index_vcf(input_file):
	index = "gatk IndexFeatureFile -I {}".format(input_file)
	os.system(index)

def main():
	# Leave commented out unless running for the first time
	split_multiallelics()
	index_vcf(output_directory + "genotype_calls_split_multiallelics.g.vcf.gz")
	
	separate_variant_nonvariant()
	filter_variants()
	select_passed_variants()
	merge_vcfs()
	bcftools_filter()
	index_vcf(output_directory + "filtered_genotype_calls_individual_genotypes.g.vcf.gz")

if __name__ == "__main__":
	main()