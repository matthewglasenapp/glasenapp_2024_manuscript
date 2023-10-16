import os 
from joblib import Parallel, delayed

output_directory = "/hb/home/mglasena/test_vcf/combined_vcf_files/"

reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Raw genotype calls file
genotype_calls = output_directory + "genotype_calls.g.vcf.gz"

# Only run once to get split multiallelic file!
def split_multiallelics():
	input_file = genotype_calls
	output_file = output_directory + "genotype_calls_split_multiallelics.g.vcf.gz"
	norm = "bcftools norm -m- -f {} -Oz -o {} {}".format(reference_genome, output_file, input_file)
	os.system(norm)

def separate_SNP_INDEL(variant_type):
	input_file = output_directory + "genotype_calls_split_multiallelics.g.vcf.gz"
	output_file = output_directory + "genotype_calls_" + variant_type + ".g.vcf.gz"
	get_variant = 'gatk SelectVariants -V {} --select-type-to-include {} --output {}'.format(input_file, variant_type, output_file)
	print(get_variant)
	os.system(get_variant)

def filter_variants(variant_type):
	#Alternate option: bcftools filter -e 'QUAL<30 || FS>60 || SOR>3 || MQ<40 || MQRankSum<-12.5 || QD<2 || ReadPosRankSum<-8' -O z -o output input
	input_file = output_directory + "genotype_calls_" + variant_type + ".g.vcf.gz"
	output_file = output_directory + "filtered_" + variant_type + ".g.vcf.gz"

	if variant_type == "SNP":
		filter_variants = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'.format(output_file, input_file)

	elif variant_type == "INDEL":
		filter_variants = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'.format(output_file, input_file)

	os.system(filter_variants)

def merge_vcfs():
	input_snp = output_directory + "filtered_SNP.g.vcf.gz"
	input_indel = output_directory + "filtered_INDEL.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	
	merge_vcfs = "gatk MergeVcfs -I {} -I {} -O {}".format(input_snp, input_indel, output_file)
	os.system(merge_vcfs)
	#os.system("rm " + input_snp)
	#os.system("rm " + input_indel)

def select_passed_variants():
	filtered_vcf = output_directory + "filtered_genotype_calls.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls_pf.g.vcf.gz"
	select_variants = "gatk SelectVariants -R {} -V {} -O {} --exclude-filtered true".format(reference_genome, filtered_vcf, output_file)
	#select_variants = "bcftools view -f .,PASS -Oz -o {} {}".format(output_file, filtered_vcf)
	os.system(select_variants)

# Filter SNPs within 3 base pairs of indel: --SnpGap 3
# Remove monomorphic SNPs where no alternative alleles are called for any of the samples: -e 'AC==0'
# Set individual genotypes with low quality or read depth to missing: -S . -e 'FMT/DP<3 | FMT/GQ<20'
def bcftools_filter():
	input_file = output_directory + "filtered_genotype_calls_pf.g.vcf.gz"
	output_file = output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz"
	#filter = '''bcftools filter --SnpGap 3 -e 'AC==0' -Ou {} | bcftools view --types snps | bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' | bcftools norm -m +any -f {} | bcftools filter -S . -e 'ALT="*"' -Oz -o {}'''.format(input_file, reference_genome, output_file)
	filter = '''bcftools filter --SnpGap 3 -e 'AC==0' -Ou {} | bcftools view --types snps | bcftools filter -e 'FMT/DP<3 | FMT/GQ<20' | bcftools norm -m +any -f {} | bcftools filter -e 'ALT="*"' -Oz -o {}'''.format(input_file, reference_genome, output_file)
	os.system(filter)
	#os.system("rm " + input_file)

def index_vcf(input_file):
	index = "gatk IndexFeatureFile -I {}".format(input_file)
	os.system(index)

def vcf_stats(input_file):
	get_samples_file = "bcftools query -l {} > samples_file.txt".format(input_file)
	os.system(get_samples_file)
	get_stats = "bcftools stats --samples-file samples_file.txt {}".format(input_file)
	os.system(get_stats)
	os.system("rm samples_file.txt")

def main():
	# Leave commented out unless running for the first time
	#split_multiallelics()
	#index_vcf(output_directory + "genotype_calls_split_multiallelics.g.vcf.gz")

	type_lst = ["SNP", "INDEL"]
	#Parallel(n_jobs=2)(delayed(separate_SNP_INDEL)(variant_type) for variant_type in type_lst)
	#Parallel(n_jobs=2)(delayed(filter_variants)(variant_type) for variant_type in type_lst)
	#merge_vcfs()
	#select_passed_variants()
	bcftools_filter()
	index_vcf(output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz")
	
	#vcf_stats(output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz")

if __name__ == "__main__":
	main()