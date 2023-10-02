import os 
from joblib import Parallel, delayed

output_directory = "/hb/scratch/mglasena/phylonet_hmm/phylonet_hmm_variant_sites_vcf/"

reference_genome = "/hb/home/mglasena/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Raw genotype calls file
#genotype_calls = "/hb/home/mglasena/raw_vcf_files/genotype_calls.g.vcf.gz"
genotype_calls_split_multiallelics = output_directory + "genotype_calls_split_multiallelics.g.vcf.gz"

samples_to_include = {
"fragilis_SRR5767279" : "QB3KMK013",
"pallidus_SRR5767285" : "QB3KMK002",
"droebachiensis_SRR5767286" : "QB3KMK014",
"pulcherrimus_SRR5767283" : "QB3KMK016",
}

# Only run once to get split multiallelic file!
def split_multiallelics():
	input_file = genotype_calls
	output_file = output_directory + "genotype_calls_split_multiallelics.g.vcf.gz"
	norm = "bcftools norm -m- -f {} -Oz -o {} {}".format(reference_genome, output_file, input_file)
	os.system(norm)

def separate_SNP_INDEL(variant_type):
	# sample_string = ""
	# for sample in samples_to_include.values():
	# 	sample_string += "-sn " + sample + " "
	# sample_string = sample_string.strip()
	samples = "QB3KMK016,QB3KMK002,QB3KMK013,QB3KMK014"

	input_file = genotype_calls_split_multiallelics
	output_file = output_directory + "genotype_calls_" + variant_type + ".g.vcf.gz"
	#get_variant = 'gatk --java-options "-Djava.io.tmpdir=/hb/scratch/mglasena/test/ -Xms40G -Xmx40" SelectVariants -V {} {} --select-type-to-include SNP --output {}'.format(input_file, sample_string, output_file)
	get_variant = 'bcftools view --samples {} --types {} -Oz -o {} {}'.format(samples, variant_type, output_file, input_file)
	print(get_variant)
	os.system(get_variant)

def filter_variants(variant_type):
	#Alternate option: bcftools filter -e 'QUAL<30 || FS>60 || SOR>3 || MQ<40 || MQRankSum<-12.5 || QD<2 || ReadPosRankSum<-8' -O z -o output input
	if variant_type == "snps":
		input_file = output_directory + "genotype_calls_snps.g.vcf.gz"
		output_file = output_directory + "filtered_snv.g.vcf.gz"
		filter_variants = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'.format(output_file, input_file)

	elif variant_type == "indels":
		input_file = output_directory + "genotype_calls_indels.g.vcf.gz"
		output_file = output_directory + "filtered_indel.g.vcf.gz"
		filter_variants = 'gatk VariantFiltration --output {} --variant {} --verbosity ERROR -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'.format(output_file, input_file)

	os.system(filter_variants)

def merge_vcfs():
	input_snp = output_directory + "filtered_snv.g.vcf.gz"
	input_indel = output_directory + "filtered_indel.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls.g.vcf.gz"
	
	merge_vcfs = "gatk MergeVcfs -I {} -I {} -O {}".format(input_snp, input_indel, output_file)
	os.system(merge_vcfs)
	os.system("rm " + input_snp)
	os.system("rm " + input_indel)

def select_passed_variants():
	filtered_vcf = output_directory + "filtered_genotype_calls.g.vcf.gz"
	output_file = output_directory + "filtered_genotype_calls_pf.g.vcf.gz"
	#select_variants = "gatk SelectVariants -R {} -V {} -O {} --exclude-filtered true".format(reference_genome, filtered_vcf, output_file)
	select_variants = "bcftools view -f .,PASS -Oz -o {} {}".format(output_file, filtered_vcf)
	os.system(select_variants)

# Filter SNPs within 3 base pairs of indel: --SnpGap 3
# Remove monomorphic SNPs where no alternative alleles are called for any of the samples: -e 'AC==0'
# Set individual genotypes with low quality or read depth to missing: -S . -e 'FMT/DP<3 | FMT/GQ<20'
def bcftools_filter():
	input_file = output_directory + "filtered_genotype_calls_pf.g.vcf.gz"
	output_file = output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz"
	filter = '''bcftools filter --SnpGap 3 -e 'AC==0' -Ou {} | bcftools view --types snps | bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' | bcftools norm -m +any -f {} | bcftools filter -S . -e 'ALT="*"' -Oz -o {}'''.format(input_file, reference_genome, output_file)
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

	#separate_SNP_INDEL()
	type_lst = ["snps", "indels"]
	#Parallel(n_jobs=2)(delayed(separate_SNP_INDEL)(variant_type) for variant_type in type_lst)
	
	#filter_variants()
	files = [output_directory + "genotype_calls_snps.g.vcf.gz", output_directory + "genotype_calls_indels.g.vcf.gz"]
	Parallel(n_jobs=2)(delayed(index_vcf)(file) for file in files)
	Parallel(n_jobs=2)(delayed(filter_variants)(variant_type) for variant_type in type_lst)
	merge_vcfs()
	#select_passed_variants()
	#bcftools_filter()
	#index_vcf(output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz")
	#vcf_stats(output_directory + "3bp_filtered_genotype_calls_pf.g.vcf.gz")

if __name__ == "__main__":
	main()