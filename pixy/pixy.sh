#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=pixy_%J.out
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-0 

vcf_file="/hb/scratch/mglasena/dxy/vcf_files/combined_vcf_files/filtered_genotype_calls_individual_genotypes.g.vcf.gz"
pop_file="popfile.txt"
bed_file="ten_kb_tracts.bed"
cores="24"
output_prefix="introgression_tracts"

pixy --stats dxy --vcf $vcf_file --populations $pop_file --bed_file $bed_file --output_prefix $output_prefix --n_cores $cores

python3 -u boostrap_pixy.py