import os
import csv
import math 

root_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/introgressed_genes/"

fasta_alignment_dir = "/hb/scratch/mglasena/phylonet_hmm/process_hmm_90/paml/introgressed_genes/single_copy_ortholog_fasta_alignments/"

stop_codon_dict = dict()
dN_dS_dict = dict()

stop_codons = ["TAA", "TAG", "TGA"]

def get_gene_paths():
    subdirectory_lst = os.listdir(fasta_alignment_dir)
    gene_paths_lst = [fasta_alignment_dir + subdir for subdir in subdirectory_lst]
    return gene_paths_lst

def find_stop_codons(dir):
    file = dir + "/consensAlign.ordered.phylip"
    gene = dir.split("/")[-1]
    lines = open(file,"r").read().splitlines()[1:]
    seq_dict = dict()
    for line in lines:
        name = line.split("  ")[0]
        seq = line.split("  ")[1]
        seq_dict[name] = seq

    for sample in seq_dict:
        seq = seq_dict[sample]
        sequence_length = len(seq)
        for i in range(0, len(seq), 3):
            if seq[i:i+3] in stop_codons:
                if i == sequence_length - 3:
                    break

                if gene in stop_codon_dict: 
                    if sample in stop_codon_dict[gene]:
                        stop_codon_dict[gene][sample].append(i+1)
                    else:
                        stop_codon_dict[gene][sample] = [i+1]

                else:
                    stop_codon_dict[gene] = {sample: [i+1]}



def extract_paml(gene):
    with open(gene + '/paml_out.txt') as f:
        lines = f.readlines() # list containing lines of file
        i = 1
        count = 0
        
        for line in lines:
            line = line.strip()
            
            if "omega (dN/dS) =" in line:
                omega = float(line.split()[-1])
                count += 1
            
            if "tree length for dN:" in line:
                dn = float(line.split()[-1])
                count += 1
            
            if "tree length for dS:" in line:
                ds = float(line.split()[-1])
                count += 1
        
            if count == 3 and math.isnan(omega) != True:
                dN_dS_dict[gene.split('/')[-1]] = [dn, ds, omega]
                count = 0
                break

    f.close()

def write_data_to_csv():
    output_tsv_file = open("dNdS.tsv", "w")
    writer = csv.writer(output_tsv_file, delimiter="\t")    
    header = ["gene", "dN", "dS", "dN/dS"]
    writer.writerow(header)

    for key, value in dN_dS_dict.items():
        if not key in stop_codon_dict: 
            data_to_write = [key, value[0], value[1], value[2]]
            writer.writerow(data_to_write)

    output_tsv_file.close()

def main():
    os.chdir(root_dir)

    gene_paths_list = get_gene_paths()

    for gene in gene_paths_list:
        find_stop_codons(gene)
        extract_paml(gene)

    write_data_to_csv()

if __name__ == "__main__":
    main()
