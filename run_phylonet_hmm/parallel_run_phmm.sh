#!/bin/bash
#SBATCH --job-name=run_phmm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mglasena@ucsc.edu
#SBATCH --output=run_phmm_%j.out
#SBATCH --nodes=1
#SBATCH --partition=128x24
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=2-0

module load java/java-8
module load python-3.6.5

python3 run_phmm.py