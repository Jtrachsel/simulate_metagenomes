#!/bin/bash

#SBATCH --job-name="simulate metagenomes"              # name of the job submitted
#SBATCH -p scavenger                                   # name of the queue you are submitting to
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=50G                                      
#SBATCH -t 48:00:00                                    # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=ENTER_EMAIL_HERE                   # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                     # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                              # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                              # optional but it prints our standard error

# ENTER COMMANDS HERE:
module load miniconda

source activate /project/fsepru/conda_envs/magsim_lite

./simulate_metagenome.py 1 original_community.tsv seq_dat output/ 2.9
./simulate_metagenome.py 2 original_community.tsv seq_dat output/ 2.9
./simulate_metagenome.py 3 original_community.tsv seq_dat output/ 2.9
./simulate_metagenome.py 4 original_community.tsv seq_dat output/ 2.9
./simulate_metagenome.py 5 original_community.tsv seq_dat output/ 2.9

parallel 'gzip {}' ::: ./output/*.fq
parallel 'gzip {}' ::: ./output/*.fna