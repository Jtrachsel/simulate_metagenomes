#!/usr/bin/env python

import os
import random
import copy
import shlex
import pandas as pd
import glob
import subprocess
from Bio import SeqIO
# from sys import argv
import sys

#### need to re-structure this project, 
# extract only the sequence data, and the simulation script
# generate a few different input tables
# then generate a handful of different samples



def calculate_effective_copy_number(row, df):
    if row['seq_type'] == 'chromosome':
        return row['copy_number']

    elif row['seq_type'] == 'plasmid':
        chr_row = df[(df['taxa'] == row['taxa']) & \
                (df['seq_type'] == 'chromosome')]
        return row['copy_number'] * chr_row.iloc[0]['copy_number']

    else:
        raise ValueError

def build_metagenome_fasta(seq_data_folder, metadata, SEED, output_dir):
    # for each taxa in our seq folder amplify and write to a file
    for ix, seq in metadata.iterrows():
        print(f"Amplifying x{seq['actual_copy_number']}: {seq['accession']}")
        in_fp = os.path.join(seq_data_folder, seq['taxa'], seq['sample'],
            seq['seq_type'], seq['accession'] + ".fasta")

        sequences = []
        record_ix = 0
        for record in SeqIO.parse(in_fp, "fasta"):
            record_ix += 1
            for copy_num in range(seq['actual_copy_number']):
                copy_record = copy.deepcopy(record)
                copy_record.id = copy_record.id + f"_{record_ix}_{copy_num}"
                copy_record.description = ""
                sequences.append(copy_record)
        # ADD sample/seed indicator to filename
        metagenome_filename = f'{output_dir}{SEED}_metagenome_fasta.fna'
        with open(metagenome_filename, 'a') as fh:
            SeqIO.write(sequences, fh, "fasta")
    # change return value to output filename
    return metagenome_filename

def select_chromosome_relative_abundance():
    """
    lognormal taxa distribution (mean 1 sigma 2)
    """
    mu = 1
    sigma = 2

    abundance = int(random.lognormvariate(mu, sigma))

    if abundance < 1:
        abundance = 1

    return abundance, 'chromosome'


def select_plasmid_copy_number():
    """
    gamma distribution biased towards lower bound for each
    """
    alpha = 4
    beta = 1
    regimes = {'low': [1, 20, 1],
               'medium': [20, 100, 10],
               'high': [500, 1000, 150]}

    copy_number_regime = random.choice(['low', 'low', 'medium', 'medium', 'high'])
    minimum, maximum, scaling = regimes[copy_number_regime]

    copy_number = random.gammavariate(alpha, beta)
    copy_number = int(copy_number * scaling)

    if copy_number < minimum:
        copy_number = minimum
    elif copy_number > maximum:
        copy_number = maximum

    return copy_number, copy_number_regime


def add_copy_numbers(metadata_fp, SEED, output_dir):
    # accomodate multiple different metadata tables?
    # this would mean you could have simulated metagenomes for different community compositions
    # ie, not all taxa present in all samples, probably more realistic
    
    #simulation_metadata_fp = 'metadata_for_simulation.tsv'
    simulation_metadata_fp = f"{output_dir}{SEED}_{metadata_fp}"
    out_fh = open(simulation_metadata_fp, 'w')
    out_fh.write("taxa\tsample\taccession\tseq_type\tcopy_number\tcopy_regime\n")
    with open(metadata_fp) as fh:
        next(fh) # skip headers
        for line in fh:
            line = line.strip().split('\t')
            if line[3] == 'chromosome':
                relative_abundance, regime = select_chromosome_relative_abundance()
                line.append(str(relative_abundance))
                line.append(regime)
            elif line[3] == 'plasmid':
                copy_number, regime = select_plasmid_copy_number()
                line.append(str(copy_number))
                line.append(regime)

            line = "\t".join(line)
            out_fh.write(line + "\n")
    out_fh.close()
    return simulation_metadata_fp

if __name__ == '__main__':


    # inputs: 
    # seq_data_folder = folder containing sequencing data corresponding to taxa and contig names in the metadata table
    # metadata = filepath to a tab delimited table containing taxa, contig names, and contig classifications
    # seed = random seed to use 
    # output_dir = path to use as an output directory

    if len(sys.argv) < 2:
        print('this script takes 5 positional arguments')
        print('1 = seed')
        print('2 = metadata path')
        print('3 = seq_data directory')
        print('4 = output directory')
        print('5 = coverage (original script used 3.9)')
        sys.exit()
    
    SEED = sys.argv[1]
    metadata = sys.argv[2]
    seq_data_folder = sys.argv[3]
    output_dir = sys.argv[4]
    fcov = sys.argv[5]

    sample_name = metadata.replace('.tsv', '')

    if not os.path.exists(output_dir):
    # Create a new directory because it does not exist 
        os.makedirs(output_dir)
        print("The new directory is created!")
    
    random.seed(SEED)

    # simulate copy numbers
    #'./MAG_gi_plasmid_analysis/data/taxa_metadata.tsv'
    # this function writes a new table with copy numbers added
    # add output folder?
    sim_metadata_fp = add_copy_numbers(metadata, SEED, output_dir)

    # calculate the effective copy number i.e. taxa abundance x plasmid copy number
    metadata = pd.read_csv(sim_metadata_fp, sep='\t')

    # this seems like it could be part of add_copy_numbers()
    # then could write out the full metadata used to build the metagenome
    # actually, it seems like add_copy_numbers() should be changed to operate on pandas dfs instead of iterating through rows....
    metadata['actual_copy_number'] = metadata.apply(\
            lambda x: calculate_effective_copy_number(x, metadata), axis=1)

    metadata.to_csv(sim_metadata_fp)

    #seq_data_folder = './MAG_gi_plasmid_analysis/data/sequences'
    
    # need to change this filepath to contain sample/seed info
    metagenome_fasta_fp = build_metagenome_fasta(seq_data_folder, metadata, SEED, output_dir)

    # fragement length is a bit arbitrary, changed to 350 to match expected frag size with nextera flex library prep
    # changed to hiseq 2x150 PE reads
    # changed fcov from 3.9 to argv[5]
    subprocess.check_call(f'art_illumina --noALN --rndSeed {SEED} --seqSys HS25 --in {metagenome_fasta_fp} --len 150 --fcov {fcov} --out {output_dir}{SEED}_{sample_name} --mflen 350 --sdev 50', shell=True)

