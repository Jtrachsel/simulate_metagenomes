# simulate_metagenomes
Code modified from https://github.com/fmaguire/MAG_gi_plasmid_analysis focusing on simulating metagenomes  

This method helps ensure that genomes carrying plasmids have realistic copy-number dynamics instead of the genome and plasmid sequences being present at 1:1 ratios (which is probably unrealistic).  

## Quickstart  

1. Clone repo  
2. Ensure requirements are met (see env.yml for conda environment specs)  
3. Decompress sequence data  
4. Run simulate_metagenome.py with desired parameters  
    - this script takes 5 positional arguments:  
        - 1 = seed  
        - 2 = metadata path  
        - 3 = seq_data directory  
        - 4 = output directory  
        - 5 = coverage (original script used 3.9)  

    git clone https://github.com/Jtrachsel/simulate_metagenomes.git  
    conda env create -n magsim_lite --file environment.yml  
    conda activate magsim_lite  
    tar -xzvf seq_data/sequences.tar.gz  
    ./simulate_metagenome.py 1 original_community.tsv ./seq_data/ ./output/ 3.9  


## Notes  

- You can change the community proflle by adding or removing organisms both in the metadata (original_community.tsv) and the seq_dat folder.  You will need to provide a plasmid/chromosome classification for each contig, and match the directory structure found in the seq_data directory.  
- Probably best to only use complete genomes.  


