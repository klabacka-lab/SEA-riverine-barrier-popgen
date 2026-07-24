The script create_alignment.py takes a csv and turns it into an aligned fasta file that can then be used in the PopulationGenetics and EcoEvolity folders.

Make sure to activate the phylo_align conda environment before running the script.

you can create the environment if you haven't yet:
conda config --set channel_priority strict
conda create -n phylo_align -c conda-forge -c bioconda python=3.12 biopython mafft

If you don't have it created, you can do so:

conda env create -f phylo_align.yml
