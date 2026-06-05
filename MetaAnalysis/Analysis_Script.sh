#!/bin/bash


#SBATCH --time=70:30:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2000000M   # memory per CPU core
#SBATCH -J "Trofimets"   # job name
#SBATCH --mail-user=cfurnes0@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniforge3/25.3.1-0
eval "$(conda shell.bash hook)"
conda activate rivers

python Riverine_Barriers_SC.py "/home/cfurnes0/MetaAnalysis/Trofimets_et_al2024/Edited_Fastas"
