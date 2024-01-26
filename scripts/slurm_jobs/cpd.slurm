#!/bin/bash
#SBATCH --job-name=cpd
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=6
#SBATCH --output=/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/cpd_log_slurm_%j.txt

source ~/.bashrc
source activate tfx_env

python /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/python_code/cpd.py
