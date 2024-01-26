#!/bin/bash
#SBATCH --job-name=tens
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --output=/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/make_tensors_log_slurm_%j.txt

source ~/.bashrc
source activate tfx_env

python /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/python_code/make_tensors.py
