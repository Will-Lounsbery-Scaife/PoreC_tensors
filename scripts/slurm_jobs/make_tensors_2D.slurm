#!/bin/bash
#SBATCH --job-name=2d
#SBATCH --partition=pe2
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --output=/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/make_tensors_2D_log_slurm_%j.txt

source ~/.bashrc
source activate cool_env

python /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/python_code/make_tensors_2D.py