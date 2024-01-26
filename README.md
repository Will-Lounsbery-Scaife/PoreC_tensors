## Tensor-Based Analysis of High-Order Chromatin Interactions in Pore-C Data

### Initial setup:
 - Use the following commands to install a Conda/Miniconda env for running TensorFox 
 - https://github.com/felipebottega/Tensor-Fox for more details
 
```bash
pip install TensorFox
conda create --name tfx_env --channel defaults jupyter numpy pandas scipy scikit-learn matplotlib numba IPython sparse_dot_mkl  
```

### Protocol for creating and factorizing high-dimensional tensors (i.e. not 2D)
Slurm scripts are in '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/slurm_jobs' \
The corresponding Python scripts are in '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/python_code' \
```bash
cd /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/slurm_jobs
```

To create tensors for Pore-C data:
```bash
sbatch make_tensors.slurm
```
 - You can check the log file here: '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/make_tensors_log_slurm_{job_number}.txt
 - The tensors are saved in '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/tensors/'
 - You can change the parameters at the top of the make_tensors.py script
 - You may have to modify the memory allocation at the top of the slurm script
    - Using a smaller window size will require more memory
    - Creating tensors for larger chromosomes will require more memory
    - Using a higher maximum cardinality will require *way* more memory

Make sure that make_tensors.slurm finishes before running the CPD job. \
To perform CPD on the tensors:
```bash
sbatch cpd.slurm
```
 - You can check the log file here: '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/cpd_log_slurm_{job_number}.txt
 - The factors will be saved in '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/factors/'

 To plot the factors, you can use the jupyter notebook 'scripts/jupyter_notebooks/cpd_results.ipynb/cpd_results.ipynb'


### Protocol for creating and decomposing 2D tensors (i.e. Pore-C matrices with max cardinality = 2)
First, you will need to create a new conda environment with PyTorch, Cooler, CoolBox, and CoolTools. \
For more information, refer to: 
 - https://gangcaolab.github.io/CoolBox/installation.html
 - https://cooler.readthedocs.io/en/latest/quickstart.html#installation
 - https://cooltools.readthedocs.io/en/latest/#

To create your conda env, use the following commands:
```bash
conda create -n cool_env python=3.10.12
conda activate cool_env
conda install -c conda-forge -c bioconda cooler cooltools bioframe coolbox numpy pandas h5py scipy
conda install pytorch=2.1.1 torchvision=0.15.2 torchaudio=2.1.1 cudatoolkit=11.3.1 -c pytorch
```
Note that these commands may take a while to run. \

To create the 2D Pore-C tensors and create .cool files from the tensors:
```bash
cd /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/slurm_jobs
sbatch make_tensors_2D.slurm
# wait till make_tensors_2D is done running
sbatch 2D_tensors_to_cools.slurm
```
Note that the 2D_tensors_to_cools script will create a balanced .cool file and and unbalanced .cool file for each chromosome \

Then, you can play with the Jupyter Notebooks to create Pore-C matrices, decompose them, and annotate them with histone/CTCF data.
```bash
cd /commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/scripts/jupyter_notebooks
```

 - coolbox_plots.ipynb: uses CoolBox to create Pore-C matrices
    - Includes tracks for histone modification and CTCF enrichment
    - You will have to change the subrange depending on which chromosome you are plotting
 - cooltools_plots.ipynb: uses CoolTools to perform Eigendecomposition on Pore-C matrices
    - You will have to modify the code depending on the chromosome you are plotting
 
