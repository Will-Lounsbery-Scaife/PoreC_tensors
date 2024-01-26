import numpy as np
import TensorFox as tfx
import pickle
import numba
import os
import sys
import time


####################
### PARAMETERS #####
####################

# The CPD ranks to use
# if you run 'ranks = [3,4], then the script will run CPD with rank 3 (creating 3 factors) as well as rank 4 (creating 4 factors)
ranks = [3,4,5,6,7] 

# The parameter sets to use. With 9 different parameter sets, the script will run CPD 9 times for each rank.
# You will probably want to experiment with these a lot -- the ones I chose were arbitrary
param_set = [1,2,3,4,5,6,7,8,9]

# A more exhausting search for the optimal configuration (rank and parameters) is probably warranted
# To change the actual parameters, see the cpd_fact function below, and the TensorFox documentation:
# https://github.com/felipebottega/Tensor-Fox/blob/master/tutorial/2-basic_options.ipynb
# https://github.com/felipebottega/Tensor-Fox/blob/master/tutorial/3-intermediate_options.ipynb
# https://github.com/felipebottega/Tensor-Fox/blob/master/tutorial/4-advanced_options.ipynb 


# Which chromosome tensors to perform CPD on
chromosomes = ["chr21"]
'''
chromosomes = [
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", 
    "chr21", "chr22"
]
'''

####################
####################
####################

print("starting script!")

# start timer
start = time.time()

tensor_directory = "/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/tensors/"
output_directory = "/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output"
factor_directory = os.path.join(output_directory, "factors/")


# Define a wrapper around the original function to always set threads to 2
original_set_num_threads = numba.set_num_threads
def patched_set_num_threads(n):
    original_set_num_threads(4)
numba.set_num_threads = patched_set_num_threads

def cpd_fact(chr, cpd_rank, param):
    class options:
        display = 1
        symm = True
    options = tfx.make_options(options)

    # Load the tensor from the saved pickle file
    with open(os.path.join(tensor_directory, f"NlaIII_GM12878_{chr}_tensor.pickle"), "rb") as f:
        T = pickle.load(f)

    # Define the rank of the tensor decomposition
    R = cpd_rank

    if True:
        if param == 1:
            print("param set = 1")

        elif param == 2:
            print("param set = 2")
            options.initialization = 'smart_random'

        elif param == 3:
            print("param set = 3")
            options.initialization = 'smart'

        elif param == 4:
            print("param set = 4")
            options.method = 'als'

        elif param == 5:
            print("param set = 5")
            options.method = 'als'
            options.inner_method = 'cg_static'

        elif param == 6:
            print("param set = 6")
            options.tol = 1e-20
            options.max_iter = 300
        
        elif param == 7:
            print("param set = 7")
            options.inner_method = 'direct'
        
        elif param == 8:
            print("param set = 8")
            options.inner_method = 'cg_static'
        
        elif param == 9:
            print("param set = 9")
            options.inner_method = 'als'

    # Decompose the tensor using CPD
    factors, output = tfx.foxit(T, R, options)

    # save the factors to a pickle file
    with open(os.path.join(factor_directory, f"{chr}_rank{cpd_rank}_paramset{param}_factors_100kb.pickle"), "wb") as f:
        pickle.dump(factors, f)

    print(f"Done with {chr} for rank {cpd_rank} with param set {param}")
    sys.stdout.flush()

for chr in chromosomes:
    for rank in ranks:
        for param in param_set:
            print(f"Starting {chr} with rank {rank}")
            cpd_fact(chr, rank, param)

            end = time.time()
            print(f"Total time: {end - start}")

print("Done!")
