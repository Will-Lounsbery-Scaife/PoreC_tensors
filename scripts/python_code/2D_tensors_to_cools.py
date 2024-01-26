import os
import pickle
import numpy as np
import math
import cooler
import pandas as pd
from coolbox.api import *
import cooler

# Define the bin size used in your 2D tensors
bin_size = 100000

# Define the chromosomes you want to make cool files for
chromosomes = [ "chr20" ]
'''
chromosomes = [
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", 
    "chr21", "chr22"
]
'''

# path to directory containing your 2D tensors
pickle_files_directory = '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/tensors_2D'

# path to directory to save your cool files
cool_directory = '/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/cools'

if not os.path.exists(cool_directory):
    os.makedirs(cool_directory)


def get_chromosome_lengths():
    chromosome_lengths = {}
    chromosome_lengths["chr1"] = 248956422
    chromosome_lengths["chr2"] = 242193529
    chromosome_lengths["chr3"] = 198295559
    chromosome_lengths["chr4"] = 190214555
    chromosome_lengths["chr5"] = 181538259
    chromosome_lengths["chr6"] = 170805979
    chromosome_lengths["chr7"] = 159345973
    chromosome_lengths["chr8"] = 145138636
    chromosome_lengths["chr9"] = 138394717
    chromosome_lengths["chr10"] = 133797422
    chromosome_lengths["chr11"] = 135086622
    chromosome_lengths["chr12"] = 133275309
    chromosome_lengths["chr13"] = 114364328
    chromosome_lengths["chr14"] = 107043718
    chromosome_lengths["chr15"] = 101991189
    chromosome_lengths["chr16"] = 90338345
    chromosome_lengths["chr17"] = 83257441
    chromosome_lengths["chr18"] = 80373285
    chromosome_lengths["chr19"] = 58617616
    chromosome_lengths["chr20"] = 64444167
    chromosome_lengths["chr21"] = 46709983
    chromosome_lengths["chr22"] = 50818468
    return chromosome_lengths

chromosome_lengths = get_chromosome_lengths()

def load_tensor_from_pickle(my_file_path):
    with open(my_file_path, 'rb') as file:
        tensor = pickle.load(file)
    return tensor.numpy()  # Convert PyTorch tensor to a NumPy array

def get_number_of_bins(chromosome_length, bin_size):
    return math.ceil(chromosome_length / bin_size) + 1

def create_cool_file(interaction_matrix, chromosome, cool_output_path, bin_size):
    # Generate a minimal BED-like dataframe for bins
    n_bins = interaction_matrix.shape[0]
    bins_data = {
        'chrom': [chromosome] * n_bins,
        'start': [i * bin_size for i in range(n_bins)],
        'end': [(i + 1) * bin_size for i in range(n_bins)],
    }
    bins_df = pd.DataFrame(bins_data)

    # Generate a dataframe for the pixel values in the matrix
    pixel_data = {
        'bin1_id': [],
        'bin2_id': [],
        'count': []
    }
    for i in range(n_bins):
        for j in range(i, n_bins):
            count = interaction_matrix[i, j]
            if count > 0:  # We only record non-zero interactions
                pixel_data['bin1_id'].append(i)
                pixel_data['bin2_id'].append(j)
                pixel_data['count'].append(count)
    pixels_df = pd.DataFrame(pixel_data)

    # Create the .cool file (unbalanced)
    cooler.create_cooler(cool_uri=cool_output_path, bins=bins_df, pixels=pixels_df)
    print("Created unbalanced .cool file:", cool_output_path)


def create_balanced_cool_file(interaction_matrix, chromosome, cool_output_path, bin_size):
    # Generate a minimal BED-like dataframe for bins
    n_bins = interaction_matrix.shape[0]
    bins_data = {
        'chrom': [chromosome] * n_bins,
        'start': [i * bin_size for i in range(n_bins)],
        'end': [(i + 1) * bin_size for i in range(n_bins)],
    }
    bins_df = pd.DataFrame(bins_data)

    # Generate a dataframe for the pixel values in the matrix
    pixel_data = {
        'bin1_id': [],
        'bin2_id': [],
        'count': []
    }
    for i in range(n_bins):
        for j in range(i, n_bins):
            count = interaction_matrix[i, j]
            if count > 0:  # only record non-zero interactions
                pixel_data['bin1_id'].append(i)
                pixel_data['bin2_id'].append(j)
                pixel_data['count'].append(count)
    pixels_df = pd.DataFrame(pixel_data)

    # Create the .cool file
    cooler.create_cooler(cool_uri=cool_output_path, bins=bins_df, pixels=pixels_df)
    clr = cooler.Cooler(cool_output_path)

    # Create a balanced version of the .cool file
    balanced_cool_output_path = cool_output_path.replace('.cool', '_balanced.cool')
    
    balanced_clr = cooler.balance_cooler(
        clr, 
        cis_only=True, 
        ignore_diags=2, 
        mad_max=5, 
        min_nnz=10, 
        min_count=0, 
        store=True, 
        store_name='weight'
    )
    print("Created balanced .cool file:", balanced_cool_output_path)



for file_name in os.listdir(pickle_files_directory):
    if file_name.endswith('.pickle'):
        chromosome = file_name.split('_')[2]
        if chromosome in chromosomes:
            file_path = os.path.join(pickle_files_directory, file_name)
            interaction_tensor = load_tensor_from_pickle(file_path)
            num_bins = get_number_of_bins(chromosome_lengths[chromosome], bin_size)

            if interaction_tensor.shape == (num_bins, num_bins):
                cool_file_name = file_name.replace('.pickle', '.cool')
                balanced_cool_file_name = cool_file_name.replace('.cool', '_balanced.cool')
                cool_file_path = os.path.join(cool_directory, cool_file_name)
                balanced_cool_file_path = os.path.join(cool_directory, balanced_cool_file_name)
                create_balanced_cool_file(interaction_tensor, chromosome, balanced_cool_file_path, bin_size)
                create_cool_file(interaction_tensor, chromosome, cool_file_path, bin_size)
                print("success")

print("Script finished successfully!")