import torch
import pandas as pd
import numpy as np
import itertools
import pickle
import os
import sys
import math
import gc

# This script creates 2D tensors from Pore-C, one tensor per chromosome. Only reads with cardinality = 2 are stored.

####################
### PARAMETERS #####
####################

# Bin size, in bp
bin_size = 100000

# The minimum frequency for a binset to be included in the tensor
# Higher min_freq --> higher quality reads, but more sparse tensor
min_freq = 2

# Which chromosomes to create tensors for
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

output_directory = "/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/output/"
tensor_directory = os.path.join(output_directory, 'tensors_2D')

####################
####################
####################


# if tensor_directory does not exist, create it
if not os.path.exists(tensor_directory):
    os.makedirs(tensor_directory)

print("starting 2D tensors!")
sys.stdout.flush()

max_cardinality = 2

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
    chromosome_lengths["chrX"] = 156040895
    chromosome_lengths["chrY"] = 57227415
    return chromosome_lengths

chromosome_lengths = get_chromosome_lengths()

def create_tensors(cardinality_data, bin_size, chrom):
    chr_length = chromosome_lengths[chrom]
    dim = math.ceil(chr_length / bin_size) + 1

    # Initialize a PyTorch tensor with zeros
    interaction_tensor = torch.zeros((dim, dim), dtype=torch.int32)

    grouped = cardinality_data.groupby('read_ID')
    aggregated = grouped.agg({'bin': list, 'midpoint': list})

    final_card = aggregated.reset_index()
    final_card = final_card.rename(columns={'bin': 'bins'})
    final_card = final_card.rename(columns={'midpoint': 'midpoints'})

    # remove duplicate bins
    final_card['bins'] = final_card['bins'].apply(lambda x: list(set(x)))

    # sort bins, ascending
    final_card['bins'] = final_card['bins'].apply(lambda x: sorted(x))

    # Delete intermediate DFs to free up memory
    del grouped
    del aggregated
    gc.collect()

    final_card = final_card[final_card['bins'].map(len) == max_cardinality]
    print(final_card.shape)
    count_dict = {}

    tik = 0
    increment_count = 0
    for _, contact in final_card.iterrows():
        tik += 1
        fragment_bins = contact.iloc[1]

        frags_with_zeros = [0] * max_cardinality
        frags_with_zeros[:len(fragment_bins)] = fragment_bins
        frags_with_zeros[len(fragment_bins):] = [0] * (max_cardinality - len(fragment_bins))

        fragment_bins_permutations = list(itertools.permutations(frags_with_zeros))

        for perm in fragment_bins_permutations:
            fragment_bins = list(perm)
            coords = [0] * max_cardinality
            coords[:len(fragment_bins)] = fragment_bins
            coords[len(fragment_bins):] = [0] * (max_cardinality - len(fragment_bins))    
            coords_tuple = tuple(coords)

            increment_count += 1

            # check if coords_tuple exists as a key in count_dict
            if coords_tuple in count_dict:
                count_dict[coords_tuple] += 1
            else:
                count_dict[coords_tuple] = 1

    # iterate over count_dict and remove any entries with a count less than 2
    for key, value in list(count_dict.items()):
        if value < 2:
            del count_dict[key]

    tensor_shape = torch.Size([dim]*max_cardinality)
    indices = list(count_dict.keys())
    values = list(count_dict.values())

    tok = 0
    # Fill the tensor using count_dict
    for coords, count in count_dict.items():
        tok +=1
        
        i, j = coords[:2]
        interaction_tensor[i, j] += count

        if tok % 10000 == 0:
            print("\ncoords:", coords)
            print("frequency so far:", count)
            print("i, j:", i, j)
            print("interaction_tensor[i, j]:", interaction_tensor[i, j])
            print("interaction_tensor[j, i]:", interaction_tensor[j, i])
            print("count_dict at i, j:", count_dict[(i, j)])
            print("count_dict at j, i:", count_dict[(j, i)])
            sys.stdout.flush()
        
    print("\n\ncreated tensor for " + chrom)
    print("2 dimensions: ", dim)
    print("total increment count: ", increment_count)
    sys.stdout.flush()
    
    return interaction_tensor

def save_tensor_to_file(tensor_list, filename):
    with open(filename, 'wb') as f:
        pickle.dump(tensor_list, f)
    print("saved", filename)
    sys.stdout.flush()
    
def create_chromosome_tensors(cardinality_data, bin_size, df_name, chromo):

    print("Creating tensors for " + chromo)
    sys.stdout.flush()

    chrom_tensor_var = create_tensors(cardinality_data, bin_size, chromo)
    
    filename = os.path.join(tensor_directory, f"{df_name}_tensor.pickle")
    save_tensor_to_file(chrom_tensor_var, filename)

    # delete chrom_tensor_var to free up memory
    del chrom_tensor_var
    gc.collect()


# path to directory containing gz files
byChr_path = "/commons/groups/gursoy_lab/wlounsberyscaife/PoreC_Tensors/data/PoreC_reads_by_Chromosome"

# get a list of all .gz files in the directory
gz_files = [f for f in os.listdir(byChr_path) if f.endswith('.gz')]

# loop over each gz file
for gz_file in gz_files:
    
    # Split the filename to determine which chromosome the file corresponds to
    current_chromosome = gz_file.split("_")[-1].split(".")[0]

    if current_chromosome not in chromosomes:
        continue

    print(f"\nReading in {gz_file}: {current_chromosome}")
    sys.stdout.flush()

    # read gz_file into a dataframe
    df = pd.read_csv(os.path.join(byChr_path, gz_file), sep='\t', header=None, usecols=[0, 1, 2, 3, 4])
    
    # rename columns: 0=chrom, 1=frag_start, 2=frag_end, 3=read_ID, 4=read_len
    df = df.rename(columns={0: "chrom", 1: "frag_start", 2: "frag_end", 3: "read_ID", 4: "read_len"})

    # Add a bin column and a midpoint column for each fragment
    df['bin'] = 0
    df['midpoint'] = (df['frag_start'] + df['frag_end']).floordiv(2)
    df.loc[:, 'bin'] = (df['midpoint'] / bin_size).apply(np.floor).astype(int) + 1

    # remove the ".gz" extension from the filename to use as the dictionary key
    key = gz_file[:-3]
    
    create_chromosome_tensors(df, bin_size, key, current_chromosome)

    del df
    gc.collect()

print("script has finished successfully!")
