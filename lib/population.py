# Population operations
import random
import numpy as np

def add_genotypes(pop_array,geno_id,times): # Used to add genotypes at the start
    empty_locs = np.where(pop_array==-1)[0]
    non_empty = np.where(pop_array>=0)[0]

    # Case I: times is smaller than (or equal to) empty_locs (prioritize replacement in empty locations)
    if times<=len(empty_locs):
        selected_indices = np.random.choice(empty_locs,times,replace=False)
        pop_array[selected_indices] = geno_id
    # Case II: times is larger than empty_locs, fill all empty_locs and fill rest randomly
    elif times>len(empty_locs):
        pop_array[empty_locs] = geno_id
        extra_indices = np.random.choice(non_empty,times-len(empty_locs),replace=False)
        pop_array[extra_indices] = geno_id

