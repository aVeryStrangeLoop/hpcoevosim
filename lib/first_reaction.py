# First reaction implementation for hpcoevosim
from lib.convert import *
import numpy as np
import random
np.seterr(divide="ignore")

def GenerateTimes(cfg,h_array,p_array,h_map,p_map):

    propensities = np.zeros((cfg.N,4))

    host_phenos = np.where(h_array==-1,-1,h_map[h_array])
    para_phenos = np.where(p_array==-1,-1,p_map[p_array])
    #print(host_phenos)
    #print(para_phenos)

    host_is_present = (host_phenos!=-1)
    para_is_present = (para_phenos!=-1)

    
    # CASE_0: No organisms
    case0_indices = np.where(np.logical_and(np.logical_not(host_is_present),np.logical_not(para_is_present)))[0]
    propensities[case0_indices] = np.array([0.0,0.0,0.0,0.0]) # [host_birth,host_death,para_birth,para_death]

    # CASE_1: Only hosts
    case1_indices = np.where(np.logical_and(host_is_present,np.logical_not(para_is_present)))[0]
    propensities[case1_indices] = np.array([cfg.r,0.0,0.0,0.0])

    # CASE_2: Only parasites
    case2_indices = np.where(np.logical_and(np.logical_not(host_is_present),para_is_present))[0]
    propensities[case2_indices] = np.array([0.0,0.0,0.0,cfg.c])

    # CASE_3: Hosts + Parasites with no phenotype match
    case3_indices = np.where(np.logical_and(np.logical_and(host_is_present,para_is_present),host_phenos!=para_phenos))[0]
    propensities[case3_indices] = np.array([cfg.r,0.0,0.0,cfg.c])

    # CASE_4: Hosts + Parasites with phenotype match
    case4_indices = np.where(np.logical_and(np.logical_and(host_is_present,para_is_present),host_phenos==para_phenos))[0]
    propensities[case4_indices] = np.array([cfg.r,cfg.d,cfg.b,cfg.c])

    times = np.random.exponential(1./propensities)
    return times

def get_newcell_id(pop_array):
    empty_indices =  np.where(pop_array==-1)[0]
    
    if len(empty_indices)!=0: # If there are empty indices, pick one of these cells
        return np.random.choice(empty_indices)
    else: # Else pick any cell
        return np.random.choice(np.where(pop_array==pop_array)[0]) 

def true_with_prob(probability):
    return random.random() < probability

def mutated(parent_gid,m,geno_len):
    as_bitstr = geno_id_to_bitstr(parent_gid,geno_len)
    as_list = list(as_bitstr)
    for idx,char in enumerate(as_list):
        is_mutated = true_with_prob(m)
        if is_mutated:
            if char=="0":
                as_list[idx] = "1"
            elif char=="1":
                as_list[idx] = "0"
    return bitstr_to_geno_id("".join(as_list))

def get_rx_times(host_gid,para_gid,hmap,pmap,cfg):
    host_is_present = (host_gid!=-1)
    para_is_present = (para_gid!=-1)
    
    # Case 0 : No organisms
    if not (host_is_present or para_is_present):
        props = np.array([0.0,0.0,0.0,0.0])
    
    # Case 1 : Only hosts
    elif (host_is_present and not(para_is_present)):
        props = np.array([cfg.r,0.0,0.0,0.0])

    # Case 2 : Only parasites
    elif (not(host_is_present) and para_is_present):
        props = np.array([0.0,0.0,0.0,cfg.c])

    # Case 3 : Both host and parasites 
    elif (host_is_present and para_is_present):
        host_pheno = hmap[host_gid]
        para_pheno = pmap[para_gid]

        if host_pheno == para_pheno:
            props = np.array([cfg.r,cfg.d,cfg.b,cfg.c])
        else:
            props = np.array([cfg.r,0.0,0.0,cfg.c])
    
    return np.random.exponential(1./props)


def ssa_step(cfg,h_array,p_array,h_map,p_map,taus):
    if h_array.min()==-1 and h_array.max()==-1 and p_array.min()==-1 and p_array.max()==-1:
        return 1.0
    min_time_idx = np.unravel_index(np.argmin(taus),taus.shape)
    min_time = np.min(taus)
    
    rxn_cell_id = min_time_idx[0]
    rxn_id = min_time_idx[1]

    taus-=min_time # Decrease all times to reactions by min

    # Now update populations and taus based on reaction ids
    if rxn_id==0: # Host birth
        parent_cellid = rxn_cell_id
        parent_gid = h_array[parent_cellid]
        parent_pheno = h_map[parent_gid]
        offspring_cellid = get_newcell_id(h_array)
        offspring_gid = mutated(parent_gid,cfg.m_h,cfg.L)
        offspring_pheno = h_map[offspring_gid]

        # Update h_array
        h_array[offspring_cellid] = offspring_gid

        # Update taus for parent and child
        taus[parent_cellid] = get_rx_times(h_array[parent_cellid],p_array[parent_cellid],h_map,p_map,cfg)
        taus[offspring_cellid] = get_rx_times(h_array[offspring_cellid],p_array[offspring_cellid],h_map,p_map,cfg)
        
    elif rxn_id==1: # Host death
        # Update h_array
        h_array[rxn_cell_id] = -1
        # Update taus for this cell
        taus[rxn_cell_id] = get_rx_times(h_array[rxn_cell_id],p_array[rxn_cell_id],h_map,p_map,cfg)

    elif rxn_id==2: # Para birth
        parent_cellid = rxn_cell_id
        parent_gid = p_array[parent_cellid]
        parent_pheno = p_map[parent_gid]
        offspring_cellid = get_newcell_id(p_array)
        offspring_gid = mutated(parent_gid,cfg.m_p,cfg.L)
        offspring_pheno = p_map[offspring_gid]

        # Update p_array
        p_array[offspring_cellid] = offspring_gid

        # Update taus for parent and child cells
        taus[parent_cellid] = get_rx_times(h_array[parent_cellid],p_array[parent_cellid],h_map,p_map,cfg)
        taus[offspring_cellid] = get_rx_times(h_array[offspring_cellid],p_array[offspring_cellid],h_map,p_map,cfg)

        
    elif rxn_id==3: # Para death
        # Update p_array
        p_array[rxn_cell_id] = -1
        taus[rxn_cell_id] = get_rx_times(h_array[rxn_cell_id],p_array[rxn_cell_id],h_map,p_map,cfg)

    return min_time
