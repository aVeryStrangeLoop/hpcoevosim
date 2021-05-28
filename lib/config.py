import numpy as np
from lib.convert import *
from lib.population import *

class Config:
    def __init__(self,param_filepath):
        self.r = 0.0 # Host-birth rate (H -> 2H)
        self.b = 0.0 # Parasite reproduction (H + P -> H + 2P)
        self.d = 0.0 # Host death rate with parasite (H + P -> P)
        self.c = 0.0 # Parasite death rate (P -> 0)
        self.m_h = 0.00 # Host mutation probability
        self.m_p = 0.00 # Parasite mutation probability
        self.L = 0 # Length of genome
        self.nP = 0 # Number of phenotypes
        self.N = 0 # World size
        self.runtime = 0.0 # Runtime
        self.save_every = 0.0 # Save_every

        with open(param_filepath) as pfile:
            for line in pfile:
                if line[0]=="#" or line[0]=="\n": # If line is empty or has comment skip
                    continue 
                row = line.split(" ")
                if row[0]=="r":
                    self.r = float(row[1])
                elif row[0]=="b":
                    self.b = float(row[1])
                elif row[0]=="d":
                    self.d = float(row[1])
                elif row[0]=="c":
                    self.c = float(row[1])
                elif row[0]=="m_host":
                    self.m_h = float(row[1])
                elif row[0]=="m_para":
                    self.m_p = float(row[1])
                elif row[0]=="L":
                    self.L = int(row[1])
                elif row[0]=="nP":
                    self.nP = int(row[1])
                elif row[0]=="N":
                    self.N = int(row[1])
                elif row[0]=="runtime":
                    self.runtime = float(row[1])
                elif row[0]=="save_every":
                    self.save_every = float(row[1])
                else:
                    print("Invalid parameter declaration : %s = %s" % (row[0],row[1]))


def GetInit(cfg,hinitfile,pinitfile):
    h_array = np.full(cfg.N,-1)
    p_array = np.full(cfg.N,-1)

    with open(hinitfile) as hifile:
        next(hifile)
        for line in hifile:
            words = line.split(",")
            if len(words)==2:
                if int(words[1])!=0:
                    add_genotypes(h_array,bitstr_to_geno_id(words[0]),int(words[1])) # Add genotypes to h_array

    with open(pinitfile) as pifile:
        next(pifile)
        for line in pifile:
            words = line.split(",")
            if len(words)==2:
                if int(words[1])!=0:
                    add_genotypes(p_array,bitstr_to_geno_id(words[0]),int(words[1])) # Add genotypes to h_array

    return h_array,p_array

def GetMaps(cfg,hmapfile,pmapfile):
    h_map = np.full(2**cfg.L,-1)
    p_map = np.full(2**cfg.L,-1)

    with open(hmapfile) as hmfile:
        next(hmfile)
        for line in hmfile:
            words = line.split(",")
            if len(words)==2:
                geno_id = bitstr_to_geno_id(words[0])
                h_map[geno_id] = int(words[1])

    with open(pmapfile) as pmfile:
        next(pmfile)
        for line in pmfile:
            words = line.split(",")
            if len(words)==2:
                geno_id = bitstr_to_geno_id(words[0])
                p_map[geno_id] = int(words[1])

    return h_map,p_map

