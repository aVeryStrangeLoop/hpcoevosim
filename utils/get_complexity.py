# Takes in a pheno_pops.csv file and returns the temporal host and parasite heterogeneity
import sys
import math
import numpy as np


# Main function
def get_complexity(filename,hmapfile,pmapfile):
    num_pheno = -1
    with open(filename) as ifile:
        header = next(ifile)
        num_pheno = (len(header.split(","))-1)/2

    host_pheno_nums = np.zeros(int(num_pheno))
    para_pheno_nums = np.zeros(int(num_pheno))    
    
    with open(hmapfile) as hfile:
        next(hfile)
        for line in hfile:
            words = line.split(",")
            if len(words)>1:
                host_pheno_nums[int(words[1])] += 1
          
    with open(pmapfile) as pfile:
        next(pfile)
        for line in pfile:
            words = line.split(",")
            if len(words)>1:
                para_pheno_nums[int(words[1])] += 1

    host_pheno_comps = math.log(np.sum(host_pheno_nums),2.) - np.log2(host_pheno_nums)
    para_pheno_comps = math.log(np.sum(para_pheno_nums),2.) - np.log2(para_pheno_nums)

    times = []
    host_comps = []
    para_comps = []
    with open(filename) as ifile:
        next(ifile)
        for line in ifile:
            words = line.split(",")
            if len(words)>1:
                time = float(words[0])
                host_pops = [int(i) for i in words[1:int(1+num_pheno)]]
                para_pops = [int(i) for i in words[1+int(num_pheno):]]
                host_comp = np.sum(np.multiply(host_pheno_comps,host_pops))/np.sum(host_pops) 
                para_comp = np.sum(np.multiply(para_pheno_comps,para_pops))/np.sum(para_pops)
                times.append(time)
                host_comps.append(host_comp)
                para_comps.append(para_comp)

    host_c0 = np.sum(np.multiply(host_pheno_comps,host_pheno_nums))/np.sum(host_pheno_nums)
    para_c0 = np.sum(np.multiply(para_pheno_comps,para_pheno_nums))/np.sum(para_pheno_nums)

    return times,host_comps,para_comps,host_c0,para_c0

def save(times,host_comps,para_comps,host_c0,para_c0):
    with open("complexity.csv","w+") as ofile:
        ofile.write("time,host_comp,para_comp,host_esc,para_esc\n")
        for idx,time in enumerate(times):
            ofile.write("%f,%f,%f,%f,%f\n" % (time,host_comps[idx],para_comps[idx],host_comps[idx]/host_c0,para_comps[idx]/para_c0))


if __name__=="__main__":
        filename = sys.argv[1]
        hmapfile = sys.argv[2]
        pmapfile = sys.argv[3]
        times,host_comps,para_comps,host_c0,para_c0 = get_complexity(filename,hmapfile,pmapfile)
        save(times,host_comps,para_comps,host_c0,para_c0)
