# Takes in a pheno_pops.csv file and returns the temporal host and parasite heterogeneity
import sys
import math
import numpy as np

# Main function
def get_dominant(filename):
    times = []
    host_doms = []
    para_doms = []
    with open(filename) as ifile:
        header = next(ifile)
        num_pheno = (len(header.split(","))-1)/2
        for line in ifile:
            words = line.split(",")
            if len(words)>1:
                time = float(words[0])
                host_pops = [int(i) for i in words[1:int(1+num_pheno)]]
                para_pops = [int(i) for i in words[1+int(num_pheno):]]
                host_dom = np.argmax(host_pops)
                para_dom = np.argmax(para_pops)
                times.append(time)
                host_doms.append(host_dom)
                para_doms.append(para_dom)
    return times,host_doms,para_doms    

def save_dominant(times,host_doms,para_doms):
    with open("output/dominant.csv","w+") as ofile:
        ofile.write("time,host_dom,para_dom\n")
        for idx,time in enumerate(times):
            ofile.write("%f,%d,%d\n" % (time,host_doms[idx],para_doms[idx]))


