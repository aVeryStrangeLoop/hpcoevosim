# Takes in a pheno_pops.csv file and returns the temporal host and parasite heterogeneity
import sys
import math

def entropy(arr):
    entr = 0.0
    tot = sum(arr)
    for i in arr:
        frac = float(i)/float(tot)
        if frac!=0.0:
            entr -= frac*math.log(frac,2.0)
    return entr

# Main function
def get_diversity(filename):
    times = []
    host_hets = []
    para_hets = []
    with open(filename) as ifile:
        header = next(ifile)
        num_pheno = (len(header.split(","))-1)/2
        for line in ifile:
            words = line.split(",")
            if len(words)>1:
                time = float(words[0])
                host_pops = [int(i) for i in words[1:int(1+num_pheno)]]
                para_pops = [int(i) for i in words[1+int(num_pheno):]]
                host_het = entropy(host_pops)
                para_het = entropy(para_pops)
                times.append(time)
                host_hets.append(host_het)
                para_hets.append(para_het)
    return times,host_hets,para_hets     

def save_diversity(times,host_hets,para_hets):
    with open("output/diversity.csv","w+") as ofile:
        ofile.write("time,host_div,para_div\n")
        for idx,time in enumerate(times):
            ofile.write("%f,%f,%f\n" % (time,host_hets[idx],para_hets[idx]))


