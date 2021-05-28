from lib.config import *
from lib.first_reaction import *
import shutil
import os
import numpy as np

def get_savestr(h_array,p_array,h_map,p_map,savetime,cfg):
    host_phenos = np.where(h_array==-1,-1,h_map[h_array])
    para_phenos = np.where(p_array==-1,-1,p_map[p_array])
    writestr = "%f," % savetime
    for host_pid in range(cfg.nP):
        this_pop = np.count_nonzero(host_phenos==host_pid)
        writestr += "%d," % this_pop
    for para_pid in range(cfg.nP):
        this_pop = np.count_nonzero(para_phenos==para_pid)
        writestr += "%d," % this_pop
    writestr = writestr[:-1]+"\n"

    return writestr
    
    

def main():

    # Get Configuration
    print("Initialising configuration and initial populations...")
    config = Config("./config/params.cfg")
    h_array,p_array = GetInit(config,"./config/host_initpop.csv","./config/para_initpop.csv") # Get initial populations
    h_map,p_map = GetMaps(config,"./config/host_gptable.csv","./config/para_gptable.csv") # Get gp-maps

    taus = GenerateTimes(config,h_array,p_array,h_map,p_map) # Generate the times to different reactions in each world cell
    
    # Generate output folder and file
    if os.path.exists('output'):
        shutil.rmtree('output')
        os.makedirs('output')
    else:
        os.makedirs('output')

    outfile = open("./output/pheno_pops.csv","w+")
    outfile_header = "time,"
    for pidx in range(config.nP):
        outfile_header += "host_p%d," % pidx
    for pidx in range(config.nP):
        outfile_header += "para_p%d," % pidx
    outfile_header = outfile_header[:-1]+"\n"
    outfile.write(outfile_header)

    T_GLOBAL = 0.0
    outfile.write(get_savestr(h_array,p_array,h_map,p_map,0.0,config))

    save_approach = 0.0
    save_counter = 0
    while T_GLOBAL < config.runtime:
        delta_t = ssa_step(config,h_array,p_array,h_map,p_map,taus)
        T_GLOBAL += delta_t
        save_approach += delta_t
        if save_approach > config.save_every:
            save_counter += 1
            save_approach -= config.save_every
            savetime = save_counter * config.save_every
            outfile.write(get_savestr(h_array,p_array,h_map,p_map,savetime,config))
            print("t=%f out of %f" % (savetime,config.runtime))

    outfile.write(get_savestr(h_array,p_array,h_map,p_map,config.runtime,config))
        

        
    



if __name__=="__main__":
    main()
