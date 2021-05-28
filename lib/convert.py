def bitstr_to_geno_id(bitstr):
    return int(bitstr,2)

def geno_id_to_bitstr(id,geno_len):
    as_str = "{0:b}".format(int(id))
    as_str = ("0"*(geno_len-len(as_str)))+as_str
    return as_str