#!/usr/bin/env python
import pandas as pd
import os, sys
from collections import defaultdict

# input file = xxx/ffn_prep/xxx.PATRIC.ffn.
tsv_ffn = sys.argv[1]
dir_w = sys.argv[2]
#tsv_ffn=
#dir_w=/blue/lagard/yuanyifeng/PA14_24190

# create output folder.
dir_out = dir_w + "/py_count"
os.makedirs(dir_out, exist_ok=True)

# parse output files.
bname = os.path.basename(tsv_ffn)
tsv_output = dir_out + "/"+ bname + "_count.txt"

# Create an empty DataFrame with column names.
dict_codon = defaultdict(list)

lst_nt = ["a", "c", "g", "t"]

with open(tsv_ffn) as ffn:
    Lines = ffn.readlines()     # the first empty line had been removed, so no need: Lines = ffn.readlines()[1:]
    for line in Lines :
        if line.startswith('>') :
            dict_codon['CDS'].append(line[:-1])
        else :
            lst_codons = line.split(" ")
            for x in lst_nt :
                for y in lst_nt :
                    for z in lst_nt :
                        codon = x + y + z
                        num_codon = lst_codons[1:].count(codon)  #add [1:] to remove the start codon.
                        dict_codon[codon].append(num_codon)

df=pd.DataFrame(dict_codon)

df.to_csv(tsv_output, sep = '\t', mode='w+',index=False, header=not os.path.isfile(tsv_output))

# END.
