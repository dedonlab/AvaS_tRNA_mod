#!/usr/bin/env python
import pandas as pd
import os, sys
from collections import defaultdict

# input file = xxx/ffn_prep/xxx.PATRIC.ffn.
tsv_ffn = sys.argv[1]
dir_w = sys.argv[2]

# create output folder.
dir_out = dir_w + "/aua_n_count"
os.makedirs(dir_out, exist_ok=True)

# parse output files.
bname = os.path.basename(tsv_ffn)
tsv_output = dir_out + "/"+ bname + "_aua_n_count.txt"

lst_nt = ["a", "c", "g", "t"]

# Create an empty DataFrame with column names.
dict_codon = defaultdict(list)

with open(tsv_ffn) as ffn:
    Lines = ffn.readlines()     # the first empty line had been removed, so no need add [1:].
    for line in Lines :
        if line.startswith('>') :
            dict_codon['CDS'].append(line[:-1])
        else :
            for x in lst_nt :
                codon = 'ata ' + x
                num_codon = line.count(codon)
                dict_codon[codon].append(num_codon)
            for x in lst_nt :
                codon = ' ' + x
                num_codon = line.count(codon)
                dict_codon[codon].append(num_codon)

df=pd.DataFrame(dict_codon)

df.to_csv(tsv_output, sep = '\t', mode='w+',index=False, header=not os.path.isfile(tsv_output))

# END.
