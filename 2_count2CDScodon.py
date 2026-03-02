#!/usr/bin/env python
import pandas as pd
import os, sys

# create output folder.
dir_out = "CDS_codon"
os.makedirs(dir_out, exist_ok=True)

tsv_count = sys.argv[1]

# parse output files.
le = len("_count.txt")
bname = os.path.basename(tsv_count)[:-le]
tsv_output = dir_out + "/"+ bname + "_CDScodon.tsv"

lst_nt = ["A", "C", "G", "U"]

lst_61codon = []
for a in lst_nt :
    for b in lst_nt :
        for c in lst_nt :
            codon = a + b + c
            lst_61codon.append(codon)

df = pd.read_table(tsv_count, sep = '\t', header=0, names=["CDS"] + lst_61codon)

df['NAU'] = df['AAU'] + df['CAU'] + df['GAU'] + df['UAU']
df['NAC'] = df['AAC'] + df['CAC'] + df['GAC'] + df['UAC']
df['CDStotalcodon'] = df[lst_61codon].sum(axis=1)
df['cdsNAUnorm'] = df['NAU'].div(df.CDStotalcodon, axis=0)
df['cdsNACnorm'] = df['NAC'].div(df.CDStotalcodon, axis=0)

gTotalNAU = df.loc[:,'NAU'].sum()
gTotalNAC = df.loc[:,'NAC'].sum()
gTotal = df['CDStotalcodon'].sum()

df['gTotalNAU'] = gTotalNAU
df['gTotalNAC'] = gTotalNAC
df['genometotalcodon'] = gTotal
df['gNAUnorm'] = df['gTotalNAU'].div(df.genometotalcodon, axis=0)
df['gNACnorm'] = df['gTotalNAC'].div(df.genometotalcodon, axis=0)

df['cdsNAU/C'] = df['NAU'].div(df.NAC, axis=0)
df['gNAU/C'] = gTotalNAU / gTotalNAC

df.to_csv(tsv_output, sep = '\t', mode='w',index=False, header=not os.path.isfile(tsv_output))

