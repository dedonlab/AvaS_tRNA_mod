## Description
This repository contains commands and scripts that haven been generated and used by the Prof. Peter Dedon’s group at MIT and the Singapore-MIT Alliance for Research and Technology (SMART)  for codon usage analysis in the manuscript:  
**“Pyridoxal phosphate-dependent biosynthesis of aminovaleramide by AvaS in tRNA”**  
Authors: *Jingjing Sun, Junzhou Wu, Yifeng Yuan, Seetharamsing Balamkundu, Agnieszka Dziergowska, Hazel Chay Suen Suen, Dwijapriya, Liang Cui, Léo Hardy, Megan En Lee, Grazyna Leszczynska, Chuan-Fa Liu, Zeynep Baharoglu, Laurence Drouard, Steven D. Bruner, Thomas J. Begley, Valérie de Crécy-Lagard, Peter C. Dedon*

## Help and Issues
Contact Yifeng Yuan at yuanyifeng@ufl.edu

## Contributors
Yifeng Yuan, Ph.D. and Peter C. Dedon (Principal Investigator)

Version History  
v0.9.0 --2026/03/02 --create the repository and migrate data.

## Dependencies
mafft v7.520 https://mafft.cbrc.jp/alignment/software/  
MASH v2.3 https://github.com/marbl/mash  
seqkit v2.8.0 https://bioinf.shenwei.me/seqkit/  
python v3.12 https://www.python.org/

## Usage
1. Dowload genome IDs in a list. For example, Select representative genomes using filters: Organism: bacteria, Reference: representative, Quality: Good, Status: Complete.
bac_rep_4245.txt (4245 representative bacterial genomes)
Repeat the step for other genus and species.
Aeromonas_gid.txt (Aeromonas)
Flavobacterium_gid.txt (Flavobacterium)
Escherichia_gid.txt (Escherichia)
Campylobacter_gid.txt (Campylobacter)
Acinetobacter_gid.txt (Acinetobacter)
other_gid.txt (other genus including Chitinophaga, Hymenobacter, Pedobacter, Pseudoalteromonas)
vibrio_gid.txt (Vibrio)
Streptomyces_gid.txt (Streptomyces)
Shewanella_gid.txt (Shewanella)
pseudo_gid.txt (Pseudomonas)
Arcobacter_gid.txt (Arcobacter)
Chryseobacterium_gid.txt (Chryseobacterium)
Abmn_gid.txt (Acinetobacter baumannii)
Ahyd_gid.txt (Aeromonas hydrophila)
Sone_gid.txt (Shewanella oneidensis)
Vcho_gid.txt (Vibrio cholerae)
Vpara_gid.txt (Vibrio parahaemolyticus)

3. Retrieve data from BV-BRC (https://www.bv-brc.org/), including CDS sequences (ffn files), genome sequences (fna files) and protein sequences (faa files)
4. 




```
mkdir fnafiles || true
mkdir ffnfiles || true
mkdir faafiles || true

for id in $(cat wget_gid.txt);  do
  wget -qN "ftp://ftp.patricbrc.org/genomes/${id}/${id}.fna" -P fnafiles
  wget -qN "ftp://ftp.patricbrc.org/genomes/${id}/${id}.PATRIC.ffn" -P ffnfiles
  wget -qN "ftp://ftp.patricbrc.org/genomes/${id}/${id}.PATRIC.faa" -P faafiles
done
```

#### 2. phylogeny tree
#### 2.1 MSA  
Retrieve selected protein 

mafft --thread 12 --maxiterate 1000 --localpair "${aln}".fasta > "${aln}"_mafft.aln

# trim alignment
module load bmge/1.12
java -jar /apps/bmge/1.12/bin/BMGE.jar -m BLOSUM30 -i "${aln}"_mafft.aln -t AA -of "${aln}"_BMGE.fasta

# convert fasta format to phy
python3 /blue/lagard/yuanyifeng/scripts/fasta2phy.py -i "${aln}"_BMGE.fasta -o "${aln}"_BMGE.phy

# tree building
module load raxml-ng/1.1.0
raxml-ng --msa "${aln}"_BMGE.fasta --model LG+G+F --prefix result --seed 123 \
         --search replicates=100 --threads 12 --tree pars{100},rand{100}

# bootstrap
raxml-ng --bootstrap --msa "${aln}"_BMGE.fasta --bs-trees 1000 --seed 123 \
        --prefix boot --model LG+G+F

# bootstrap comparison
raxml-ng --support --tree result.raxml.bestTree --threads 12 \
         --bs-trees boot.raxml.bootstraps --prefix support


# Blastp
# performed on Hipergator at hpg.rc.ufl.ed

module load ncbi_blast/2.15.0
# build blastDB
dir_db=/blue/lagard/yuanyifeng/BV-DB/faafiles

ls ${dir_db}/*faa | xargs -i makeblastdb -in {} -dbtype prot -parse_seqids -out {}_blastdb

dir_o=outfile_cog
qfaa= # name query sequences here.

mkdir ${dir_o} || true

for db in $(ls ${dir_db}/*faa); do
  dbfaaname=$(basename $db)
  dbname=${dbfaaname%%.PATRIC.faa}
  out=${dir_o}/${dbname}.out
  blastp -query "$qfaa" -db "$db"_blastdb -out "$out" \
         -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" \
         -num_threads 3 -max_target_seqs 10000
done

# tBlastn
# prepare a list of genome id for genomes of each genus selected. Name them "${org}"_gid.txt

dir_o=outfiles_tblastn_"${org}"
mkdir ${dir_o} || true

qfaa= # name query sequences here.
org= # name the genus for searching here.
dir_db=/blue/lagard/yuanyifeng/PA14_24190/bv_"${org}"/fnafiles

## build blastDB
ls ${dir_db}/*fna | xargs -i makeblastdb -in {} -dbtype nucl -out {}_blastdb

for gid in $(cat "${org}"_gid.txt); do
  genome="${dir_db}"/"${gid}".fna
  makeblastdb -in "${genome}" -dbtype nucl -out "${genome}"_blastdb
  out=${dir_o}/${gid}.aln
  tblastn -query "$qfaa" -db "${genome}"_blastdb -out "$out" \
         -num_threads 6 -max_target_seqs 10 \
         -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles"
         #-outfmt 0
done

# MASH
module load mash/2.3
mash sketch -p 6 -o abs_sketch -l avaS_abs_gid.txt
mash sketch -p 6 -o pre_sketch -l avaS_pre_gid.txt
mash dist -p 6 abs_sketch.msh pre_sketch.msh > mash_result.txt

# extract mash result
for g in $(cat gid.txt | sed 's/^g//') ; do
  grep "${g}.fna" mash_result.txt | sed -E 's#/blue/lagard/yuanyifeng/PA14_24190/bv_[A-Za-z]+/fnafiles/#g#g' | sort -k3,3n | sed 's/\.fna//g' | head -1
done

# Codon analysis
# Run the attached scripts below
1_pycodon_count.sh
2_count2CDScodon.sh
3_AUA_Nstats.sh
3_AUAstats.sh
4_CDS_AUA_N_summary.sh

# GESA

# 1. prepare GO annotation for each organism in EXCEL and save as xxx_go_bvbrc.txt
# convert to .gmt format

python convert_to_gmt.py Abmn_go_bvbrc.txt Abmn_go.gmt
python convert_to_gmt.py Aero_go_bvbrc.txt Aero_go.gmt
python convert_to_gmt.py PAO1_go_bvbrc.txt PAO1_go.gmt
python convert_to_gmt.py Sone_go_bvbrc.txt Sone_go.gmt
python convert_to_gmt.py Vcho_go_bvbrc.txt Vcho_go.gmt

# 2. prepare Gene list sorted with aua synonmus ratio.
for file in $(ls *.PATRIC.ffn_CDScodon.tsv); do
  output=${file%%.PATRIC.ffn_CDScodon.tsv}
  awk -F '\t' 'FNR>1{if ($14+$15+$17>0) {print $1"\t"$14/($14+$15+$17)} else {print $1"\t"0}}' "${file}" | sort -k2,2nr > "${output}"_CDS_aua_syn_sort.txt
done

for file in $(ls *_CDS_aua_syn_sort.txt); do
  # remove some |foo bar following fig ID
  sed -i 's/^>fig|/#/' "${file}"
  sed -i 's/|.*\t/\t/' "${file}"
  sed -i 's/^#/fig|/' "${file}"
  # replace | and . with _ in fig ID.
  sed -i 's/|/_/' "${file}"
  sed -i 's/\./_/' "${file}"
  sed -i 's/\.peg\./_peg_/' "${file}"
done

# 3. prepare Gene list sorted with aua number
for file in $(ls *.PATRIC.ffn_CDScodon.tsv); do
  output=${file%%.PATRIC.ffn_CDScodon.tsv}
  awk -F '\t' 'FNR>1{print $1"\t"$14}' "${file}" | sort -k2,2nr > "${output}"_aua_num_sort.txt
done

for file in $(ls *_aua_num_sort.txt); do
  # remove some |foo bar following fig ID
  sed -i 's/^>fig|/#/' "${file}"
  sed -i 's/|.*\t/\t/' "${file}"
  sed -i 's/^#/fig|/' "${file}"
  # replace | and . with _ in fig ID.
  sed -i 's/|/_/' "${file}"
  sed -i 's/\./_/' "${file}"
  sed -i 's/\.peg\./_peg_/' "${file}"
done

#set qt path to the result above
export QT_QPA_PLATFORM_PLUGIN_PATH=/usr/lib/x86_64-linux-gnu/qt5/plugins/platforms

# 3 Run gseapy preranked analysis
# run the python3 script below

import gseapy as gp

# ranked gene list. column 1 gene id , column 2 value such as Foldchange, aua number
# rnk_gene = '470.11567_CDS_aua_syn_sort.txt'

rnk_gene = '666.4624_aua_num_sort.txt'
all_go_gmt = 'Vcho_go.gmt' # all go annotation in gmt format
out_dir = 'aua_num_Vcho'

pre_res = gp.prerank(
    rnk=rnk_gene,
    gene_sets=all_go_gmt,
    outdir= out_dir,
    permutation_num=20000,
    min_size=10,
    max_size=500,
    seed=123
)

# View top results
#print(pre_res.res2d.head())

# Extract the results table
df = pre_res.res2d.reset_index()
# df.head()

# plot as the traditional GO plot

# Plot top 5 or 10 enriched GO terms using seaborn
# Filter and sort top 5 by lowest FDR q-value
df = df.sort_values("FDR q-val").head(5)

# Add GO descriptions from the gmt file
desc_map = {}
with open(all_go_gmt) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            desc_map[parts[0]] = parts[1]  # GO ID : description

# Add a new column to dataframe
df["Description"] = df["Term"].map(desc_map)
df["Label"] = df["Term"] + " — " + df["Description"]  # Optional: both ID and desc
df["logFDR"] = df["FDR q-val"].apply(lambda x: -np.log10(float(x) + 1e-10))

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import math

def round_down_to_one_decimal(x):
    return math.floor(x * 10) / 10

def generate_steps(x, n=3):
    # Find the max value ≤ x that has only 1 decimal place
    x_rounded = round(np.floor(x * 10) / 10, 1)
    # Create possible values with 1 decimal place
    candidates = [round(i / 10, 1) for i in range(1, int(x_rounded * 10) + 1)]
    # Select evenly spaced values (3 values)
    if len(candidates) < n:
        raise ValueError("Not enough values to generate steps")
    idxs = np.linspace(0, len(candidates) - 1, n, dtype=int)
    return [candidates[i] for i in idxs]


# check the range of q-value
print(df["FDR q-val"])

# set 0.05 as the middle color of hue (grey) using TwoSlopeNorm.
# add ticks on the colorbar that smaller than 0.05 (on the red side).
# they should be in the range of q-val

from matplotlib.colors import TwoSlopeNorm

vmin = df["FDR q-val"].min()
vmax = df["FDR q-val"].max()
vcenter = 0.05

if vmin <= 0.05 and vmax >= 0.5 :
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    ticks = [0.01, 0.03, 0.05, 0.25, 0.5, 0.75, 0.1]
elif vmin <= 0.05 and vmax > 0.1 :
    norm = TwoSlopeNorm(vmin=0.01, vcenter=vcenter, vmax=round_down_to_one_decimal(vmax))
    ticks = [0.01, 0.03, 0.05] + generate_steps(vmax, n=3)[1:]
elif vmin <= 0.05 and vmax < 0.1 :
    norm = TwoSlopeNorm(vmin=0.01, vcenter=vcenter, vmax=0.1)
    ticks = [0.01, 0.03, 0.05, 0.1]
else :
    norm = TwoSlopeNorm(vmin=0.04, vcenter=vcenter, vmax=1)
    ticks = [0.04, 0.05, 0.25, 0.5, 0.75, 1]

# mapping q-val to color
import matplotlib.cm as cm

cmap = cm.get_cmap("coolwarm_r")
colors = [cmap(norm(val)) for val in df["FDR q-val"]]

# Plot GO description
fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(
    df["NES"],
    df["Term"],
    s=df["logFDR"] * 500+5,  # the larger the more significant (small q-val)
    c=colors)

# add colorbar and legend
sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.5) # shrink=0.5, make it shorter
cbar.set_label("FDR q-value")
cbar.set_ticks(ticks)
#cbar.ax.invert_yaxis()  # option

# reverse the y axis
ax.invert_yaxis()

ax.set_xlabel("Normalized Enrichment Score (NES)")
ax.set_ylabel("GO Term")
ax.set_title("Top Enriched GO Terms")

plt.tight_layout()
plt.show()

