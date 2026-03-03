import gseapy as gp
import sys

rnk_gene = sys.argv[1]
all_go_gmt = sys.argv[2]
out_dir = sys.argv[3]

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
