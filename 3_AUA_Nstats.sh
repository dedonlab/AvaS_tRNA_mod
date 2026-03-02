#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanyifeng@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#############################################

dir_w=$1

for file in $(ls ${dir_w}/ffn_prep/*\.ffn) ; do
  python3 3_count_CDS_AUA_N.py "${file}" "${dir_w}"
done
