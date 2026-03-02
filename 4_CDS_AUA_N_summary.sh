#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=10-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
###############################################################

output=AUA_N_stats_05202025.txt

echo -e "gid\ttotal_nnn_a\ttotal_nnn_c\ttotal_nnn_g\ttotal_nnn_t\taua_a\taua_c\taua_g\taua_t\taua_a/aua_h" > "${output}"

for file in $(ls /blue/lagard/yuanyifeng/PA14_24190/codon/aua_n_count/*.PATRIC.ffn_aua_n_count.txt); do
  gid=g"$(basename "$file" | sed 's/.PATRIC.ffn_aua_n_count.txt//')"
  total_aua_a=$(awk -F '\t' '{print $2}' "${file}" | paste -sd+ | bc)
  total_aua_c=$(awk -F '\t' '{print $3}' "${file}" | paste -sd+ | bc)
  total_aua_g=$(awk -F '\t' '{print $4}' "${file}" | paste -sd+ | bc)
  total_aua_t=$(awk -F '\t' '{print $5}' "${file}" | paste -sd+ | bc)
  total_nnn_a=$(awk -F '\t' '{print $6}' "${file}" | paste -sd+ | bc)
  total_nnn_c=$(awk -F '\t' '{print $7}' "${file}" | paste -sd+ | bc)
  total_nnn_g=$(awk -F '\t' '{print $8}' "${file}" | paste -sd+ | bc)
  total_nnn_t=$(awk -F '\t' '{print $9}' "${file}" | paste -sd+ | bc)
  
  total_aua_h=$(echo "${total_aua_c} + ${total_aua_g} + ${total_aua_t}" | bc)
  ratio_aua_a_over_aua_h=$(echo "scale=6; ${total_aua_a} / ${total_aua_h} " | bc)

  echo -e "${gid}\t${total_aua_a}\t${total_aua_c}\t${total_aua_g}\t${total_aua_t}\t${total_nnn_a}\t${total_nnn_c}\t${total_nnn_g}\t${total_nnn_t}\t${ratio_aua_a_over_aua_h}" >> "${output}"
done

