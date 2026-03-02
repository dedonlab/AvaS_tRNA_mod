#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=10-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
###############################################################

output=AUA_stats_05202025.txt

echo -e "gid\ttotal_aua\ttotal_auh\ttotal_nnn\tratio_aua_syn_fold\taua_syn_gt_ratio\tratio_aua_fold\taua_gt_ratio\tnum_cds\tnum_cds_gt_aua_syn_ratio_per1000\tnum_cds_gt_aua_ratio_per100" > "${output}"

for file in $(ls /blue/lagard/yuanyifeng/PA14_24190/codon/CDS_codon/*_CDScodon.tsv); do
  gid=g"$(basename "$file" | sed 's/.PATRIC.ffn_CDScodon.tsv//')"
  total_aua=$(awk -F '\t' 'FNR>1{print $14}' "${file}" | paste -sd+ | bc)
  total_auc=$(awk -F '\t' 'FNR>1{print $15}' "${file}" | paste -sd+ | bc)
  total_aut=$(awk -F '\t' 'FNR>1{print $17}' "${file}" | paste -sd+ | bc)
  total_auh=$(echo "${total_aua} + ${total_auc} + ${total_aut}" | bc)
  ratio_aua_syn_fold=$(echo "scale=6; ${total_aua} / ${total_auh}" | bc)
  total_nnn=$(awk -F '\t' 'FNR==2{print $73}' "${file}")
  ratio_aua_fold=$(echo "scale=6; ${total_aua} / ${total_nnn}" | bc)
  aua_syn_gt_ratio=$(awk -F '\t' -v r=${ratio_aua_syn_fold} 'FNR>1&&($14+$15+$17)>0&&$14/($14+$15+$17)>r{print $1}' "${file}" | wc -l)
  aua_gt_ratio=$(awk -F '\t' -v r=${ratio_aua_fold} 'FNR>1&&($14/$68)>r{print $1}' "${file}" | wc -l)
  num_cds=$(awk -F '\t' 'FNR>1{print $1}' "${file}" | wc -l)
  num_cds_gt_aua_syn_ratio_per1000=$(echo "scale=6; ${aua_syn_gt_ratio} / ${num_cds} *1000" | bc)
  num_cds_gt_aua_ratio_per100=$(echo "scale=6; ${aua_gt_ratio} / ${num_cds} *100" | bc)
  echo -e "${gid}\t${total_aua}\t${total_auh}\t${total_nnn}\t${ratio_aua_syn_fold}\t${aua_syn_gt_ratio}\t${ratio_aua_fold}\t${aua_gt_ratio}\t${num_cds}\t${num_cds_gt_aua_syn_ratio_per1000}\t${num_cds_gt_aua_ratio_per100}" >> "${output}"
done
