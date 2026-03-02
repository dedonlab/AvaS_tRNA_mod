#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single core
#SBATCH --cpus-per-task=1             #number of CPUs
#SBATCH --time=7-00:00:00              # Time limit hrs:min:sec
#SBATCH --mem=4gb
###############################################################

dir_w=$1

<<'####'
# 1. covert to lowercase.
# create ffn prep folder.

#mkdir ${dir_w}/ffn_prep || true

module load seqkit/2.4.0

# -l, print sequences in lower case.
# -m, print sequences >= 10 aa (33 nts).

for raw_ffn in $(ls ${dir_w}/ffnfiles/666*.ffn); do
  bname=$(basename ${raw_ffn})
  seqkit seq -l -m 33 ${raw_ffn} > ${dir_w}/ffn_prep/${bname}
done

####

# 2. prepare ffn files.
for ffn in $(ls ${dir_w}/ffn_prep/666*ffn); do
  # prepare sequences.
  sed -i '/^>/ s/ .*$/#/g' $ffn      # replace space after fig number with #.
  sed -i ':a;N;$!ba;s/\n//g' $ffn    # delet all line breaks.
  sed -i 's/>/\n>/g' $ffn            # make each > a new line.
  sed -i 's/#/\n/g' $ffn             # make # a new line.
  ##sed -i '/^>/! s/^...//' $ffn       # delete start codon.
  sed -i '/^>/! s/...$//' $ffn       # delete stop codon.
  sed -i '/^>/! s/.\{3\}/& /g' $ffn  # insert space every triplet.
  sed -i '/^$/d' $ffn                # remove empty lines.
done

####
#<<'####'

# call py script to count codons.
for ffn in $(ls ${dir_w}/ffn_prep/666*ffn); do
  python 1_pycodon_count.py ${ffn} ${dir_w}
done


# call py script to count ata_n.
for ffn in $(ls ${dir_w}/ffn_prep/666*ffn); do
  python 1_count_CDS_AUA_N.py ${ffn} ${dir_w}
done

#### END ####
