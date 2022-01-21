#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="expression_quant"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=%j_exprq.output
#SBATCH --error=%j_exprq.error

# add stringtie module
module add UHTS/Analysis/kallisto/0.46.0

MERGEDFASTA=/data/courses/rnaseq/lncRNAs/Project2/borisgusev/expression_quant/merged_exons.fa

kallisto index -i kallisto_index $MERGEDFASTA
