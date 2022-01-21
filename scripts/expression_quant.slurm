#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="exprquant"
#SBATCH --nodes=1
#SBATCH --array=0-11
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH --output=%j_exprq.output
#SBATCH --error=%j_exprq.error

# add stringtie module
module add UHTS/Analysis/kallisto/0.46.0

K_INDEX=/data/courses/rnaseq/lncRNAs/Project2/borisgusev/expression_quant/kallisto_index
FWD=(/data/courses/rnaseq/lncRNAs/Project2/fastq/*R1*)
BWD=(/data/courses/rnaseq/lncRNAs/Project2/fastq/*R2*)

mkdir -p kallisto_output/$SLURM_ARRAY_TASK_ID

kallisto quant -i $K_INDEX -o kallisto_output/$SLURM_ARRAY_TASK_ID -t 2 -b 10 --rf-stranded ${FWD[$SLURM_ARRAY_TASK_ID]} ${BWD[$SLURM_ARRAY_TASK_ID]}
