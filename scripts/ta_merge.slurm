#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="ta_merge"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=%j_merge.output
#SBATCH --error=%j_merge.error

# add stringtie module
module add UHTS/Aligner/stringtie/1.3.3b;

# Path Variables
REF_ANNOTATION='/data/courses/rnaseq/lncRNAs/Project2/references/annotations/gencode.v38.annotation.gtf.gz'
GTF_FILES='/data/courses/rnaseq/lncRNAs/Project2/borisgusev/t_assembly/stringtie_output/gtf_list.txt'

# Unzip Ref
gunzip < $REF_ANNOTATION > annotation.gtf

stringtie --merge -G annotation.gtf -o stringtie_output/merged.gtf $GTF_FILES

rm annotation.gtf

