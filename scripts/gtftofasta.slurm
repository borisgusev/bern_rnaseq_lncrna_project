#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name='gtf_to_fasta'
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=%j_gtf.output
#SBATCH --error=%j_gtf.error

# add stringtie module
module add UHTS/Assembler/cufflinks/2.2.1

MERGEDGFF=/data/courses/rnaseq/lncRNAs/Project2/borisgusev/t_assembly/stringtie_output/merged.gtf
GENOMEFASTA=/data/courses/rnaseq/lncRNAs/Project2/references/genome_sequence/GRCh38.primary_assembly.genome.fa.gz

gunzip < $GENOMEFASTA > genome.fa

gffread $MERGEDGFF -w merged_exons.fa -g genome.fa

rm genome.fa
