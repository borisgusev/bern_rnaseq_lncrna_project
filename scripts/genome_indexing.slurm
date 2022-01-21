#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="genome_indexing"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2:00:00
#SBATCH --mem=40GB
#SBATCH --output=indexing_%j.output
#SBATCH --error=indexing_%j.error

# add STAR module
module add UHTS/Aligner/STAR/2.7.9a

# declare path variables
OUTPUT='genome_index'
GENOME='/data/courses/rnaseq/lncRNAs/Project2/references/genome_sequence/GRCh38.primary_assembly.genome.fa'

# unzip files whilst keeping the zipped version
gunzip < $GENOME.gz > $GENOME

# run STAR to generate the genome indicies
# see https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf#section.2
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $OUTPUT --genomeFastaFiles $GENOME

# remove the unzipped files
rm $GENOME

echo 'done'
