#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="map_reads"
#SBATCH --nodes=1
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=map_%j.output
#SBATCH --error=map_%j.error

# add STAR module
module add UHTS/Aligner/STAR/2.7.9a

# assign directory variables
G_INDEX='genome_index'
FWD=(data/*R1*)
BWD=(data/*R2*)

# run STAR to 
# see https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf#section.3
STAR --runThreadN 4 \
     --genomeDir $G_INDEX \
     --readFilesIn ${FWD[$SLURM_ARRAY_TASK_ID]} ${BWD[$SLURM_ARRAY_TASK_ID]}\
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ./mapping_output/$SLURM_ARRAY_TASK_ID/ \
     --outSAMstrandField intronMotif 

echo 'done'
