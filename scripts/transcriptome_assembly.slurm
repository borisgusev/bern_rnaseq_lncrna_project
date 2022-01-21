#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="transcriptome_assembly"
#SBATCH --nodes=1
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=ta_%j.output
#SBATCH --error=ta_%j.error

# add stringtie module
module add UHTS/Aligner/stringtie/1.3.3b;

# gtf reference
REFERENCE='/data/courses/rnaseq/lncRNAs/Project2/references/annotations/gencode.v38.annotation.gtf'
BAM_FILE=(/data/courses/rnaseq/lncRNAs/Project2/borisgusev/mapping/star_output/$SLURM_ARRAY_TASK_ID/*.bam)

# unzip the zipped reference, keeping the original zipped
gunzip < $REFERENCE.gz > temp_$SLURM_ARRAY_TASK_ID

stringtie -o stringtie_output/$SLURM_ARRAY_TASK_ID.gtf -p 1 --rf -G temp_$SLURM_ARRAY_TASK_ID $BAM_FILE

# delete unzipped reference
rm temp_$SLURM_ARRAY_TASK_ID
