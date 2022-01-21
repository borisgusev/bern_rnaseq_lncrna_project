#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="qc"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --output=qc_%j.output
#SBATCH --error=qc_%j.error

# add fastqc module
module add UHTS/Quality_control/fastqc/0.11.9

mkdir output

# fastqc -t 8 -o output data/*

singularity exec /software/singularity/containers/MultiQC-1.11-1.ubuntu20.sif multiqc `realpath output`

