#!/bin/bash

#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="count_reads"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB

INPUT_FILES="$@"

for FILE in $INPUT_FILES
do
    FILENAME=`basename $FILE`
    LINES=`zcat $FILE | wc -l`
    READS=$(($LINES / 4))
    echo $FILENAME $READS >> count_reads.output
done

echo 'done'