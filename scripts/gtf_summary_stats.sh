#!/bin/zsh

FILE=$1
OUTPUT=summary_stats.txt

echo 'Summary Stats for ' $1 > $OUTPUT
echo 'made file'

awk 'BEGIN {exons=0; transcripts=0} {if ($3 == "exon") {++exons}; if ($3 == "transcript") {++transcripts}} END {print "Exons\t" exons "\nTranscripts\t" transcripts}' $FILE >> $OUTPUT
echo 'counted exons and transcripts'

echo 'Genes\t' `awk -F '\t' '{print $9}' $FILE | awk -F ';' '{print $3}' | awk -F ' ' '$1=="gene_name" {print $2}' | awk -F '.' '{print $1}' | sort | uniq | wc -l` >> $OUTPUT
echo 'counted genes'

echo 'Novel genes\t' `awk -F '\t' '{print $9}' $FILE | awk -F ';' '{print $2}' | uniq | grep 'MSTRG' | wc -l` >> $OUTPUT
echo 'counted novel genes'

echo 'Single exon count\t' `awk '{print $3}' $FILE | uniq -c | awk '$2=="exon" && $1=="1" {print $0}' | wc -l` >> $OUTPUT
echo 'counted single exons'
