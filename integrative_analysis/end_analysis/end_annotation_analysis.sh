#!/bin/zsh

# find lsit of transcripts that have a 
TRANSCRIPTS=$1 
CAGECLUSTERS=$2
POLYASITES=$3

bedtools window -l 100 -r 0 -sw -u -a $TRANSCRIPTS -b $CAGECLUSTERS > $TRANSCRIPTS.withcage

bedtools window -l 0 -r 100 -sw -u -a $TRANSCRIPTS -b $POLYASITES > $TRANSCRIPTS.withpolya

echo 'done'
