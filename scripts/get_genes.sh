#!/usr/bin/env bash

# script to geta file with gene_id <-> transcript_id mapping

awk '{if ($3=="exon") print}' < $1 | cut -f 9 | perl -ne 'chomp; $gname=""; $gid=""; $tid=""; if ($_ =~ /gene_name\s+\"(\S+)\"\;/){$gname=$1}; if ($_ =~ /gene_id\s+\"(\S+)\"\;/){$gid=$1}; if ($_ =~ /transcript_id\s+\"(\S+)\"\;/){$tid=$1} print "$gname\t$gid\t$tid\n"' | sort | uniq > genename_gid_tid.tsv
