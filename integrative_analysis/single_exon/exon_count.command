awk '{print }' merged.gtf | uniq -c | grep 'exon' | awk '{print }' > exon_count.txt
