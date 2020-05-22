#!/bin/bash
filename="$1"
bedtools intersect -a <(cut -f1,2,3,6 "$filename" | sort-bed -) -b <(grep "gene" ~/ownCloud/GMAP_vs_dd_v6_pcf_idn_0.7_cov_0.7.gff3 | awk -F"=" '$1=$1' OFS="\t" | cut  -f1,4,5,11) -wa -wb | sort -k4,4 | awk -F"\t" '!seen[$4,$8]++' - |  awk '{sub(/_[0-9]*$/, "", $4); sub(/_[0-9]*$/, "", $8); print}' OFS="\t" | awk '{if($4!=$8) print $0}' | cut -f4,8 | awk '!seen[$0]++' | cut -f1 | uniq | wc -l
