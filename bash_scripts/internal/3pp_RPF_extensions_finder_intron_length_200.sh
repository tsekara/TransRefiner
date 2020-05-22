#!/bin/bash


# bash 3pp_extensions_finder_intron_length_100.sh <(grep "gene" /Users⁩/sekaran⁩/ownCloud⁩/conda_xtender⁩/gmap_dir⁩/GMAP_internal.gff3 | cut -f1,4,5,7,9| awk -F"=" '$4=$1' OFS="\t" | cut -f1,2,3,4,7 | awk '{print $5,$1,$2,$3,$4}' OFS="\t" |  awk 'NR==FNR {a[$1] = $2; next} { $6=a[$1] ;print $0}' OFS="\t" 3pp_ids_strand.txt  - | awk '{print $2,$3,$4,$5,$1,$6}' OFS="\t") > Extensions_3pp_intron_length_200.bed



#time_stamp=$(date +%Y-%m-%d-%T)
#mkdir Transcript_Scaffold_bed_files

mkdir -p Transcript_Scaffold_bed_files

filename="$1"

while read line;

do

scaffold=$(echo $line | awk '{print $1}')
transcript=$(echo $line | awk '{print $5}')
transcript_start=$(echo $line | awk '{print $2}')
transcript_end=$(echo $line | awk '{print $3}')
genome_strand=$(echo $line | awk '{print $4}')
transcriptome_strand=$(echo $line | awk '{print $6}')

#transcript="dd_Smed_v6_10468_0_1"
#scaffold="scaffold16385"

samtools view -@ 40 -F 4 -bS -o ./Transcript_Scaffold_bed_files/$scaffold.bam  /Users/sekaran/ownCloud/conda_xtender/Tophat2_alignment/RPF_direct_alignment_non-splice.bam $scaffold

bam2bed < ./Transcript_Scaffold_bed_files/$scaffold.bam | awk '!($8 ~ /N/)' | sort-bed -| cut -f 1,2,3,6 > ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed

rm -rf ./Transcript_Scaffold_bed_files/$scaffold.bam


grep $scaffold /Users/sekaran/ownCloud/conda_xtender/gmap_dir/GMAP_internal.gff3| grep "gene" | grep $transcript | cut -f1,4,5,7,9  | awk '{OFS=FS="\t"}  {if($4 == "+") {print $1,$3-1,$3,$4,$5} else if($4 =="-") {print $1,$2,$2+1,$4,$5}}' - > ./Transcript_Scaffold_bed_files/$transcript.bed


if grep -q "-" ./Transcript_Scaffold_bed_files/$transcript.bed

then

paste <(bedops --range 100 -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold /Users/sekaran/ownCloud/conda_xtender/gmap_dir/GMAP_internal.gff3 | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - ./Transcript_Scaffold_bed_files/$transcript.bed | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold /Users/sekaran/ownCloud/conda_xtender/gmap_dir/GMAP_internal.gff3 | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - | bedtools coverage -counts -a - -b ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed | awk '{ if ($4 > 0 ) {print $0} }') <(echo $genome_strand) <(echo $transcript) <(echo $transcript_start) <(echo $transcript_end) <(echo $transcriptome_strand) | awk -F"\t" -v OFS="\t" '$6!="" {Genome_strand=$5; transcript=$6; trans_start=$7; trans_end=$8; transcriptome_strand=$9} $0!="" {$5=Genome_strand; $6=transcript; $7=trans_start; $8=trans_end; $9=transcriptome_strand; print $0}' | awk '{if($2 < $7) print $0}' | awk 'NR>1{print last} {last=$0} END{$0=last;if($3>$7)$3=$7-1;print}' OFS='\t' - | awk '{ if($2-$3==0) $2=$3-1; print $0}' OFS='\t' | grep "scaffold" -


else

paste <(bedops --range 100 -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold /Users/sekaran/ownCloud/conda_xtender/gmap_dir/GMAP_internal.gff3 | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - ./Transcript_Scaffold_bed_files/$transcript.bed | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold /Users/sekaran/ownCloud/conda_xtender/gmap_dir/GMAP_internal.gff3 | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - | bedtools coverage -counts -a - -b ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed | awk '{ if ($4 > 0 ) {print $0} }') <(echo $genome_strand) <(echo $transcript) <(echo $transcript_start) <(echo $transcript_end) <(echo $transcriptome_strand) | awk -F"\t" -v OFS="\t" '$6!="" {Genome_strand=$5; transcript=$6; trans_start=$7; trans_end=$8; transcriptome_strand=$9} $0!="" {$5=Genome_strand; $6=transcript; $7=trans_start; $8=trans_end; $9=transcriptome_strand; print $0}' | awk '{if($2 > $7) print $0}'| awk 'NR==1{if($2<$8)$2=$8+1} 1' OFS='\t' - |awk '{ if($2-$3==0) $2=$3-1; print $0}' OFS='\t'| grep "scaffold" -

fi

done < "$filename"
