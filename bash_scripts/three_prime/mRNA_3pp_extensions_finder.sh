#!/bin/bash

function show_usage ()
{
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -gff3|--annotation, path to GMAP gf3 file for 3pp\n"
    printf " -b|--bam  Topat2 alignment file\n"
    printf " -o|--orientation, Transcript orientation in genome and transcriptome\n"
    printf " -l|--search_length, Standard search length\n"
    printf " -h|--help, Print help\n"
return 0
}

#Checks for usage
if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]] || [[ "$1" == "" ]] ;then
    show_usage
fi

#getting the repective arguments from the user

while [ ! -z "$1" ]; do
  case "$1" in
       --annotation|-gff3)
         shift
         annotation="$1"
         #echo "$bed_file"
         ;;
       --bam|-b)
         shift
         bam="$1"
        #echo "$genome_fasta"
         ;;
       --orientation|-o)
       shift
       orientation_file="$1"
       #echo "$orientation"
           ;;
       --search_length|-l)
       shift
       search_length="$1"
           ;;
     *)
      show_usage
      ;;
  esac
shift

mkdir -p Transcript_Scaffold_bed_files


while IFS= read -r line;

do

  scaffold=$(echo $line | awk '{print $1}')
  transcript=$(echo $line | awk '{print $5}')
  transcript_start=$(echo $line | awk '{print $2}')
  transcript_end=$(echo $line | awk '{print $3}')
  genome_strand=$(echo $line | awk '{print $4}')
  transcriptome_strand=$(echo $line | awk '{print $6}')

  #transcript="dd_Smed_v6_10468_0_1"
  #scaffold="scaffold16385"

  samtools view -@ 40 -F 4 -bS -o ./Transcript_Scaffold_bed_files/$scaffold.bam $bam $scaffold

  bam2bed < ./Transcript_Scaffold_bed_files/$scaffold.bam |  sort-bed - | cut -f 1,2,3,6 > ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed

  bamToBed -i ./Transcript_Scaffold_bed_files/$scaffold.bam -split | cut -f1,2,3,4,6 | sort-bed - > ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed


  rm ./Transcript_Scaffold_bed_files/$scaffold.bam


  grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9  | awk '{OFS=FS="\t"}  {if($4 == "+") {print $1,$3-1,$3,$4,$5} else if($4 =="-") {print $1,$2,$2+1,$4,$5}}' - > ./Transcript_Scaffold_bed_files/$transcript.bed


  if grep -q "-" ./Transcript_Scaffold_bed_files/$transcript.bed

  then

  paste <(bedops --range $((search_length / 2)) -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - ./Transcript_Scaffold_bed_files/$transcript.bed | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - )  > non_spliced_extensions.bed

  paste <(bedops --range $((search_length / 2)) -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - non_spliced_extensions.bed  | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - | bedtools coverage -counts -a - -b ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed | awk '{ if ($4 > 0 ) {print $0} }') <(echo $genome_strand) <(echo $transcript) <(echo $transcript_start) <(echo $transcript_end) <(echo $transcriptome_strand) | awk -F"\t" -v OFS="\t" '$6!="" {Genome_strand=$5; transcript=$6; trans_start=$7; trans_end=$8; transcriptome_strand=$9} $0!="" {$5=Genome_strand; $6=transcript; $7=trans_start; $8=trans_end; $9=transcriptome_strand; print $0}' | awk '{if($2 < $7) print $0}' | awk 'NR>1{print last} {last=$0} END{$0=last;if($3>$7)$3=$7-1;print}' OFS='\t' - | awk '{ if($2-$3==0) $2=$3-1; print $0}' OFS='\t' | grep "scaffold" -



  else

  paste <(bedops --range $((search_length / 2)) -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - ./Transcript_Scaffold_bed_files/$transcript.bed | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Sorted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - )  > non_spliced_extensions.bed

  paste <(bedops --range $((search_length / 2)) -m <(bedops -n 100% ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) | bedops -e 1 - non_spliced_extensions.bed | bedops -e <(bedops -n 100% ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed <(grep $scaffold $annotation | grep "gene" | grep $transcript | cut -f1,4,5,7,9)) - | bedops -m - | bedtools coverage -counts -a - -b ./Transcript_Scaffold_bed_files/Splitted_$scaffold.bed | awk '{ if ($4 > 0 ) {print $0} }') <(echo $genome_strand) <(echo $transcript) <(echo $transcript_start) <(echo $transcript_end) <(echo $transcriptome_strand) | awk -F"\t" -v OFS="\t" '$6!="" {Genome_strand=$5; transcript=$6; trans_start=$7; trans_end=$8; transcriptome_strand=$9} $0!="" {$5=Genome_strand; $6=transcript; $7=trans_start; $8=trans_end; $9=transcriptome_strand; print $0}' |  awk '{if($2 > $7) print $0}' | awk 'NR==1{if($2<$8)$2=$8+1} 1' OFS='\t' - |awk '{ if($2-$3==0) $2=$3-1; print $0}' OFS='\t'| grep "scaffold" -

  fi


done < "$orientation_file"

done
