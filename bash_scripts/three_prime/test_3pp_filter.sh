#!/bin/bash

function show_usage ()
{
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -b|--bed  Extensions from "Extensions finder" script\n"
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
       --bed|-b)
         shift
         bed="$1"
        #echo "$bed"
         ;;
      *)
      show_usage
      ;;
  esac
shift

while IFS= read -r line;
do

scaffold=$(echo $line | awk '{print $1}')
echo -e "$line" > Transcript_Scaffold_bed_files/temp.txt
plus=$(awk '$4=="+"' Transcript_Scaffold_bed_files/Sorted_"$scaffold".bed | bedmap --count Transcript_Scaffold_bed_files/temp.txt - )
minus=$(awk '$4=="-"' Transcript_Scaffold_bed_files/Sorted_"$scaffold".bed | bedmap --count Transcript_Scaffold_bed_files/temp.txt  -)
paste <(echo -e $line) <(echo $plus) <(echo $minus) | awk -v OFS="\t" '$1=$1' - | awk '{ if ($5 == "+" && $10 > $11) print $0 ; else if ($5 == "-" && $11 > $10 ) print $0 }' - | awk '{if($4>5) print $0}' | awk '{if($5=="+") $2=$2-1}1' OFS="\t" - | sort -u -t$'\t' -k 6,6 -
done < "$bed"
done
