#!/bin/bash

function show_usage ()
{
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -b|--bed  Extensions \n"
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
         bed_file="$1"
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
echo -e "$line" > ./Transcript_Scaffold_bed_files/temp.txt
awk '{if($4>=5) print $0}'| awk '{if($5=="-") $2=$2-1}1' OFS="\t" -
done < "$bed_file"
done
