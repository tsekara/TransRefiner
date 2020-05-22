#!/bin/bash

#The input for this script is the output of "5pp_extensions_finder.sh".

filename="$1"
while read line;
do

scaffold=$(echo $line | awk '{print $1}')
echo -e "$line" > Transcript_Scaffold/temp.txt
plus=$(awk '$4=="+"' Transcript_Scaffold/Sorted_"$scaffold".bed | bedmap --count Transcript_Scaffold/temp.txt - )
minus=$(awk '$4=="-"' Transcript_Scaffold/Sorted_"$scaffold".bed | bedmap --count Transcript_Scaffold/temp.txt  -)
paste <(echo -e $line) <(echo $plus) <(echo $minus) | awk -v OFS="\t" '$1=$1' - | awk '{ if ($5 == "+" && $10 > $11) print $0 ; else if ($5 == "-" && $11 > $10 ) print $0 }' - | awk '{if($4>5) print $0}' | awk '{if($5=="-") $2=$2-1}1' OFS="\t" -
done < "$filename"

# scaffold=$(echo $line | awk '{print $1}') ----> Prints the current line of the input and later extracts the "scaffold" and stores in the variable called "scaffold"

#echo -e "$line" > ./Transcript_Scaffold/temp.txt --------> Prints the current line into a file called "./Transcript_Scaffold/temp.txt" in the working directory.

#plus=$(awk '$4=="+"' Sorted_"$scaffold".bed | bedmap --count ./Transcript_Scaffold/temp.txt - ) ---------------> Extracts all the reads in "+" orientation in that particular "scaffold" and later counts the number of reads within in the range specifed in ./Transcript_Scaffold/temp.txt and stores in a varibale called "plus"


#minus=$(awk '$4=="-"' Sorted_"$scaffold".bed | bedmap --count ./Transcript_Scaffold/temp.txt  -) ---------------> Extracts all the reads in "-" orientation in that particular "scaffold" and later counts the number of reads within in the range specifed in ./Transcript_Scaffold/temp.txt and stores in a varibale called "minus"

#paste <(echo -e $line) <(echo $plus) <(echo $minus) -------------> Prints the current line alogn with plus and minus in a spcae delim format


#sed 's/ /\t/g' - --------------------> As the above command generates the output in non tab delim format, the sed command replaces spcaes with tab delim format


# awk '{ if ($5 == "+" && $10 > $11) print $0 ; else if ($5 == "-" && $11 > $10 ) print $0 }' - ----------------> So if the genome orientation($5) is "+" and the "plus" column($10) is greater than "minus" column($11), exclude that line OR if the genome orientation($5) is "-" and the "plus" column($10) is less than "minus" column($11) exclude that line


#awk '{if($4>5) print $0}' -------------------------> retain ranges which has counts greater than 5
