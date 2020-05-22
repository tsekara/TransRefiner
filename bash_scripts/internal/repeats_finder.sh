#!/bin/bash

filename="$1"
while read line; 
do
transcript=$(echo $line | awk '{print $1}')
bl2seq -i internal_2018-04-03-09:53:11/"$transcript"_extension.fasta -j internal_2018-04-03-09:53:11/"$transcript"_original.fasta -D 1 -p blastn | grep "#" -v -

done < "$filename"


