#!/bin/bash



# The input for the "5_pp_extensions_appender_v3.sh" is the output from the script "5_prime_partial_transcript_extender_version3.sh" which contains the "Scaffold" "Extension_start" "Extensiopn_end" "Counts" "Genome_strand" "Transcript" "Transcript_start" "Transcript_End" "Transcriot_Orientation". SO we make use of the file and retain only unique entries based on the transcripts along with its Genome orientation and Transcript orientation like following

##########  bash 5pp_extensions_appender_v3.sh <(sort -u -t$'\t' -k 6,6 Cleaned_Extensions_mRNA_5pp_intron_length_200.bed| cut -f1,5,6,9)       

time_stamp=$(date +%Y-%m-%d-%T)

mkdir 5pp_extension_and_original_fasta_$time_stamp

filename="$1"

while read line; 

do

transcript=$(echo $line | awk '{print $3}')
genome_orientation=$(echo $line | awk '{print $2}')
transcript_orientation=$(echo $line | awk '{print $4}')

if [[ ($genome_orientation == "-") && ($transcript_orientation == "-") ]] # comparing  genome and transcriptome orientation
then

# Extracting the transcript stored in "$transcript" from extensions and rearranging the bed file so that transcript_id comes in the fourth column which manadate for "bedtools getfasta -name" paramater.
# "sed 's/:.*//'" removes the unwanted charatcers after the fasta header. 
# "awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' -" Concatanates the fasta sequecnes from same id into single fasta
# "sed '/^\s*$/d'" - removes the blank lines created above, if any




grep -w "$transcript" Cleaned_Extensions_mRNA_5pp_intron_length_200.bed| awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /home/tsekara/ownCloud/SmedAsxl_genome_v1.1.fa -bed - | sed 's/:.*//' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | sed '/^\s*$/d' -> ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta

grep ">" ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | sed 's/>//' - | seqtk subseq -l 100000 /home/tsekara/ownCloud/dd_Smed_v6.pcf.contigs.fasta - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta

paste ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | awk '{if (NR%2==0) {print $1 $2} else {print $1}}' 

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "+") && ($transcript_orientation == "+") ]]
then



grep -w "$transcript" Cleaned_Extensions_mRNA_5pp_intron_length_200.bed| awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /home/tsekara/ownCloud/SmedAsxl_genome_v1.1.fa -bed - | sed 's/:.*//' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | sed '/^\s*$/d' -> ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta

grep ">" ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | sed 's/>//' - | seqtk subseq -l 100000 /home/tsekara/ownCloud/dd_Smed_v6.pcf.contigs.fasta - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta

paste ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta | awk '{if (NR%2==0) {print $1 $2} else {print $1}}' 

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "-") && ($transcript_orientation == "+") ]]
then



grep -w "$transcript" Cleaned_Extensions_mRNA_5pp_intron_length_200.bed| awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /home/tsekara/ownCloud/SmedAsxl_genome_v1.1.fa -bed - | sed 's/:.*//' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | sed '/^\s*$/d' - | fastx_reverse_complement -i - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta

grep ">" ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | sed 's/>//' - | seqtk subseq -l 100000 /home/tsekara/ownCloud/dd_Smed_v6.pcf.contigs.fasta - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta

paste ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta | awk '{if (NR%2==0) {print $1 $2} else {print $1}}' 

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "+") && ($transcript_orientation == "-") ]]
then



grep -w "$transcript" Cleaned_Extensions_mRNA_5pp_intron_length_200.bed| awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /home/tsekara/ownCloud/SmedAsxl_genome_v1.1.fa -bed - | sed 's/:.*//' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | sed '/^\s*$/d' -| fastx_reverse_complement -i - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta

grep ">" ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | sed 's/>//' - | seqtk subseq -l 100000 /home/tsekara/ownCloud/dd_Smed_v6.pcf.contigs.fasta - > ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta

paste ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_original.fasta ./5pp_extension_and_original_fasta_$time_stamp/"$transcript"_extension.fasta | awk '{if (NR%2==0) {print $1 $2} else {print $1}}' 

######################################################################################################################################################################################################
else

echo Invlaid Entry:"$transcript"

fi
done < "$filename"
