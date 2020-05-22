#!/bin/bash



# The input for the "5_pp_extensions_appender_v3.sh" is the output from the script "5_prime_partial_transcript_extender_version3.sh" which contains the "Scaffold" "Extension_start" "Extensiopn_end" "Counts" "Genome_strand" "Transcript" "Transcript_start" "Transcript_End" "Transcriot_Orientation". SO we make use of the file and retain only unique entries based on the transcripts along with its Genome orientation and Transcript orientation like following

##########  bash internal_extensions_appender.sh <(awk 'NR==FNR { seen[$6]++; next } seen[$6]' /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_five_prime_extensions.bed /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_three_prime_extensions.bed | cut -f1,5,6,9 | sort -u -t$'\t' -k3,3) > internal_common_5pp_3pp.fasta

time_stamp=$(date +%Y-%m-%d-%T)

mkdir internal_$time_stamp

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

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_five_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_three_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

grep ">" ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 /Users/sekaran/ownCloud/conda_xtender/transcriptome_fasta/dd_Smed_v6.pcf.contigs.fasta - > ./internal_$time_stamp/"$transcript"_original.fasta

paste ./internal_$time_stamp/"$transcript"_3pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

#paste ./internal_$time_stamp/"$transcript"_3pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "+") && ($transcript_orientation == "+") ]]
then

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_five_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_three_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

grep ">" ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 /Users/sekaran/ownCloud/conda_xtender/transcriptome_fasta/dd_Smed_v6.pcf.contigs.fasta - > ./internal_$time_stamp/"$transcript"_original.fasta

paste ./internal_$time_stamp/"$transcript"_5pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_3pp_extension.fasta  | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

#paste ./internal_$time_stamp/"$transcript"_5pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "-") && ($transcript_orientation == "+") ]]
then

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_five_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2  | fastx_reverse_complement -i - > ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_three_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -   | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2 | fastx_reverse_complement -i - > ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

grep ">" ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 /Users/sekaran/ownCloud/conda_xtender/transcriptome_fasta/dd_Smed_v6.pcf.contigs.fasta - > ./internal_$time_stamp/"$transcript"_original.fasta

paste ./internal_$time_stamp/"$transcript"_5pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_3pp_extension.fasta  | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

#paste ./internal_$time_stamp/"$transcript"_5pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

######################################################################################################################################################################################################
elif [[ ($genome_orientation == "+") && ($transcript_orientation == "-") ]]
then

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_five_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2 | fastx_reverse_complement -i - > ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

grep -w "$transcript" /Users/sekaran/ownCloud/conda_xtender/internal_extensions/Cleaned_three_prime_extensions.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi /Users/sekaran/ownCloud/conda_xtender/Genome_fasta/SmedAsxl_genome_v1.1.fa -bed - | perl -p -e 's/:.*//g' - | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2 | fastx_reverse_complement -i - > ./internal_$time_stamp/"$transcript"_3pp_extension.fasta

grep ">" ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 /Users/sekaran/ownCloud/conda_xtender/transcriptome_fasta/dd_Smed_v6.pcf.contigs.fasta - > ./internal_$time_stamp/"$transcript"_original.fasta

paste ./internal_$time_stamp/"$transcript"_3pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_5pp_extension.fasta | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

#paste ./internal_$time_stamp/"$transcript"_3pp_extension.fasta ./internal_$time_stamp/"$transcript"_original.fasta ./internal_$time_stamp/"$transcript"_5pp_extension.fasta

else

echo Invalid Entry:"$transcript"

fi
done < "$filename"
