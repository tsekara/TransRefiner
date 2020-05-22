function show_usage ()
{
    printf "Usage: $0 [options [parameters]]\n"
    printf "\n"
    printf "Options:\n"
    printf " -b_5pp|--bed_5pp, path to 5pp filtered extensions\n"
    printf " -b_3pp|--bed_3pp, path to 5pp filtered extensions\n"
    printf " -gfa|--genome_fasta  Genome fasta\n"
    printf " -tfa|--transcriptome_fasta, Transcriptome fasta\n"
    printf " -o|--orientation, Filtered 5pp and 3pp transcripts orientation in genome and transcriptome \n"
    printf " -h|--help, Print help\n"
return 0
}


#Checks for usage
if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]] || [[ "$1" == "" ]] ;then
    show_usage
fi

#f [ "$#" -ne 3 ]; then
#    echo "Please ensure that all the params are provided"
#    #exit 1
#fi

#getting the repective arguments from the user

while [ ! -z "$1" ]; do
  case "$1" in
     --bed_5pp|-b_5pp)
         shift
         bed_5pp_file="$1"
         #echo "$bed_5pp_file"
         ;;
      --bed_3pp|-b_3pp)
         shift
         bed_3pp_file="$1"
         ;;
     --genome_fasta|-gfa)
         shift
         genome_fasta="$1"
        #echo "$genome_fasta"
         ;;
     --transcriptome_fasta|-tfa)
        shift
        transcriptome_fasta="$1"
        #echo "$transcriptome_fasta"
         ;;
     --orientation|-o)
       shift
       filtered_trans_orientation_file="$1"
       #echo "$orientation"
        ;;
      *)
        show_usage
        ;;
  esac
shift
# The input for the "3_prime_partial_extensions_appender.sh" is the output from the script "3_prime_partial_transcript_extender_version3.sh" which contains the "Scaffold" "Extension_start" "Extensiopn_end" "Counts" "Genome_strand" "Transcript" "Transcript_start" "Transcript_End" "Transcriot_Orientation". SO we make use of the file and retain only unique entries using "sort" command based on the transcripts along with its Genome orientation and Transcript orientation like following

##########  bash 3pp_extension_appender_v2.sh <(sort -u -t$'\t' -k 6,6 $bed_file_file | cut -f1,5,6,9)

#time_stamp=$(date +%Y-%m-%d-%T)

mkdir -p internal_extension_and_original_fasta

  while IFS= read -r line; do

transcript=$(echo $line | awk '{print $3}')
genome_orientation=$(echo $line | awk '{print $2}')
transcript_orientation=$(echo $line | awk '{print $4}')

if [[ ($genome_orientation == "-") && ($transcript_orientation == "-") ]] # comparing  genome and transcriptome orientation
then

# Extracting the transcript stored in "$transcript" from extensions and rearranging the bed file so that transcript_id comes in the fourth column which manadate for "bedtools getfasta -name" paramater.
# "perl -p -e 's/:.*//g'" removes the unwanted charatcers after the fasta header.
# "awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' -" Concatanates the fasta sequecnes from same id into single fasta
# "sed '/^\s*$/d'" - removes the blank lines created above, if any

grep -w "$transcript"  $bed_5pp_file | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi $genome_fasta -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta
#grep -w "$transcript"  $bed_5pp_file | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}'

#| bedtools getfasta -name -fi $genome_fasta -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta

grep -w "$transcript" $bed_3pp_file | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi $genome_fasta -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_extension_and_original_fasta/"$transcript"_3pp_extension.fasta

grep ">" ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 $transcriptome_fasta - > ./internal_extension_and_original_fasta/"$transcript"_original.fasta

paste ./internal_extension_and_original_fasta/"$transcript"_3pp_extension.fasta ./internal_extension_and_original_fasta/"$transcript"_original.fasta ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

#paste ./internal_extension_and_original_fasta/"$transcript"_3pp_extension.fasta ./internal_extension_and_original_fasta/"$transcript"_original.fasta ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta

######################################################################################################################################################################################################

elif [[ ($genome_orientation == "+") && ($transcript_orientation == "+") ]]
then

grep -w "$transcript"  $bed_5pp_file | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi $genome_fasta -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta

grep -w "$transcript" $bed_3pp_file | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9}' | bedtools getfasta -name -fi $genome_fasta -bed - | perl -p -e 's/:.*//g' -  | awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' - | tail -n+2> ./internal_extension_and_original_fasta/"$transcript"_3pp_extension.fasta

grep ">" ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta | perl -p -e 's/>//g' - | seqtk subseq -l 100000 $transcriptome_fasta - > ./internal_extension_and_original_fasta/"$transcript"_original.fasta

paste ./internal_extension_and_original_fasta/"$transcript"_5pp_extension.fasta ./internal_extension_and_original_fasta/"$transcript"_original.fasta ./internal_extension_and_original_fasta/"$transcript"_3pp_extension.fasta  | awk '{if (NR%2==0) {print $1 $2 $3} else {print $2}}'

else

echo Invalid Entry:"$transcript"

fi
    done < "$filtered_trans_orientation_file"
done
