# Run samtools view on to filter out the reads which match "specified" read length
samtools view ./Novel_ORFs_Prediction/RPF_raw_fastq/vs_complete_smed_ddV6_transcriptome/Sorted_GFP_RPF_linker_free_trimmed_1.fastq.bam  | gawk '(! and(4, $2))' | awk '{if (length($10)==30) print $3"\t",$4)}' | uniq -c > output.txt
### The above command provides output as follows:

6 dd_Smed_v6_10001_0_1	 240
3 dd_Smed_v6_10003_0_1	 355
11 dd_Smed_v6_10003_0_1	 869
9 dd_Smed_v6_10005_0_4	 845
11 dd_Smed_v6_10006_0_1	 134
14 dd_Smed_v6_10007_0_1	 991
13 dd_Smed_v6_10010_0_1	  98
3 dd_Smed_v6_10014_0_1	  51
1 dd_Smed_v6_1001_0_1	    56
# Later rearrange the columns like following: 
Id	Position	Counts
dd_Smed_v6_10001_0_1	240	 6
dd_Smed_v6_10003_0_1	355	 3
dd_Smed_v6_10003_0_1	869	11
dd_Smed_v6_10005_0_4	845	 9
dd_Smed_v6_10006_0_1	134	11
dd_Smed_v6_10007_0_1	991	14
dd_Smed_v6_10010_0_1	98	13
dd_Smed_v6_10014_0_1	51	 3
dd_Smed_v6_1001_0_1	  56	 1

#The modified above table is stored as text file and later imported into R.

#### To complement the above table, CDS ranges for every transcript provided in gff3 file is obtained using following command:

grep "CDS" smed_dd_v6.fasta.transdecoder.gff3 | awk '{$6=$5-$4; print $1"\t"$4"\t"$5"\t"$6}' | awk 'FNR != 1 && prev != $1 {print maxRecord;prev=$1;maxCalc=0;}{calc=$4;if(calc>maxCalc){maxRecord = $0;maxCalc = calc;}}END{print maxRecord}' > only_cds_gff3.txt

################## The above generated only_cds_gff3.txt and output.txt should be merged based on the 


################### merge the above output with only_cds_gff3.txt file(file contains CDS ranges) based on "Id" column ############################
cds_counts_len=merge(output.txt,only_cds_gff3.txt,by="Id")
##############################################################################################################################################

# Input data frame must have following columns in following order
# 1st Column= Id
# 2nd Column= Start_position
# 3rd Column= End_position
# 4th Column= Nucleotide_position
# 5th Column= Counts
sub_codon_phasing=function(df)
{
  #indexing=((df$Pos - df$Start)%%3 + 1) * (df$Start < df$Pos) * (df$End > df$Pos)
  indexing=((df[,4] - df[,2])%%3 + 1) * (df[,2] < df[,4]) * (df[,3] > df[,4])
  tab=matrix(0,nrow(df),3)
  #for(i in 1:3) tab[indexing==i,i] <- df$Hits[indexing==i]
  for(i in 1:3) tab[indexing==i,i] <- df[,5][indexing==i]
 # results=aggregate(tab,list(df$Id),FUN=sum)
  results=aggregate(tab,list(df[,1]),FUN=sum)
  return(results)
}

############################### generating bar plots fror frames ####################################

plot_ly(
  x = c("Frame 1", "Frame 2", "Frame 3"),
  y = c(92290, 88314, 73030),
  name = "Thilli_Custom_Codon_periodicity",
  type = "bar",color = c("Frame 1","Frame 2","Frame 3")
)
