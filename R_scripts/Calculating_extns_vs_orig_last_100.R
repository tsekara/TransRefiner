
## in bash
grep "gene" GMAP_Extended_last_100nts_deleted_intron_length_200.gff3 | awk -F"=" '$1=$1' OFS="\t" | cut  -f1,4,5,7,11 
> GMAP_Orig_last_100nts_deleted_intron_length_200.bed

grep "gene" GMAP_Last_100nts_deleted_Smed_complete_633_cov_0.9_idn_0.9.gff3 | awk -F"=" '$1=$1' OFS="\t" | cut  -f1,4,5,7,11
> GMAP_Orig_last_100nts_deleted_intron_length_200.bed


orig_last100nts=read.delim("GMAP_Orig_last_100nts_deleted_intron_length_200.bed",header=F)
colnames(orig_last100nts)=c("Chr","Orig_Start","Orig_End","Orig_Strand","Gene")

extn_last_100=read.delim("GMAP_Extn_last_100nts_deleted_intron_length_200.bed",header=F)
colnames(extn_last_100)=c("Chr","Extn_Start","Extn_End","Extn_Strand","Gene")

#Merging Original and Extended regions
merged_orig_extn=merge(orig_last100nts,extn_last_100,by="Gene")
#Missing codns_predicted
missing_codon_predicted=dim(merged_orig_extn[(merged_orig_extn$Extn_Strand == "+" & merged_orig_extn$Extn_End >  merged_orig_extn$Orig_End
                     |merged_orig_extn$Extn_Strand == "-" & merged_orig_extn$Extn_Start < merged_orig_extn$Orig_Start), ])[1] 

##Mising codon unpredicted
missing_codon_unpredicted=dim(merged_orig_extn[(merged_orig_extn$Extn_Strand == "+" & merged_orig_extn$Extn_End <  merged_orig_extn$Orig_End
                      |merged_orig_extn$Extn_Strand == "-" & merged_orig_extn$Extn_Start > merged_orig_extn$Orig_Start), ] )[1]