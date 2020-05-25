
## in bash
grep "gene" GMAP_Extended_first_100nts_deleted_intron_length_200.gff3 | awk -F"=" '$1=$1' OFS="\t" | cut  -f1,4,5,7,11 
> GMAP_Extended_first_100nts_deleted_intron_length_200.bed

grep "gene" GMAP_First_100nts_deleted_Smed_complete_633_cov_0.9_idn_0.9.gff3 | awk -F"=" '$1=$1' OFS="\t" | cut  -f1,4,5,7,11
> GMAP_Orig_first_100nts_deleted_intron_length_200.bed


orig=read.delim("GMAP_Orig_first_100nts_deleted_intron_length_200.bed",header=F)
colnames(orig)=c("Chr","Orig_Start","Orig_End","Orig_Strand","Gene")
extn=read.delim("GMAP_Extended_first_100nts_deleted_intron_length_200.bed",header=F)
colnames(extn)=c("Chr","Extn_Start","Extn_End","Extn_Strand","Gene")
first_100nts=merge(orig,extn,by="Gene")
first_100_extn_greater_than_orig=first_100nts[(first_100nts$Extn_Strand == "+" & first_100nts$Extn_Start <  first_100nts$Orig_Start
                                            |first_100nts$Extn_Strand == "-" & first_100nts$Extn_End > first_100nts$Orig_End), ]
dim(first_100_extn_greater_than_orig)

first_100_extn_lesser_than_orig=dim(first_100nts[(first_100nts$Extn_Strand == "+" & first_100nts$Extn_Start >=  first_100nts$Orig_Start
 |first_100nts$Extn_Strand == "-" & first_100nts$Extn_End <= first_100nts$Orig_End), ])
dim(first_100_extn_lesser_than_orig)
