
library(dplyr)
sample=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/Cleaned_Extensions_first_50nts_deleted_intron_length_70.bed",header=F)
colnames(sample)=c("Gene","Start","End","Counts","Chromosome", "Transcript_Start","Transcript_End","Strand")
sa=group_by(sample,Gene)
dim(as.data.frame(summarize(sa,Extn_start=min(Start),Extn_end=max(End))))
first_50nts=as.data.frame(summarize(sa,Extn_start=min(Start),Extn_end=max(End)))
 
# in linux #

grep "gene" GMAP_Smed_CDS_position_cov_.8_idn_.8.gff3 | awk -F"=" '$1=$1' OFS="\t"  - | cut -f4,5,7,11 | awk '{print $4,$1,$2,$3}' OFS="\t" -


dmel_cds_gene_start_end=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/dmel_CDS_start_end_strand.txt",header=F)
colnames(dmel_cds_gene_start_end)=c("Gene","Orig_Start","Orig_end","Strand")
dim(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
head(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
tail(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
first_50nts=merge(first_50nts,dmel_cds_gene_start_end,by="Gene")
dim(first_50nts[(first_50nts$Strand == "+" & first_50nts$start <  first_50nts$Orig_Start |first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end ), ])