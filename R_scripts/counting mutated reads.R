
######################################## Aggragating based on the varaiant calling results.txt from GATK ############################################

library(data.table)
sam=read.delim("/home/tsekara/David_Schmidt/Bowtie2_GATK/Sample_59_negcon/Sample_59_negon_variant_calling_results.txt")
setDT(sam)[, list(AtoG_counts=sum(REF=='A' & ALT=='G'),
                  GtoA_counts = sum(REF=='G' & ALT=='A'),
                  TtoC_counts = sum(REF=='T' & ALT=='C'),
                  CtoT_counts = sum(REF=='C' & ALT=='T')), by =  CHROM]


################################## Reformating bases from bam2read counts ####################################### 


bam_read=read.delim("sample_bamread_counts.csv",header=F)
colnames(bam_read)=c("ID","Pos","Ref","Depth","A","C","G","T","N")
bam_read$A_counts = gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=bam_read$A)
bam_read$C_counts = gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=bam_read$C)
bam_read$G_counts = gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=bam_read$G)
bam_read$T_counts = gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=bam_read$T)
bam_read$A=NULL
bam_read$G=NULL
bam_read$T=NULL
bam_read$C=NULL
bam_read$N=NULL
bam_read$X=NULL


################################# Counting mutated bases on the reads from Varscan2 output #############################

sample_60=read.delim("Sample_60_bam_readcounts.txt",header=F)
sample_60=sample_60[c(-1,-2,-3),]
sample_60_org=sample_60
sample_60$V7 <-as.factor(sub("^(DEL|INS).*", "", sample_60$V7))
sample_60$V6 <-as.factor(sub("^(DEL|INS).*", "", sample_60$V6))
sapply(sample_60,class)
sample_60$Ref_base=as.factor(gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\1",x=sample_60$V6))
sample_60$Ref_base_counts=as.factor(gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=sample_60$V6))
sample_60$Alt_base=as.factor(gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\1",x=sample_60$V7))
sample_60$Alt_base_counts=as.factor(gsub(pattern="(^[A|G|T|C]+):([0-9]+):(.+$)",replacement="\\2",x=sample_60$V7))
sapply(sample_60,class)
sample_60$V3=NULL
sample_60$V4=NULL
sample_60$V5=NULL
sample_60$V6=NULL
sample_60$V7=NULL
colnames(sample_60)[1]="ID"
colnames(sample_60)[2]="Pos"
head(sample_60)
sapply(sample_60,class)
sample_60$Alt_base_counts=as.numeric(as.character(sample_60$Alt_base_counts))
aggregate(Alt_base_counts ~ ID + Ref_base + Alt_base, data=sample_60, sum)
sample_60$Move <- with(sample_60, paste0(Ref_base,"_to_",Alt_base))
sample_60_summary=spread(aggregate(Alt_base_counts ~ ID + Move, data=sample_60, sum), Move, Alt_base_counts)
head(sample_60_summary)
write.table(sample_60_summary,file="sample_60_summary.txt",sep="\t",row.names=F,quote=F)



















