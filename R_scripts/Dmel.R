main = "egr6_3h vs gfp_3h, correct null model", xlab = "CORRECTED p-values")
table(egr6_3h_vs_gfp_3h[,"padj"] < 0.1)
table(egr6_3h_vs_gfp_3h[,"padj"] < 0.05)
plotMA(egr6_3h_vs_gfp_3h)
head(egr6_3h_vs_gfp_3h)
sig_egr6_3h_vs_gfp_3h=as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.1,])
write.table(sig_egr6_3h_vs_gfp_3h,"~/owncloud.gwdg.de/Frances's_Thilli_DESeq2_standard_pipeline/Sig_egr6_3h_vs_gfp_3h.txt",
sep="\t",row.names = T,quote=F)
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,]))
egr6_3h_vs_gfp_3h=results(dds, contrast = c("conditions","egr6_3h","gfp_3h"))
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,]))
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,]))
egr6_3h_vs_gfp_3h <- egr6_3h_vs_gfp_3h[ !is.na(egr6_3h_vs_gfp_3h$padj), ]
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,]))
as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,])[1]
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,])[1])
dim(as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,])[,1])
as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,])[,1]
as.data.frame(egr6_3h_vs_gfp_3h[egr6_3h_vs_gfp_3h$padj<.05,])[1]
egr6_18h_vs_gfp_18h=results(dds, contrast = c("conditions","egr6_18h","gfp_18h"))
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$padj), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$pvalue), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[, -which(names(egr6_18h_vs_gfp_18h) == "padj")]
FDR_egr6_18h_vs_gfp_18h <- fdrtool(egr6_18h_vs_gfp_18h$stat, statistic= "normal", plot = T)
sig_egr4_18h_vs_gfp_18h=as.data.frame(DESeq2Res[DESeq2Res$padj<.1,])
write.table(sig_egr4_18h_vs_gfp_18h,"~/owncloud.gwdg.de/Frances's_Thilli_DESeq2_standard_pipeline/Sig_egr4_18h_vs_gfp_18h.txt",
sep="\t",row.names = T,quote=F)
egr6_18h_vs_gfp_18h=results(dds, contrast = c("conditions","egr6_18h","gfp_18h"))
egr6_18h_vs_gfp_18h=results(dds, contrast = c("conditions","egr6_18h","gfp_18h"))
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$padj), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$pvalue), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[, -which(names(egr6_18h_vs_gfp_18h) == "padj")]
FDR_egr6_18h_vs_gfp_18h <- fdrtool(egr6_18h_vs_gfp_18h$stat, statistic= "normal", plot = T)
sig_egr6_18h_vs_gfp_18h=as.data.frame(egr6_18h_vs_gfp_18h[egr6_18h_vs_gfp_18h$padj<.1,])
write.table(sig_egr6_18h_vs_gfp_18h,"~/owncloud.gwdg.de/Frances's_Thilli_DESeq2_standard_pipeline/Sig_egr6_18h_vs_gfp_18h.txt",
sep="\t",row.names = T,quote=F)
sig_egr6_18h_vs_gfp_18h
FDR_egr6_18h_vs_gfp_18h
rm(egr6_18h_vs_gfp_18h)
rm(FDR_egr6_18h_vs_gfp_18h)
egr6_18h_vs_gfp_18h=results(dds, contrast = c("conditions","egr6_18h","gfp_18h"))
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$padj), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[ !is.na(egr6_18h_vs_gfp_18h$pvalue), ]
egr6_18h_vs_gfp_18h <- egr6_18h_vs_gfp_18h[, -which(names(egr6_18h_vs_gfp_18h) == "padj")]
FDR_egr6_18h_vs_gfp_18h <- fdrtool(egr6_18h_vs_gfp_18h$stat, statistic= "normal", plot = T)
egr6_18h_vs_gfp_18h
egr6_18h_vs_gfp_18h[,"padj"]  <- p.adjust(FDR_egr6_18h_vs_gfp_18h$pval, method = "BH")
sig_egr6_18h_vs_gfp_18h=as.data.frame(egr6_18h_vs_gfp_18h[egr6_18h_vs_gfp_18h$padj<.1,])
dim(sig_egr6_18h_vs_gfp_18h)
write.table(sig_egr6_18h_vs_gfp_18h,"~/owncloud.gwdg.de/Frances's_Thilli_DESeq2_standard_pipeline/Sig_egr6_18h_vs_gfp_18h.txt", sep="\t",row.names = T,quote=F)
source("https://bioconductor.org/biocLite.R")
biocLite("maSigPro")
library(maSigPro)
vignette("maSigPro")
browseVignettes("maSigPro")
log2(1)
log2(1.5)
e-0.5
log2(0.5)
ls()
test=read.delim("~/test.txt")
test
library(data.table)
setDT(df1)[,.(start = (start+end - start +1)[-.N],
end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1)[-.N],
end = (End +shift(Start, type='lead')-End-1)[-.N] ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1)[-.N],
end = (End +shift(Start, type='lead')-End-1)[-.N], transcript=Transcript_Id ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1)[-.N],end = (End +shift(Start, type='lead')-End-1)[-.N], transcript=Transcript_Id[-.N] ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1),end = (End +shift(Start, type='lead')-End-1)[-.N], transcript=Transcript_Id[-.N] ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1),end = (End +shift(Start, type='lead')-End-1)[-.N], ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1),end = (End +shift(Start, type='lead')-End-1), ) , by = Scaffolds]
setDT(test)[,.(start = (Start+End - Start +1),end = (End +shift(Start, type='lead')-End-1)) , by = Scaffolds]
setDT(df1)[,.(start = (end+1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group]
setDT(test)[,.(start = (End+1)[-.N],end = (End +shift(Start, type='lead')-End-1)[-.N] ) , by = group]
setDT(test)[,.(start = (End+1)[-.N],end = (End +shift(Start, type='lead')-End-1)[-.N] ) , by = Scaffolds]
test
source("https://bioconductor.org/biocLite.R")
biocLite("RiboProfiling")
load("/home/tsekara/R_Scripts/RiboProfiling.RData")
ls()
aln <- readGAlignments(
BamFile("~/owncloud.gwdg.de/Remaping_with_extended_transcriptome/RPF/entire_extended_dd_Smed_v6/Extended.bam")
)
library(GenomicAlignments#)
)
aln <- readGAlignments(
BamFile("~/owncloud.gwdg.de/Remaping_with_extended_transcriptome/RPF/entire_extended_dd_Smed_v6/Extended.bam")
)
alnGRanges=readsToReadStart(aln)
alnGRanges=readsToStartOrEnd(aln, what="start")
library(RiboProfiling)
alnGRanges=readsToStartOrEnd(aln, what="start")
alnGRanges=  alnGRanges[which(!is.na(match(alnGRanges$score,28:31)))]
library(GenomicRanges)
library(GenomicFeatures)
mytxdb=makeTxDbFromGFF("~/owncloud.gwdg.de/Remaping_with_extended_transcriptome/RPF/entire_extended_dd_Smed_v6/Extended_RPFs.fasta.transdecoder.gff3")
cds=cdsBy(mytxdb,by="tx",use.names= TRUE)
exonGRanges=exonsBy(mytxdb,by="tx",use.names= TRUE)
cdsPosTransc=orfRelativePos(cds, exonGRanges)
countsDataCtrl1=countShiftReads(exonGRanges=exonGRanges[names(cdsPosTransc)],cdsPosTransc=cdsPosTransc,alnGRanges=alnGRanges,shiftValue=-14)
listCountsPlots=countsPlot(list(countsDataCtrl1[[1]]),grep("_counts$",colnames(countsDataCtrl1[[1]])),1)
cdsBy()
?cdsBy
head(Cds)
head(cds)
source("https://bioconductor.org/biocLite.R")
biocLite("derfinder")
three=read.delim("~/Transcriptome_Enhancement/RPFs/Pooled/pooled_master_wo_H2Bs/three_prime_counts.txt")
five=read.delim("~/Transcriptome_Enhancement/RPFs/Pooled/pooled_master_wo_H2Bs/five_prime_counts.txt")
cds=read.delim("~/Transcriptome_Enhancement/RPFs/Pooled/pooled_master_wo_H2Bs/CDS_counts.txt")
head(three)
head(five)
head(cds)
head(merge(three,cds,by="Id"))
head(merge(merge(three,cds,by="Id"),five,by="Id"))
mRNA_features=merge(merge(three,cds,by="Id"),five,by="Id")
a=as.data,frame(mRNA_features,row.names=1)
a=as.data.frame(mRNA_features,row.names=1)
head(mRNA_features)
write.table(mRNA_features,"~/owncloud.gwdg.de/mRNA_feature_counts.txt",sep="\t",row.names = F,quote=F)
res=read.delim("~/owncloud.gwdg.de/mRNA_feature_counts.txt")
head(res)
res$three_rime_UTR_start=NULL
colnames(res)[1]
res=read.delim("~/owncloud.gwdg.de/mRNA_feature_counts.txt",row.names = 1)
res=read.delim("~/owncloud.gwdg.de/mRNA_feature_counts.txt")
dim(mRNA_features)
unique(mRNA_features)
dim(unique(mRNA_features))
dim(unique(mRNA_features$Id))
res[!duplicated(res), ]
library(data.table)
res_dt=data.table(res)
unique(res_dt, by = "Id")
as.data.frame(unique(res_dt, by = "Id"))
dim(as.data.frame(unique(res_dt, by = "Id")))
mRNA_features=as.data.frame(unique(res_dt, by = "Id"))
write.table(mRNA_features,"~/owncloud.gwdg.de/mRNA_feature_counts.txt",sep="\t",row.names = F,quote=F)
res=read.delim("~/owncloud.gwdg.de/mRNA_feature_counts.txt",row.names = 1)
head(res)
boxplot(res)
boxplot(logs(res))
boxplot(log2(res))
boxplot(log10(res))
boxplot(log2(res))
res[,c(1,3,2)]
res=res[,c(1,3,2)]
head(res)
?boxplot
boxplot(res,col=c("green","blue","red"))
boxplot(log2(res),col=c("green","blue","red"))
dim(res)
dmel_fasta_header=read.delim("~/Drosophila_melanogaster/Dunn_et_al/dmel-CDS_fasta_headers.txt", header=F)
head(dmel_fasta_header)
colnames(dmel_fasta_header)
colnames(dmel_fasta_header)=c("Gene_names","g_ID_tr_ID")
RPF_counts=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/target_counts_Dunn.csv")
head(RPF_counts)
RPF_counts=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/target_counts_Dunn.csv",sep = ",")
head(RPF_counts)
RPF_counts$Experiment=NULL
head(RPF_counts)
colnames(RPF_counts)[1]=c("g_ID_tr_ID")
head(RPF_counts)
head(merge(dmel_fasta_header,RPF_counts,by="g_ID_tr_ID"))
dim(merge(dmel_fasta_header,RPF_counts,by="g_ID_tr_ID"))
dim(RPF_counts)
res=merge(dmel_fasta_header,RPF_counts,by="g_ID_tr_ID")
setdiff(res$g_ID_tr_ID,RPF_counts$g_ID_tr_ID)
setdiff(as.character(res$g_ID_tr_ID),as.character(RPF_counts$g_ID_tr_ID))
setdiff(res$g_ID_tr_ID,RPF_counts$g_ID_tr_ID)
dim(res[res$Reads> 500,])
sample(res$Gene_names,1000)
sort(sample(res$Gene_names,1000))
dim(sort(sample(res$Gene_names,1000)))
length(sort(sample(res$Gene_names,1000)))
length(sort(sample(res$Gene_names[res$Reads>500,],1000)))
length(sort(sample(res[res$Reads>500,],1000)))
res[res$Reads>500,]
ls()
head(res)
dim(res)
dim(res[res$Reads>500,])
dim(res[res$Reads>300,])
res_g_300=res[res$Reads>300,]
random_1000_gene_names=sort(sample(res_g_300$Gene_names,1000))
random_1000_gene_names
write.table(random_1000_gene_names,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Random_1000_genes.txt",sep="\t"row.names=F,quote=F)
write.table(random_1000_gene_names,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Random_1000_genes.txt",sep="\t",row.names=F,quote=F)
intron_len=read.delim("~/Dmel_intronseq.fasta.fai",header=F)
head(intron_len)
intron_len$V3=NULL
intron_len$V4=NULL
intron_len$V5=NULL
plot(log(intron_len))
rownames(intron_len)=intron_len$V1
intron_len$V1=NULL
head(intron_len)
plot(log2(intron_len$V2))
plot(density(log2(intron_len$V2)))
log2(70)
plot(density(log2(intron_len$V2)),main = "Drosophila_Intron_length_Distribution",xlab = "Log2_intron_len")
random_1000_gene_names_last_50_100nts_deleted=sort(sample(res_g_300$Gene_names,1000))
write.table(random_1000_gene_names_last_50_100nts_deleted,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/Random_gene_names_last_100nts_deleted.txt",sep="\t",row.names=F,quote=F)
write.table(random_1000_gene_names_last_50_100nts_deleted,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/Random_1000_gene_names_last_1000nts_deleted.txt",sep="\t",row.names=F,quote=F)
log2(200)
log2(250)
sample=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/sample.bed",header=F)
samples
sample
colnames(sample)
colnames(sample)=c("gene","start","end")
group_by(sample, gene)
library(data.table)
library(dplyr)
install.packages("dplyr")
library(dplyr)
group_by(dt, group)
dt
group_by(sample, sample$gene)
?group_by
sample=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/sample.bed",header=F)
sample
colnames(sample)=c("gene","start","end")
sample
sample=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/sample.bed",header=F)
colnames(sample)=c("gene","start","end")
sample
colnames(sample)[4]="Strand"
library(data.table)
sample_dt=as.data.table(sample)
sample_dt
DT[ , .SD[which.min(start)], by = gene]
sample_dt[ , .SD[which.min(start)], by = gene]
sample_dt[ , .SD[which.min(start) & which.max(end())], by = gene]
sample_dt[ , .SD[which.min(start) & which.max(end())], by = gene]
sample_dt[ , .SD[which.min(start) & which.max(end)], by = gene]
summarize(sample_dt,start=min(start),max=max(end))
group_by(sample)
group_by(sample,gene)
sa=group_by(sample,gene)
summarize(sa,start=min(start),max=max(end))
sample=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/Cleaned_Extensions_first_50nts_deleted_intron_length_70.bed",header=F)
colnames(sample)
colnames(sample)=c("Gene","Start","End","Counts","Chromosome", "Transcript_Start","Transcript_End","Strand")
sa=group_by(sample,Gene)
summarize(sa,start=min(start),max=max(end))
summarize(sa,start=min(Start),max=max(End))
as.data.frame(summarize(sa,start=min(Start),max=max(End)))
first_50nts=as.data.frame(summarize(sa,start=min(Start),max=max(End)))
head(first_50nts)
first_50nts=as.data.frame(summarize(sa,start=min(Start),End=max(End)))
head(first_50nts)
dmel_cds_gene_start_end=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/dmel_CDS_start_end_strand.txt",header=F)
head(dmel_cds_gene_start_end)
colnames(dmel_cds_gene_start_end)=c("Gene","Orig_Start","Orig_end","Strand")
dim(dmel_cds_gene_start_end)
dim(first_50nts)
dim(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
head(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
tail(merge(first_50nts,dmel_cds_gene_start_end,by="Gene"))
first_50nts=merge(first_50nts,dmel_cds_gene_start_end,by="Gene")
head(first_50nts)
dim(first_50nts)
first_50nts[if(first_50nts$Strand == "-" && first_50nts$End > first_50nts$Orig_end), ]
first_50nts[if(first_50nts$Strand == "-" && first_50nts$End > first_50nts$Orig_end) ]
first_50nts[(first_50nts$Strand == "-" && first_50nts$End > first_50nts$Orig_end), ]
dim(first_50nts[(first_50nts$Strand == "-" && first_50nts$End > first_50nts$Orig_end), ])
dim(first_50nts[(first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end), ])
first_50nts[(first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end), ]
first_50nts[(first_50nts$Strand == "+" & first_50nts$start <  first_50nts$Orig_Start), ]
dim(first_50nts[(first_50nts$Strand == "+" & first_50nts$start <  first_50nts$Orig_Start), ])
dim(first_50nts[(first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end), ])
dim(first_50nts[(first_50nts$Strand == "+" & first_50nts$start <  first_50nts$Orig_Start |first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end) ), ])
dim(first_50nts[(first_50nts$Strand == "+" & first_50nts$start <  first_50nts$Orig_Start |first_50nts$Strand == "-" & first_50nts$End > first_50nts$Orig_end ), ])
write.table(first_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/Extensions_vs_original_start_end.txt",sep="\t",quote=F,row.names(F))
write.table(first_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/Extensions_vs_original_start_end.txt",sep="\t",quote=F,row.names=F)
first_100=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
colnames(first_100)=c("Gene","Start","End","Counts","Chromosome", "Transcript_Start","Transcript_End","Strand")
first_100=group_by(first_100,Gene)
dim(as.data.frame(summarize(first_100,Extn_start=min(Start),Extn_end=max(End))))
as.data.frame(summarize(first_100,Extn_start=min(Start),Extn_end=max(End)))
first_100=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
colnames(first_100)=c("Gene","Start","End","Counts","Chromosome", "Transcript_Start","Transcript_End","Strand")
dim(first_100)
first_100_group=group_by(sample,Gene)
first_100_group
head(first_100)
colnames(first_100)=c("Chromosome","Start","End","Counts","Gene", "Transcript_Start","Transcript_End","Strand")
first_100_group=group_by(first_100,Gene)
first_100_group
dim(as.data.frame(summarize(first_100_group,start=min(Start),Extn_end=max(End))))
as.data.frame(summarize(first_100_group,start=min(Start),Extn_end=max(End)))
head(first_100)
first_100=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
head(first_100)
colnames(first_100)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
head(first_100)
rm(first_100_group)
as.data.frame(summarize(first_100,start=min(Extn_Start),Extn_end=max(Extn_End)))
first_100_group=group_by(first_100,Gene)
as.data.frame(summarize(first_100_group,start=min(Extn_Start),Extn_end=max(Extn_End)))
first_100nts=as.data.frame(summarize(first_100_group,start=min(Extn_Start),Extn_end=max(Extn_End)))
dim(first_100)
dim(first_100nts)
first_100nts=as.data.frame(summarize(first_100_group,start=min(Extn_Start),Extn_end=max(Extn_End),Genome_strand="Genome_strand"))
first_100nts
first_100nts=as.data.frame(summarize(first_100_group,start=min(Extn_Start),Extn_end=max(Extn_End),Genome_strand=Genome_strand))
first_100nts=as.data.frame(summarize(first_100_group,start=min(Extn_Start),Extn_end=max(Extn_End))
)
first_100nts
dim(merge(first_100nts,dmel_cds_gene_start_end,by="Gene"))
first_100nts=merge(first_100nts,dmel_cds_gene_start_end,by="Gene")
head(first_100nts)
colnames(first_100nts)[2]
colnames(first_100nts)[2]="Extn_start"
dim(first_100nts[(first_100nts$Strand == "+" & first_100nts$Extn_start <  first_100nts$Orig_Start |first_100nts$Strand == "-" & first_100nts$Extn_end > first_50nts$Orig_end ), ])
dim(first_100nts[(first_100nts$Strand == "+" & first_100nts$Extn_start <  first_100nts$Orig_Start |first_100nts$Strand == "-" & first_100nts$Extn_end > first_100nts$Orig_end ), ])
Extn_vs_orig_first_100nts=first_100nts[(first_100nts$Strand == "+" & first_100nts$Extn_start <  first_100nts$Orig_Start |first_100nts$Strand == "-" & first_100nts$Extn_end > first_100nts$Orig_end ), ]
write.table(Extn_vs_orig_first_100nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/Extn_vs_orig_first_100bts.txt",sep="\t",quote=F,row.names = F)
dim(first_100)
first_100
head(first_100)
dim(first_100_group)
first_100_group
dim(first_100nts)
ls()
dim(first_50nts)
last_50=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
head(last_50)
colnames(first_100)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
colnames(last_50)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
last_50_group=group_by(last_50,Gene)
head(last_50_group)
dim(as.data.frame(summarize(last_50_group,Extn_start=min(Extn_start),Extn_end=max(Extn_End))))
dim(as.data.frame(summarize(last_50_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End))))
last_50=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/Cleaned_Extensions_last_50nts_deleted_intron_length_70.bed",header=F)
colnames(last_50)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
last_50_group=group_by(last_50,Gene)
dim(as.data.frame(summarize(last_50_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End))))
last_50nts=as.data.frame(summarize(last_50_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End))
)
dim(first_50nts[(first_50nts$Strand == "-" & first_50nts$start <  first_50nts$Orig_Start |first_50nts$Strand == "+" & first_50nts$End > first_50nts$Orig_end ), ])
last_50nts=merge(last_50nts,dmel_cds_gene_start_end,by="Gene")
head(last_50nts)
dim(last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
dim(last_50nts[(last_50nts$Strand == "-" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "+" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
last_50nt_extn_vs_orig=last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ]
write.table(last_50nt_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/last_50nt_extn_vs_orig.txt",sep="\t",quote=F,row.names = F)
last_100=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/Cleaned_Extensions_last_100nts_deleted_intron_length_70.bed",header=F)
colnames(last_50)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
colnames(last_100)=c("Chromosome","Extn_Start","Extn_End","Counts","Genome_strand", "Gene","Transcript_Start","Transcript_End","Strand")
last_100_group=group_by(last_100,Gene)
dim(as.data.frame(summarize(last_100_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End))))
last_100nts=as.data.frame(summarize(last_100_group,Extn_start=min(Extn_Start),Extn_end=max(Extn_End)))
last_100nts=merge(last_100nts,dmel_cds_gene_start_end,by="Gene")
head(last_100nts)
dim(last_100nts[(last_100nts$Strand == "-" & last_100nts$Extn_start <  last_100nts$Orig_Start | last_100nts$Strand == "+" & last_100nts$Extn_end > last_100nts$Orig_end ), ])
last_100nts-extn_vs_orig=last_100nts[(last_100nts$Strand == "-" & last_100nts$Extn_start <  last_100nts$Orig_Start | last_100nts$Strand == "+" & last_100nts$Extn_end > last_100nts$Orig_end ), ]
last_100nts_extn_vs_orig=last_100nts[(last_100nts$Strand == "-" & last_100nts$Extn_start <  last_100nts$Orig_Start | last_100nts$Strand == "+" & last_100nts$Extn_end > last_100nts$Orig_end ), ]
write.table(last_100nts_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/Last_100nts_extn_vs_orig.txt",sep="\t",quote=F,row.names = F)
ls()
dim(last_100nts)
dim(last_50nts)
dim(first_50nts)
dim(first_100nts)
write.table(last_100nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/LAST_100nts_deleted_Extensions_vs_original.txt",sep="\t",quote=F,row.names = F)
write.table(last_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/LAST_50nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
write.table(first_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted//FIRST_50nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
write.table(first_100nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/FIRST_100nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
write.table(first_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/FIRST_50nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
write.table(last_50nt_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/last_50nts_deleted_extensions_greater_than_Orig.txt",sep="\t",quote=F,row.names = F)
dim(last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
head(last_50nts)
dim(last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
dim(last_50nts[(last_50nts$Strand == "-" & last_50nts$Extn_start <  last_50nts$Orig_start | last_50nts$Strand == "+" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
dim(last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
dim(last_50nts[(last_50nts$Strand == "+" & last_50nts$Extn_start <  last_50nts$Orig_start | last_50nts$Strand == "-" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
head(first_50nts)
first_100_chromosome=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/Cleaned_Extensions_first_100nts_deleted_intron_length_70.bed",header=F)
head(first_100_chromosome)
first_100_chromosome$V2=NULL
first_100_chromosome$V3=NULL
first_100_chromosome$V4=NULL
first_100_chromosome$V5=NULL
first_100_chromosome$V7=NULL
first_100_chromosome$V8=NULL
first_100_chromosome$V9=NULL
head(first_100_chromosome)
dim(merge(first_100_chromosome,first_50nts,by.x="V2",by.y="Gene"))
colnames(first_100_chromosome)=c("Chromosome","Gene")
dim(merge(first_100_chromosome,first_50nts,by="Gene"))
head(merge(first_100_chromosome,first_50nts,by="Gene"))
head(merge(first_50nts,first_100_chromosome,by="Gene"))
dim(merge(first_50nts,first_100_chromosome,by="Gene"))
dim(first_50nts)
dim(unique[,first_100_chromosome$Chromosome])
head(first_100_chromosome)
head(unique(first_100_chromosome$Gene))
head(unique [,first_100_chromosome$Gene])
head(unique[first_100_chromosome$Gene,])
first_100_chromosome[!duplicated(first_100_chromosome)]
first_100_chromosome[!duplicated(first_100_chromosome),]
first_100_chromosome=first_100_chromosome[!duplicated(first_100_chromosome),]
head(first_100_chromosome)
dim(first_100_chromosome)
first_50_chromosome=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/Cleaned_Extensions_first_50nts_deleted_intron_length_70.bed",header = F)
colnames(first_50_chromosome)
dim(first_50_chromosome)
first_500_chromosome=first_50_chromosome[!duplicated(first_50_chromosome),]
dim(first_50_chromosome)
dim(first_500_chromosome)
head(first_50_chromosome)
rm(first_500_chromosome)
first_50_chromosome$V2=NULL
first_50_chromosome$V3=NULL
first_50_chromosome$V4=NULL
first_50_chromosome$V6=NULL
first_50_chromosome$V7=NULL
first_50_chromosome$V8=NULL
head(first_50_chromosome)
first_50_chromosome=first_50_chromosome[!duplicated(first_50_chromosome),]
dim(first_50_chromosome)
colnames(first_50_chromosome)
head(first_50_chromosome)
colnames(first_50_chromosome)=c("Gene","Chromosome")
dim(merge(first_50nts,first_50_chromosome,by="Gene"))
head(merge(first_50nts,first_50_chromosome,by="Gene"))
first_50nts=merge(first_50nts,first_50_chromosome,by="Gene")
write.table(first_50nts,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_50nts_deleted/FIRST_50nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
ls()
head(first_100_chromosome)
head(first_100)
dim(first_100)
dim(first_100_chromosome)
head(first_100_chromosome)
first_100_deleted_extn_vs_orig=read.delim("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/FIRST_100nts_deleted_extensions_vs_Orig.txt",header = T)
head(first_100_deleted_extn_vs_orig)
dim(first_100_deleted_extn_vs_orig)
dim(merge(first_100_deleted_extn_vs_orig,first_100_chromosome,by="Gene"))
head(merge(first_100_deleted_extn_vs_orig,first_100_chromosome,by="Gene"))
first_100_deleted_extn_vs_orig=merge(first_100_deleted_extn_vs_orig,first_100_chromosome,by="Gene")
head(first_100_deleted_extn_vs_orig)
dim(first_100_deleted_extn_vs_orig)
write.table(first_100_deleted_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/First_100nts_deleted/FIRST_100nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
dim(last_50nt_extn_vs_orig)
last_50nt_extn_vs_orig
last_50nt_extn_vs_orig=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/LAST_50nts_deleted_extensions_vs_Orig.txt")
head(last_100nts_extn_vs_orig)
dim(last_100nts_extn_vs_orig)
ls()
last_50_chromosome=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/Cleaned_Extensions_last_50nts_deleted_intron_length_70.bed",header=F)
head(last_50_chromosome)
head(last_50_group)
last_50_chromosome$V2=NULL
last_50_chromosome$V3=NULL
last_50_chromosome$V4=NULL
last_50_chromosome$V5=NULL
last_50_chromosome$V7=NULL
last_50_chromosome$V8=NULL
last_50_chromosome$V9=NULL
head(last_50nt_extn_vs_orig)
colnames(last_50_chromosome)
head(last_50_chromosome)
colnames(last_50_chromosome)=c("Chromosome","Gene")
head(last_50_chromosome)
last_50_chromosome[!duplicated(last_50_chromosome)]
last_50_chromosome[!duplicated(last_50_chromosome),]
dim(last_50_chromosome[!duplicated(last_50_chromosome),])
last_50_chromosome=last_50_chromosome[!duplicated(last_50_chromosome),]
head(last_50_chromosome)
head(last_50nt_extn_vs_orig)
dim(merge(last_50nt_extn_vs_orig,last_50_chromosome,by="Gene"))
dim(last_50nt_extn_vs_orig)
dim(last_50nts[(last_50nts$Strand == "-" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "+" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
head(last_50nts[(last_50nts$Strand == "-" & last_50nts$Extn_start <  last_50nts$Orig_Start | last_50nts$Strand == "+" & last_50nts$Extn_end > last_50nts$Orig_end ), ])
dim(merge(last_50nt_extn_vs_orig,last_50_chromosome,by="Gene"))
head(merge(last_50nt_extn_vs_orig,last_50_chromosome,by="Gene"))
last_50nt_extn_vs_orig=merge(last_50nt_extn_vs_orig,last_50_chromosome,by="Gene")
write.table(last_50nt_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_50nts_deleted/LAST_50nts_deleted_extensions_vs_Orig.txt",sep="\t",quote=F,row.names = F)
last_100nt_extn_vs_orig=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/LAST_100nts_deleted_Extensions_vs_original.txt")
head(last_100nts_extn_vs_orig)
dim(last_100nt_extn_vs_orig)
dim(last_100nts_extn_vs_orig)
head(last_100nts_extn_vs_orig)
last_100nts_chromosome=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/Cleaned_Extensions_last_100nts_deleted_intron_length_70.bed")
head(last_100nts_chromosome)
last_100nts_chromosome=read.delim("/home/tsekara/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/Cleaned_Extensions_last_100nts_deleted_intron_length_70.bed",header=F)
head(last_100nts_chromosome)
last_100nts_chromosome$V2=NULL
last_100nts_chromosome$V3=NULL
last_100nts_chromosome$V4=NULL
last_100nts_chromosome$V5=NULL
last_100nts_chromosome$V7=NULL
last_100nts_chromosome$V8=NULL
last_100nts_chromosome$V9=NULL
colnames(last_100nts_chromosome)
head(last_100nts_chromosome)
colnames(last_100nts_chromosome)=c("Chromosome","Gene")
last_100nts_chromosome=last_100nts_chromosome[!duplicated(last_100nts_chromosome),]
head(last_100nts_chromosome)
dim(last_100nts_chromosome)
head(merge(last_100nt_extn_vs_orig,last_100nts_chromosome,by="Gene"))
dim(merge(last_100nt_extn_vs_orig,last_100nts_chromosome,by="Gene"))
last_100nt_extn_vs_orig=merge(last_100nt_extn_vs_orig,last_100nts_chromosome,by="Gene")
dim(last_100nt_extn_vs_orig)
write.table(last_100nt_extn_vs_orig,"~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Last_100nts_deleted/LAST_100nts_deleted_Extensions_vs_original.txt",sep="\t",quote=F,row.names = F)
getwd()
save.image("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Dmel.RData")
savehistory("~/Drosophila_melanogaster/Dunn_et_al/trimmed/single_replicates/Dmel.R")
