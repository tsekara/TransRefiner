mytxdb=makeTxDbFromGFF("../GMAP_vs_dd_v6_pcf_idn_0.7_cov_0.7.gff3")
introns_within_transcript=as.data.frame(intronsByTranscript(mytxdb,use.names=T))
head(introns_within_transcript)
introns_within_transcript$group=NULL
introns_within_transcript$seqnames=NULL
introns_within_transcript$start=NULL
introns_within_transcript$end=NULL
introns_within_transcript$strand=NULL
head(introns_within_transcript)
dim(introns_within_transcript)
introns_length_within_transcript_median_aggregate=aggregate(introns_within_transcript$width,by=list(introns_within_transcript$group_name),FUN=median)
dim(introns_length_within_transcript_median_aggregate)
head(introns_length_within_transcript_median_aggregate)
colnames(introns_length_within_transcript_median_aggregate)=c("Transcript","Median-Intron-length")
head(introns_length_within_transcript_median_aggregate)
plot(density(log2(introns_length_within_transcript_median_aggregate$`Median-Intron-length`)))
plot(density(introns_length_within_transcript_median_aggregate$`Median-Intron-length`))
plot(density(introns_length_within_transcript_median_aggregate$`Median-Intron-length`),main="Intron_length_distribution_within_transcripts")
plot(density(log2(introns_length_within_transcript_median_aggregate$`Median-Intron-length`)),main="Intron_length_distribution_within_transcripts")
plot(density(log2(introns_length_within_transcript_median_aggregate$`Median-Intron-length`)),main="Intron_length_distribution_entire_transcripts")
plot(density(log2(introns_length_within_transcript_median_aggregate$`Median-Intron-length`)),main="Intron length distribution entire transcripts")
three_prime_UTR=as.data.frame(threeUTRsByTranscript(mytxdb,use.names=T))
three_prime_UTR$group=NULL
three_prime_UTR$strand=NULL
three_prime_UTR$exon_id=NULL
three_prime_UTR$exon_name=NULL
three_prime_UTR$exon_rank=NULL
head(three_prime_UTR)
three_prime_UTR$seqnames=NULL
three_prime_UTR$start=NULL
three_prime_UTR$end=NULL
head(three_prime_UTR)
library(data.table)
dim(three_prime_UTR)
three_prime_UTR_intron=setDT(three_prime_UTR)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
head(three_prime_UTR_intron)
three_prime_UTR_intron$length=abs(three_prime_UTR_intron$start-three_prime_UTR_intron$end)
head(three_prime_UTR_intron)
ls()
rm(three_prime_UTR)
three_prime_UTR=threeUTRsByTranscript(mytxdb,use.names=T)
three_prime_UTR=as.data.frame(threeUTRsByTranscript(mytxdb,use.names=T))
three_prime_UTR$group=NULL
three_prime_UTR$strand=NULL
three_prime_UTR$exon_id=NULL
three_prime_UTR$exon_name=NULL
three_prime_UTR$exon_rank=NULL
head(three_prime_UTR)
three_prime_UTR$seqnames=NULL
head(three_prime_UTR)
dim(three_prime_UTR)
three_prime_UTR$width=NULL
library(data.table)
head(three_prime_UTR)
three_prime_UTR_intron=setDT(three_prime_UTR)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
head(three_prime_UTR_intron)
three_prime_UTR_intron$length=abs(three_prime_UTR_intron$start-three_prime_UTR_intron$end)
head(three_prime_UTR_intron)
plot(density(log2(three_prime_UTR_intron$length)))
plot(density(log2(three_prime_UTR_intron$length)),main="Intron length distribution 3' UTR")
plot(density(log2(three_prime_UTR_intron$length)),main="Intron length distribution in 3' UTR")
dim(aggregate(three_prime_UTR_intron$length,by=list(three_prime_UTR_intron$group_name),FUN=median))
dim(three_prime_UTR)
head(three_prime_UTR)
dim(three_prime_UTR_intron)
head(three_prime_UTR_intron)
length(unique(three_prime_UTR_intron$group_name))
three_prime_UTR_intron=aggregate(three_prime_UTR_intron$length,by=list(three_prime_UTR_intron$group_name),FUN=median)
dim(three_prime_UTR_intron)
head(three_prime_UTR_intron)
colnames(three_prime_UTR_intron)=c("Transcript","Intron_length_in_three_prime_UTR")
plot(density(log2(three_prime_UTR_intron$Intron_length_in_three_prime_UTR)),main="Intron length distribution in three prime UTR")
cds=as.data.frame(cdsBy(mytxdb,use.names=T))
cds$group=NULL
cds$seqnames=NULL
cds$strand=NULL
cds$cds_id=NULL
cds$cds_name=NULL
cds$exon_rank=NULL
cds_intron=setDT(cds)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
cds_intron$length=abs(cds_intron$end-cds_intron$start)
cds_introns_median=aggregate(cds_intron$length,by=list(cds_intron$group_name),FUN=median)
head(cds_intron)
head(cds_introns_median)
dim(cds_introns_median)
colnames(cds_introns_median)
colnames(cds_introns_median)=c("Transcript","Intron_length_CDS")
plot(density(log2(cds_introns_median$Intron_length_CDS)),main="Intron length distribution in CDS")
five_prime_UTR=fiveUTRsByTranscript(mytxdb,use.names=T)
five_prime_UTR=as.data.frame(fiveUTRsByTranscript(mytxdb,use.names=T))
five_prime_UTR$group=NULL
five_prime_UTR$strand=NULL
five_prime_UTR$exon_id=NULL
five_prime_UTR$exon_name=NULL
five_prime_UTR$exon_rank=NULL
five_prime_UTR$seqnames=NULL
five_prime_UTR$width=NULL
five_prime_UTR_intron=setDT(five_prime_UTR)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
head(five_prime_UTR_intron)
five_prime_UTR_intron$length=abs(five_prime_UTR_intron$start-five_prime_UTR_intron$end)
five_prime_UTR_intron=aggregate(five_prime_UTR_intron$length,by=list(five_prime_UTR_intron$group_name),FUN=median)
plot(density(log2(five_prime_UTR_intron$length)),main="Intron length distribution in five prime UTR")
plot(density(log2(five_prime_UTR_intron$x)),main="Intron length distribution in five prime UTR")
a$length=abs(five_prime_UTR_intron$start-five_prime_UTR_intron$end)
dim(five_prime_UTR_intron)
a=five_prime_UTR_intron
three_prime_UTR_intron$length=abs(three_prime_UTR_intron$start-three_prime_UTR_intron$end)
save.image("Dec2017/Intron_length_distribution.RData")
savehistory("Dec2017/Intron_length_distribution.R")
