setwd("~/Transcriptome_Enhancement/Intron_length/")
dat=read.delim("5prime_UTR_intron_length.txt")
head(dat)
dat$start=NULL
dat$end=NULL
head(dat)
plot(log2(dat$width))
plot(density(log2(dat$width)))
plot(density((dat$width)))
plot(density(log2(dat$width)))
plot(density(log2(dat$width)),xlim=c(0,1,2,3,4,5,6,7,8,9,10))
library(ggplot2)
ggplot(dat, aes(width)) + geom_density()
ggplot(dat, aes(log2(width))) + geom_density()
ggplot(dat, aes(log2(width))) + geom_density(adjust=1/10)
ggplot(dat, aes(log2(width))) + geom_density(adjust=10)
library(GenomicFeatures)
?makeTxDbFromGFF
mytxdb=makeTxDbFromGFF("gmap_smed_dd_v6_vs_Asxl_genome_no_chimeric_alignments.gff3")
?fiveUTRsByTranscript
head(intronsByTranscript(mytxdb,use.names=T))
head(intronsByTranscript(mytxdb,use.names=F
))
head(as.data.frame(intronsByTranscript(mytxdb,use.names=T)))
introns_within_transcript=as.data.frame(intronsByTranscript(mytxdb,use.names=T))
head(introns_within_transcript)
introns_within_transcript$group=NULL
introns_within_transcript$seqnames=NULL
introns_within_transcript$start=NULL
introns_within_transcript$end=NULL
introns_within_transcript$strand=NULL
plot(density(log2(introns_within_transcript$width)))
?plot
plot(density(log2(introns_within_transcript$width)),col=red)
plot(density(log2(introns_within_transcript$width)),col="red")
plot(density(log2(introns_within_transcript$width)),col="red",main="Intron_length_distribution_within_transcripts")
head(aggregate(introns_within_transcript$width, data = introns_within_transcript, median))
head(aggregate(introns_within_transcript, data = introns_within_transcript$group_name, FUN=median))
head(aggregate(introns_within_transcript, by = introns_within_transcript$group_name, FUN=median))
head(aggregate(introns_within_transcript, by=as-list(introns_within_transcript$group_name), FUN=median))
head(aggregate(introns_within_transcript, by=as.list(introns_within_transcript$group_name), FUN=median))
aggregate(introns_within_transcript$width,by=list(introns_within_transcript$group_name),FUN=median)
head(aggregate(introns_within_transcript$width,by=list(introns_within_transcript$group_name),FUN=median))
introns__median_lenght_within_transcript=aggregate(introns_within_transcript$width,by=list(introns_within_transcript$group_name),FUN=median)
plot(density(log2(introns__median_lenght_within_transcript$x)))
three_prime_UTR=threeUTRsByTranscript(mytxdb,use.names=T)
three_prime_UTR
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
introns_median_length_within_transcript=aggregate(three_prime_UTR_intron$length,by=list(three_prime_UTR_intron$group_name),FUN=median)
head(introns__median_lenght_within_transcript)
plot(density(log2(introns_median_length_within_transcript$x)))
plot(density(log2(introns_median_length_within_transcript$x)),main="Three prime UTR intron length distribution")
plot(density(log2(three_prime_UTR_intron$length)),main="Intron length distribution in three prime UTR")
head(as.data.frame(cdsBy(mytxdb,use.names=T)))
cds=as.data.frame(cdsBy(mytxdb,use.names=T))
cds$group=NULL
cds$seqnames=NULL
cds$strand=NULL
cds$cds_id=NULL
cds$cds_name=NULL
cds$exon_rank=NULL
cds_intron=abs(setDT(cds)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name])
cds_intron=setDT(cds)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
cds_intron$length=abs(cds_intron$end-cds_intron$start)
head(cds)
head(cds_intron)
plot(density(log2(cds_intron$length)))
plot(density(log2(cds_intron$length)),main="Intron length distribution in CDS")
cds_introns_median=aggregate(cds_intron$length,by=list(cds_intron$group_name),FUN=median)
head(cds_introns_median)
plot(density(log2(cds_introns_median$x)))
plot(density(log2(dat$width)))
plot(density(log2(dat$width)),main="Intron length distribution in 5 prine UTR")
savehistory("Intron_length_distribution_analysis.R")
