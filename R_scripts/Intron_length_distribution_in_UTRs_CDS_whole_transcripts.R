library(GenomicFeatures)
library(data.table)

mytxdb=makeTxDbFromGFF("gmap_smed_dd_v6_vs_Asxl_genome_no_chimeric_alignments.gff3")

####################### Whole transcripts ################################################
introns_within_transcript=as.data.frame(intronsByTranscript(mytxdb,use.names=T))
introns_within_transcript$group=NULL
introns_within_transcript$seqnames=NULL
introns_within_transcript$start=NULL
introns_within_transcript$end=NULL
introns_within_transcript$strand=NULL
plot(density(log2(introns_within_transcript$width)),col=red)
introns_median_length_within_transcript=aggregate(introns_within_transcript$width,by=list(introns_within_transcript$group_name),FUN=median)
plot(density(log2(introns_median_length_within_transcript)))
########################### Intron length in 3'UTR #######################################
three_prime_UTR=threeUTRsByTranscript(mytxdb,use.names=T)
three_prime_UTR=as.data.frame(threeUTRsByTranscript(mytxdb,use.names=T))
three_prime_UTR$group=NULL
three_prime_UTR$strand=NULL
three_prime_UTR$exon_id=NULL
three_prime_UTR$exon_name=NULL
three_prime_UTR$exon_rank=NULL
three_prime_UTR$seqnames=NULL
three_prime_UTR$width=NULL
three_prime_UTR_intron=setDT(three_prime_UTR)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
three_prime_UTR_intron$length=abs(three_prime_UTR_intron$start-three_prime_UTR_intron$end)
plot(density(log2(three_prime_UTR_intron$length)))
########################################### Median #######################################
introns_median_length_three_prime_UTR=aggregate(three_prime_UTR_intron$length,by=list(three_prime_UTR_intron$group_name),FUN=median)
plot(density(log2(three_prime_UTR_intron$length)),main="Intron length distribution in three prime UTR")

########################### Intron length in 5'UTR #######################################
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
five_prime_UTR_intron$length=abs(five_prime_UTR_intron$end-five_prime_UTR_intron$start)
plot(density(log2(five_prime_UTR_intron$length)))
########################################### Median ########################################
introns_median_length_five_prime_UTR=aggregate(five_prime_UTR_intron$length,by=list(five_prime_UTR_intron$group_name),FUN=median)
plot(density(log2(five_prime_UTR_intron$length)),main="Intron length distribution in five prime UTR")


################################# Intron length in CDS ####################################

cds=as.data.frame(cdsBy(mytxdb,use.names=T))
cds$group=NULL
cds$seqnames=NULL
cds$strand=NULL
cds$cds_id=NULL
cds$cds_name=NULL
cds$exon_rank=NULL
cds_intron=setDT(cds)[,.(start = (start+end - start +1)[-.N],end = (end +shift(start, type='lead')-end-1)[-.N] ) , by = group_name]
cds_intron$length=abs(cds_intron$end-cds_intron$start)
plot(density(log2(cds_intron$length)),main="Intron length distribution in CDS")

################################ median ###################################################

cds_introns_median=aggregate(cds_intron$length,by=list(cds_intron$group_name),FUN=median)
plot(density(log2(cds_introns_median$x)))







