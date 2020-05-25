library(rtracklayer)
library(GenomicFeatures)
myaln=readGAlignments(BamFile("C:/Users/tsekaran/Documents/Sorted_GFP_RPF_linker_free_trimmed_2.fastq.bam"))
mytxdb=makeTxDbFromGFF("C:/Users/tsekaran/Documents/smed_dd_v6.fasta.transdecoder.gff3",organism = "Schmidtea mediterranea", dataSource = "Transdecoder")
matchLenDistr <- histMatchLength(myaln, 0)
matchLenDistr[[2]]
seqlevels(myaln)
seqlevels(mytxdb)
alnGRanges <- readsToReadStart(myaln)
oneBinRanges <- aroundPromoter(mytxdb, alnGRanges, flankSize = 20)
oneBinRanges
readStartCov(alnGRanges,oneBinRanges,matchSize =c(29:31),fixedInterval =c(-20,20),renameChr = "aroundTSS",charPerc = "perc")
plotSummarizedCov(readStartCov(alnGRanges,oneBinRanges,matchSize =c(29:31),fixedInterval =c(-20,20),renameChr = "aroundTSS",charPerc = "perc"))
