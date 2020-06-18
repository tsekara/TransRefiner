
################################## Read Length Dstrubution #############################################
aln=readGAlignments(BamFile("http://genomique.info/data/public/RiboProfiling/ctrl.bam"))
aln=ctrlGAlignments
matchLenDistr=histMatchLength(aln,0)
matchLenDistr[[2]]

################################### read start Coverage plot around the TSS ##########################################
alnGRanges=readsToReadStart(aln)
#txdb object with annotations
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
oneBinRanges<-aroundPromoter(txdb, alnGRanges,percBestExpressed=0.001)
#the coverage in the TSS flanking region for the reads with match sizes 29:31
listPromoterCov<-readStartCov(alnGRanges,oneBinRanges,matchSize=c(29:31),fixedInterval=c(-20,20),renameChr="aroundTSS",charPerc="perc")
plotSummarizedCov(listPromoterCov)

########################### Count Reads on Features ##########################################

alnGRanges=  alnGRanges[which(!is.na(match(alnGRanges$score,30:33)))]
cds=cdsBy(mytxdb,by="tx",use.names= TRUE)
exonGRanges=exonsBy(mytxdb,by="tx",use.names= TRUE)
cdsPosTransc=orfRelativePos(cds, exonGRanges)
countsDataCtrl1=countShiftReads(exonGRanges=exonGRanges[names(cdsPosTransc)],cdsPosTransc=cdsPosTransc,alnGRanges=alnGRanges,shiftValue=-14)
listCountsPlots=countsPlot(list(countsDataCtrl1[[1]]),grep("_counts$",colnames(countsDataCtrl1[[1]])),1)

######################## Count Reads on Codons #########################################

listReadsCodon=countsDataCtrl1[[2]]
cds=cdsBy(mytxdb,use.names=TRUE)
orfCoord=cds[names(cds) %in% names(listReadsCodon)]
library(BSgenome.Hsapiens.UCSC.hg19)
#chromosome names should correspond between the BAM,
#the annotations, and the genome
genomeSeq=BSgenome.Hsapiens.UCSC.hg19
#codon frequency, coverage, and annotation
codonData=codonInfo(listReadsCodon, genomeSeq, orfCoord)


codonUsage=codonData[[1]]
codonCovMatrix=codonData[[2]]
#keep only genes with a minimum number of reads
nbrReadsGene=apply(codonCovMatrix,1, sum)
ixExpGenes=which(nbrReadsGene>=50)
codonCovMatrix=codonCovMatrix[ixExpGenes, ]
#get the PCA on the codon coverage
codonCovMatrixTransp=t(codonCovMatrix)
rownames(codonCovMatrixTransp)=colnames(codonCovMatrix)
colnames(codonCovMatrixTransp)=rownames(codonCovMatrix)
listPCACodonCoverage=codonPCA(codonCovMatrixTransp,"codonCoverage")
listPCACodonCoverage[[2]]








