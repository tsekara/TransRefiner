genomecoverage=read.delim("/home/tsekara/owncloud.gwdg.de/Xtender/Start_site_RPF_reads/Original/Complete_transcripts/Complete_transcripts_80nts/GenomeCoverage_counts.bed")
colnames(genomecoverage)=c("Target","Position","Counts")
sorted_genomecoverage=genomecoverage[order(genomecoverage$Position),]
sorted_genomecoverage$Target=NULL
median_counts=aggregate(sorted_genomecoverage$Counts ~ sorted_genomecoverage$Position, data=sorted_genomecoverage, FUN=median)
write.table(median_counts,"median_counts.txt",sep="\t",quote=F)
