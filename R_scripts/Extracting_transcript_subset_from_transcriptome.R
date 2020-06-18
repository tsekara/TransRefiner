# Extracting transcript subset from transcriptome 

library(seqinr)
ids=as.character(read.delim("path/to/ids/file.txt"))
dd=read.fasta("path/to/transcriptome.fasta",seqtype="DNA",as.string=T)
fasta_seq = dd[names(dd) %in% ids]
write.fasta(fasta_seq, names(fasta_seq), file=paste(dir,"/",name,sep=""))