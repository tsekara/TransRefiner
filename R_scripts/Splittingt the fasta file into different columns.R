###### Splittingt the fasta file into different columns ########################

lines <- readLines('/home/tsekara/Transcriptomes/trinity5_frameDPcorrected.fa')
indx <- grepl('>', lines)
Sequence <- tapply(seq_along(indx),cumsum(indx), FUN=function(x) 
  paste(lines[tail(x,-1)], collapse=""))
d1 <- data.frame(names=lines[indx], Sequence, stringsAsFactors=FALSE)
head(d1,2)



##################### Splitting the fasta file into two columns using Seqinr package ############################

sep=data.frame(Transcript_ID=names(tr5_fasta), Seqs=unlist(getSequence(tr5_fasta, as.string=T)))