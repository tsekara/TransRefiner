dresden=function(ids,file_name,directory)
{
  library(seqinr)
  ids=as.character(ids)
  name=file_name
  dir=directory
  dd=read.fasta("/home/tsekara/Transcriptomes/smed_dd_v6.fasta",seqtype ="DNA",as.string = T)
  fasta_seq=dd[names(dd) %in% ids]
  write.fasta(fasta_seq,file=paste(dir,"/",name,sep=""))
}


library(seqinr)
ids=as.character(read.delim("path/to/ids/file.txt"))
dd=read.fasta("path/to/transcriptome.fasta",seqtype="DNA",as.string=T)
fasta_seq = dd[names(dd) %in% ids]
write.fasta(fasta_seq, names(fasta_seq), file=paste(dir,"/",name,sep=""))