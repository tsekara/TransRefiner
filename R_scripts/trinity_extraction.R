trinity=function(ids,file_name,directory)
{
  library(seqinr)
  ids=as.character(ids)
  name=file_name
  dir=directory
  dd=read.fasta("/home/tsekara/Transcriptomes/trinity5.fasta",seqtype ="DNA",as.string = T)
  fasta_seq=unlist(dd[names(dd) %in% ids])
  write.table(fasta_seq,file=paste(dir,"/",name,sep=""),sep="\t")
}