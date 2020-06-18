######### If the file names do not have ".txt" at the end, do append it before giving it as a parameter in custom-codon_periodicity function below ############# 

file_names=paste(file_names,".txt",sep="")

########### riboDat file is genrated by providing the sorted bam file in read #######################################

custom_codon_periodicity=function(file_names,riboDat,Frame_counts_path,Plots_path)
{
  library("parallel")
  library("foreach")
  library("doParallel")
  library("riboSeqR")
  library("Rsamtools")
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl, cores = detectCores() - 1)
  
  foreach(i=1:length(file_names),.packages = c("riboSeqR","Rsamtools"),.combine=rbind) %dopar%
try({
    cds=read.delim(file_names[i],as.is=T,header=F)
    cds=GRanges(cds$V1,IRanges(start=cds$V2,end=cds$V3))
    cds$frame=(start(cds)-1)%%3
    fCs=frameCounting(riboDat,cds,lengths = 27:32)
    fS=readingFrame(rC=fCs,lengths = 27:32)
    write.table(fS,file=paste(Frame_counts_path,file_names[i],sep=""),sep="\t",quote=F,row.names=F)
    pdf(file=paste(Plots_path,gsub("txt","pdf",file_names[i]),sep=""))
    plotFS(fS)
    dev.off()
  })
}
  