
##################################### GFP_1 ##############################
file_names=list.files()

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

#for(i in 1:length(file_names_3))
foreach(i=1:length(file_names),.packages = c("riboSeqR"),.combine=rbind) %dopar%
{
  cds=read.delim(file_names[i],as.is=T,header=F)
  
  cds=GRanges(cds$V1,IRanges(start=cds$V2,end=cds$V3))
  cds$frame=(start(cds)-1)%%3
  fCs=frameCounting(riboDat,cds,lengths = 27:32)
  fS=readingFrame(rC=fCs,lengths = 27:32)
  write.table(fS,file=paste("/home/tsekara/Novel_ORFs_Prediction/Codon_periodicity/RiboseqR/Frame_counting_matrix/",file_names[i],sep=""),sep="\t",quote=F,row.names=F)
  pdf(file=paste("/home/tsekara/Novel_ORFs_Prediction/Codon_periodicity//RiboseqR/Plots/",gsub("txt","pdf",file_names[i]),sep=""))
  plotFS(fS)
  dev.off()
}

##################################### GFP_2 ##############################
file_names=list.files()
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

#for(i in 1:length(gfp2))
foreach(i=1:length(gfp2),.packages = c("riboSeqR"),.combine=rbind) %dopar%
{
  cds=read.delim(gfp2[i],as.is=T,header=F)
  
  cds=GRanges(cds$V1,IRanges(start=cds$V2,end=cds$V3))
  cds$frame=(start(cds)-1)%%3
  fCs=frameCounting(riboDat_gfp_2,cds,lengths = 27:32)
  fS=readingFrame(rC=fCs,lengths = 27:32)
  write.table(fS,file=paste("/home/tsekara/Novel_ORFs_Prediction/Codon_periodicity/RiboseqR/GFP_2_Frame_counting_matrix/",gfp2[i],sep=""),sep="\t",quote=F,row.names=F)
  pdf(file=paste("/home/tsekara/Novel_ORFs_Prediction/Codon_periodicity//RiboseqR/GFP_2_Plots/",gsub("txt","pdf",gfp2[i]),sep=""))
  plotFS(fS)
  dev.off()
}