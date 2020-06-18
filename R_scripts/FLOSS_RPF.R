library("ggplot2")



source("mouse-annotation.R")
minMatch=26
maxMatch=34
fileRainer="Ctrl1"

funcFloss=function(fileRainer){
  dataRainer=read.table(paste(fileRainer,"_covMatchLen.txt",sep=""),h=T,sep="\t")
  
  #ixColFloss=grep("^Floss",colnames(dataRainer))
  ixMin=which(colnames(dataRainer)==paste("Floss.",minMatch,sep=""))
  ixMax=which(colnames(dataRainer)==paste("Floss.",maxMatch,sep=""))
  ixColFloss=ixMin:ixMax
  
  dataPC_chosen=dataRainer[na.omit(match(refCodingIds,dataRainer$Transcript)),]
  
  #the proteinCoding of reference
  dataRef=dataPC_chosen[which(dataPC_chosen$Transcript.Type=="protein_coding"),ixColFloss]
  
  #these are genes that have annotation conflicts between my version and Ensembl new version
  conflictGenes=dataPC_chosen$Transcript[which(dataPC_chosen$Transcript.Type!="protein_coding")]
  
  rownames(dataRef)=as.character(dataPC_chosen$Transcript[which(dataPC_chosen$Transcript.Type=="protein_coding")])
  dataRef_num=matrix(as.numeric(as.character(unlist(dataRef))),nrow=nrow(dataRef))
  rownames(dataRef_num)=rownames(dataRef)
  colnames(dataRef_num)=colnames(dataRef)
  
  
  allData=dataRainer[-match(conflictGenes,dataRainer$Transcript),]
  keepData=allData[,ixColFloss]
  names_AllData=apply(allData[,c(2,10,11)],1,paste,sep="",collapse="")
  rownames(keepData)=names_AllData
  keepData_num=matrix(as.numeric(as.character(unlist(keepData))),nrow=nrow(keepData))
  rownames(keepData_num)=rownames(keepData)
  colnames(keepData_num)=colnames(keepData)
  
  ## L1 norm length distribution difference between two length counts
  ldistDiff <- function(d1, d2) { sum(abs(colSums(rbind(d1/sum(d1), - d2/sum(d2))))) / 2.0 }
  
  refSums=colSums(na.omit(dataRef_num))
  dataGGP_MLength=cbind(gsub("Floss.","", names(refSums)),refSums)
  colnames(dataGGP_MLength)=c("matchSize","counts")
  dataGGP_MLength=data.frame(dataGGP_MLength)
  dataGGP_MLength$counts=as.numeric(as.character(dataGGP_MLength$counts))
  #plot the distribution of the reference
  barPlotMatchSize = ggplot(dataGGP_MLength, aes(matchSize,counts)) + 
    geom_bar(stat = "identity", position = "dodge",color = "steelblue",fill = "steelblue")+
    labs(x = "Size of the match", y = paste(fileRainer,"Ref coverage per match size",sep=""),fill = NULL)+
    theme(axis.title.y = element_text(family = "serif",size = 15,face = "bold"))+
    theme(axis.title.x = element_text(family = "serif",size = 15,face = "bold"))
  print(barPlotMatchSize)
  
  floss <- apply(keepData_num, 1, function(l) { ldistDiff(l,refSums) })
  
  dataRainer_floss=cbind(allData,floss)
  
  plot(density(na.omit(floss)),main=paste(fileRainer,"Floss density",sep=""),bty="n",col="black",lwd=2,ylim=c(0,7))
  lines(density(na.omit(dataRainer_floss$floss[which(dataRainer_floss$Transcript.Type=="protein_coding")])),lwd=2,col="red")
  lines(density(na.omit(dataRainer_floss$floss[which(dataRainer_floss$Transcript.Type=="lincRNA")])),lwd=2,col="blue")
  lines(density(na.omit(dataRainer_floss$floss[which(dataRainer_floss$Transcript.Type!="protein_coding")])),lwd=2,col="darkgreen")
  lines(density(na.omit(dataRainer_floss$floss[which(dataRainer_floss$Transcript.Type=="snoRNA")])),lwd=2,col="orange")
  #the sno don't work.
  legend("topright",legend=c("all","protein_coding","lincRNA","nonPC","snoRNA"),lwd=2,bty="n",col=c("black","red","blue","darkgreen","orange"))
  
  #plot(rowSums(na.omit(dataRainer_floss[,ixColFloss])),dataRainer_floss$floss)
  return(dataRainer_floss)
}

pdf("LengthDist_FlossScore.pdf")
dataRainer1=funcFloss("Ctrl1")
dataRainer2=funcFloss("Ctrl2")
dev.off()

mergeDataRainer=merge(dataRainer1,dataRainer2,by=c(2,10,11))

meanFloss=apply(cbind(mergeDataRainer$floss.x,mergeDataRainer$floss.y),1,mean)

mergeDataRainer=cbind(mergeDataRainer,meanFloss)
