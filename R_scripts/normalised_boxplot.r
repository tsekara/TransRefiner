combine=function(file)
{

split_list <- unlist(strsplit(file,split=","))
setwd("D:/combine")
dataset <- do.call("cbind",lapply(split_list,FUN=function(files)
{
read.table(files,header=TRUE, sep="\t")
}
)
)
names(dataset)[1]=paste("Probe_ID")
 drop=c("ProbeID")
dataset=dataset[,!(names(dataset)%in%drop)]

write.table(dataset,file ="D:/output/illumina.txt",sep="\t",col.names=T,row.names=F)

      raw.Data=read.delim(file="D:/output/illumina.txt",row.names = 1, dec = ".")
      raw.expression <- raw.Data[,seq(1,dim(raw.Data)[2],2)]
      dim(raw.expression)   
      raw.calls <- raw.Data[,seq(2,dim(raw.Data)[2],2)]
      dim(raw.calls)

      library (limma)
      data.quantiles<- normalizeBetweenArrays (as.matrix (raw.expression), method = "quantile")
      colnames (data.quantiles) <- colnames (raw.expression)
      rownames (data.quantiles) <- rownames (raw.expression)   
      
   
      setwd("D:/output")
      pdf(file=paste("_boxplot.pdf",sep="")) 
      boxplot(data.quantiles,log="y", pch='.',col=rainbow(ncol(data.quantiles)))
      dev.off()  
      return("_boxplot.pdf")


}