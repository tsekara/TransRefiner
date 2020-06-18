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
      
      correlation.matrix=cor( raw.expression )
      library("gplots")
      library("RColorBrewer")
      col.pal<-colorRampPalette (brewer.pal (10, "RdYlBu")) (256)
      setwd("D:/output")
       pdf(file=paste("Heatmap.pdf",sep="")) 
      heatmap.2(correlation.matrix, col = col.pal, key =TRUE, symkey= FALSE, density.info = "none", trace = "none", margins = c(10,10))
      dev.off()
      return("Heatmap.pdf")
      }