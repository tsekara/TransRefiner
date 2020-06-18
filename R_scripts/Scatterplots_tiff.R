scatterplots_tiff=function(file)
{
  
  IDs <- colnames(file)
  for (i in 1:(dim(file)[2]-1))
  {
    for( j in i:(dim(file)[2]) )
    {
      
      if (i != j)
      {          
        tiff(file=paste(IDs[i],"_vs_",IDs[j],".tiff",sep=""),height = 12, width = 17, units = 'cm',compression = "lzw", res = 400) 
        correlation <- cor(file[,i],file[,j])
        maximum <- max(log2(file[,i]))
        minimum <- min(log2(file[,i]))
        plot(log2(file[,i]),log2(file[,j]),xlab=IDs[i],ylab=IDs[j],pch='.',main =paste0(IDs[i],"_vs_",IDs[j]) )
        title(sub=paste0("R= ", cor(file[,i],file[,j])))
        dev.off() 
      }
    }  
  }  
  
}