# Small function that takes data.frame as input to generate scatterplots. 
# Run the command has scatterplots(data.frame.filename)


scatterplots_inonepdf=function(file)
{

  IDs <- colnames(file)
 pdf(file="Scatterplots.pdf",onefile=TRUE) 
  for (i in 1:(dim(file)[2]-1))
  {
    for( j in i:(dim(file)[2]) )
    {
    
      if (i != j)
      {          
       
        correlation <- cor(file[,i],file[,j])
        maximum <- max(log2(file[,i]))
        minimum <- min(log2(file[,i]))
        plot(log2(file[,i]),log2(file[,j]),xlab=IDs[i],ylab=IDs[j],pch='.',main =paste0(IDs[i],"_vs_",IDs[j]) )
        title(sub=paste0("R= ", cor(file[,i],file[,j])))
      }
    }  
  }  
   dev.off()     
}