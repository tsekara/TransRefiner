

ProbeIDs_Comparisions=function(file)
{        
my.file=file     
a=0
b=0
str=unlist(strsplit(file,split=","))
 for(i in 1:length(str))
  
  {       
          i=1
        setwd="D:/combine"
          for(j in 2:length(str))
          {
          j=2
          f=read.delim(str[i],row.names=1,dec=".")
          g=read.delim(str[j],row.names=1,dec=".")        
          a[i]=nrow(f)  
          b[j]=nrow(g)                
           }
           if(a[i]!=b[j])
           {
           # stop("Ids Mismatch")
           break
           }
                   
  }
     if(a[i]==b[j])
         {
           print ("IDs match")
         }
     else
         {
           print ("IDs does not match")
         }
}


  









