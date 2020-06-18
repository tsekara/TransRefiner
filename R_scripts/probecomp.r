fn=function(fnames)
{
strp=unlist(strsplit(fnames,split=","))
n= numeric(length(strp))
setwd("D:/combine")
for(i in seq_along(strp))
{
    d = read.delim(strp[i])
    n[i] = dim(d)[1]
}
plot(n)
return(n)
}
