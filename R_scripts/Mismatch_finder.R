mismatch_finder=function(wt,mt,df)
{
  for (i in 1:length(names(wt))) {
  df$ID= print(names(wt)[i])
  df$AG=  print(length(intersect(which(wt[[i]] == "a"), which(mt[[i]] == "g"))))
  df$GA=  print(length(intersect(which(wt[[i]] == "g"), which(mt[[i]] == "a"))))
  df$CT=  print(length(intersect(which(wt[[i]] == "t"), which(mt[[i]] == "c"))))
  df$TC=  print(length(intersect(which(wt[[i]] == "c"), which(mt[[i]] == "t"))))
                    
  }

}



