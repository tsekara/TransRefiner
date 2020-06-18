res <- do.call(rbind.data.frame,lapply(strsplit(as.character(t5_GO$GOs), "; "),
  function(x){
  x1 <- tapply(x, sub(':.*', '', x), FUN=paste, collapse=";")
  x1[match(c('C', 'F', 'P'),  names(x1))]}))
res1 <-  data.frame(Name=t5_GO[,1], setNames(res, c('Cellular_Component','Molecular_Function', 'Biological_Process')), stringsAsFactors=FALSE)