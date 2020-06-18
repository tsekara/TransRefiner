function (data,annotation, p_method = "correlation", p_link = "average", p_zscore = FALSE, p_heatmap = TRUE, main_annot = "Clustering", classes = colnames(data), col_size = 0.3, row_size = 0.3, genes = rownames(data), p_classes = colnames(data), c_margins = c (6, 6))
{
require (amap)
require(marray)
zscore<-function (row)
{
mw<-mean (row, na.rm= TRUE)
std<-sd (row, na.rm= TRUE)
row<- (row-mw)/std
return (row)
}
col.pal<-maPalette (low = "darkblue", high = "darkred", mid = "white", k = 1000)
#standardize rows
if(p_zscore==TRUE)
{
clusterdata<- t(apply (data, 1, zscore))
rownames (clusterdata) <- annotation[as.character(rownames (data)), 1]
}
else
{
clusterdata<- as.matrix(data)
rownames(clusterdata) <- annotation[as.character(rownames(data)),1]
}
# cluster data
clustersamples<-as.dendrogram(hcluster (t (clusterdata), method = p_method, link = p_link))
clustergenes<-as.dendrogram(hcluster (clusterdata, method = p_method, link = p_link))
col.pal<-maPalette (low = "darkblue", high = "darkred", mid = "white", k = 1000)
csc<-rainbow (max (max (as.integer(factor (p_classes)))))[as.integer( factor (p_classes))]
heatmap (clusterdata, Rowv = clustergenes, Colv = clustersamples, col = col.pal, main = main_annot, labRow = rownames (clusterdata), labCol = classes, cexCol = col_size, cexRow = row_size, ColSideColors = csc, margins = c_margins)    
detach("package:amap")
}
