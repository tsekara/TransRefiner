################################ Loading the Blast hits for every specices #######################################
c_elegans=read.delim("/home/tsekara/Annotation_Database/October_2015/C_elegans/BLATSp_C_elegans/C_elegans_cov_evalue_annotations_subset.txt")
smansoni=read.delim("/home/tsekara/Annotation_Database/October_2015/S_mansoni/BLASTp_S_mansoni/S_mansoni_cov_evlaue_anno_subset.txt")
human=read.delim("/home/tsekara/Annotation_Database/October_2015/Human/BLASTp_Human/Human_cov-evalue_annotations_subset.txt")
dmel=read.delim("/home/tsekara/Annotation_Database/October_2015/D_mel/BLASTp_Dmel/dmel_cov_evalue_anno_subset.txt")
mouse=read.delim("/home/tsekara/Annotation_Database/October_2015/Mouse/BLASTp_Mouse/Mouse_cov_evlaue_annotations_subset.txt")
planarains=read.delim("/home/tsekara/Annotation_Database/October_2015/Planarians/BLASTp_Planarians/planarians_cov_evalue_anno_subset.txt")


############################# Combining the DEGs and blast hits ###############################
# DEGs list along with their log2FoldChange is created either seperatley as txt file and loaded into R or within R #####################

deg=read.delim("/path/to/DEG_file_with_log2foldchange")
deg_blast_hits=merge(merge(merge(merge(merge(merge(deg,human,by="Trinity5_Id",all.x=T),mouse,by="Trinity5_Id",all.x=T),c_elegans,by="Trinity5_Id",all.x=T),dmel,by="Trinity5_Id",all.x=T),planarains,by="Trinity5_Id",all.x=T),smansoni,by="Trinity5_Id",all.x=T)
write.table(deg_blast_hits,".../DEG_BLAST_hits.txt",sep="\t",row.names=F,quote=F)

################################### Load the list of trasncripts which has multiple ORFs. This can be found by comparing the list of DEGs for which 
### BLAST hits has to appended and list of trasncripts which has multiple ORFs in trinity5_framedp_corrected_pepdb ###############################

multiple_ORFs_ids=read.delim("/home/tsekara/Annotation_Database/October_2015/Multiple_ORFs_transcripts_Ids.txt")
dim(multiple_ORFs_ids)
deg_multiple_orfs=as.character(intersect(deg$Trinity5_Id,multiple_ORFs_ids$Trinity5_Id))
head(deg_multiple_orfs)

################################# Species specific mulitple ORFs integration #######################################
########### MOUSE ###################         
mouse_orfs=NULL
for(i in deg_multiple_orfs)
{
  mouse_orfs=rbind(mouse_orfs,mouse[grepl(i,mouse$Trinity5_Id),])
}
########### HUMAN ###################
human_orfs=NULL
for(i in deg_multiple_orfs)
{
  human_orfs=rbind(human_orfs,human[grepl(i,human$Trinity5_Id),])
}
########## D_mel ##################
dmel_orfs=NULL
for(i in deg_multiple_orfs)
{
  dmel_orfs=rbind(dmel_orfs,d_mel[grepl(i,d_mel$Trinity5_Id),])
}
####### C_elegans #############
celegans_orfs=NULL
for(i in deg_multiple_orfs)
{
  celegans_orfs=rbind(celegans_orfs,c_elegans[grepl(i,c_elegans$Trinity5_Id),])
}
######## Planarians ################
planarians_orfs=NULL
for(i in multiple_orfs)
{
  planarians_orfs=rbind(planarians_orfs,planarains[grepl(i,planarains$Trinity5_Id),])
}

############## S_mansoni #########
smansoni_orfs=NULL
for(i in multiple_orfs)
{
  smansoni_orfs=rbind(smansoni_orfs,smansoni[grepl(i,smansoni$Trinity5_Id),])
}

############################### Merging all species ORFs into singel dataframe #############

multiple_orfs=merge(merge(merge(merge(merge(human_orfs,mouse_orfs,by="Trinity5_Id",all=T),celegans_orfs,by="Trinity5_Id",all=T),dmel_orfs,by="Trinity5_Id",all=T),planarians_orfs,by="Trinity5_Id",all=T),smansoni_orfs,by="Trinity5_Id",all=T)

############################# creating a empty column, Log2FoldChange in multiple_orfs data.frame ####################

multiple_orfs$log2FoldChange=NA

##########################  Rearranging columns suitable for binding to the exisitng deg_blast_hits dataframe #################

new_multiple_orfs=multiple_orfs[,c(1,26,2:25)]

############################ Row-binding the new_multipleorfs with existing deg_blast_hits ###########################

deg_blast_hits_multiple_orfs=rbind(deg_blast_hits,new_multiple_orfs)


write.table(deg_blast_hits_multiple_orfs,"DEG_blast_hits_multiple_orfs.txt",sep="\t",row.names=F,quote=F)






