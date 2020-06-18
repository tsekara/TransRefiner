mouse_orfs=NULL
for(i in multiple_orfs)
{
mouse_orfs=rbind(mouse_orfs,mouse[grepl(i,mouse$Trinity5_Id),])
}
human_orfs=NULL
for(i in multiple_orfs)
{
human_orfs=rbind(human_orfs,human[grepl(i,human$Trinity5_Id),])
}
dmel_orfs=NULL
for(i in multiple_orfs)
{
dmel_orfs=rbind(dmel_orfs,d_mel[grepl(i,d_mel$Trinity5_Id),])
}
celegans_orfs=NULL
{
celegans_orfs=rbind(celegans_orfs,c_elegans[grepl(i,c_elegans$Trinity5_Id),])
}
planarians_orfs=NULL
for(i in multiple_orfs)
{
planarians_orfs=rbind(planarians_orfs,planarains[grepl(i,planarains$Trinity5_Id),])
}
smansoni_orfs=NULL
{
smansoni_orfs==rbind(smansoni_orfs,smansoni[grepl(i,smansoni$Trinity5_Id),])
}