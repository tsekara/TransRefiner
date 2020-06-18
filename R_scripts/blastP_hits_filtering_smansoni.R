smansoni_evalue=do.call(rbind,lapply(split(smansoni,smansoni$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(smansoni_evalue)

smansoni_cov_grt_50=do.call(rbind,lapply(split(smansoni,smansoni$Trinity5_Id),function(x) x[(x$Coverage_percentage >= 50),]))
dim(smansoni_cov_grt_50)
smansoni_cov_grt_50_evalue=do.call(rbind,lapply(split(smansoni_cov_grt_50,smansoni_cov_grt_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(smansoni_cov_grt_50_evalue)

smansoni_cov_less_50=do.call(rbind,lapply(split(smansoni,smansoni$Trinity5_Id),function(x) x[(x$Coverage_percentage < 50),]))
dim(smansoni_cov_less_50)
smansoni_cov_less_50_evalue=do.call(rbind,lapply(split(smansoni_cov_less_50,smansoni_cov_less_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(smansoni_cov_less_50_evalue)

smansoni_cov_evalue=rbind(smansoni_cov_grt_50_evalue,smansoni_cov_less_50_evalue)
smansoni_cov_evalue=do.call(rbind,lapply(split(smansoni_cov_evalue,smansoni_cov_evalue$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(smansoni_cov_evalue)


write.table(blastp_planarians_best_hits,"blastP_planarians_best_hits.txt",sep="\t",row.names=F,quote=F)
