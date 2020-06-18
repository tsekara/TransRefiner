
human_evalue=do.call(rbind,lapply(split(human,human$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(human_evalue)
human_cov_grt_50=do.call(rbind,lapply(split(human,human$Trinity5_Id),function(x) x[(x$Coverage_percentage >= 50),]))
dim(human_cov_grt_50)
human_cov_grt_50_evalue=do.call(rbind,lapply(split(human_cov_grt_50,human_cov_grt_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(human_cov_grt_50_evalue)
human_cov_less_50=do.call(rbind,lapply(split(human,human$Trinity5_Id),function(x) x[(x$Coverage_percentage < 50),]))
dim(human_cov_less_50)
human_cov_less_50_evalue=do.call(rbind,lapply(split(human_cov_less_50,human_cov_less_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(human_cov_less_50_evalue)
human_cov_evalue=rbind(human_cov_grt_50_evalue,human_cov_less_50_evalue)
human_cov_evalue=do.call(rbind,lapply(split(human_cov_evalue,human_cov_evalue$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(human_cov_evalue)

 