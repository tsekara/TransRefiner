planarians_evalue=do.call(rbind,lapply(split(planarians,planarians$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(planarians_evalue)
planarians_cov_grt_50=do.call(rbind,lapply(split(planarians,planarians$Trinity5_Id),function(x) x[(x$Coverage_percentage >= 50),]))
dim(planarians_cov_grt_50)
planarians_cov_grt_50_evalue=do.call(rbind,lapply(split(planarians_cov_grt_50,planarians_cov_grt_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(planarians_cov_grt_50_evalue)
planarians_cov_less_50=do.call(rbind,lapply(split(planarians,planarians$Trinity5_Id),function(x) x[(x$Coverage_percentage < 50),]))
dim(planarians_cov_less_50)
planarians_cov_less_50_evalue=do.call(rbind,lapply(split(planarians_cov_less_50,planarians_cov_less_50$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(planarians_cov_less_50_evalue)
planarians_cov_evalue=rbind(planarians_cov_grt_50_evalue,planarians_cov_less_50_evalue)
planarians_cov_evalue=do.call(rbind,lapply(split(planarians_cov_evalue,planarians_cov_evalue$Trinity5_Id),function(x) x[which.min(x$E_value),]))
dim(planarians_cov_evalue)