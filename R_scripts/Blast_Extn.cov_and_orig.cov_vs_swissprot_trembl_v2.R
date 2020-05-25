# Load the Grouped blast hits into R and provide the respective objects into following function

orig_extn_comp=function(original, extended)
{
  library(data.table)
  library(dplyr)
  
  
  original=as.data.frame(setDT(original)[, .SD[which.min(min_Evalue)], by=X.query_acc]) 
  #extended=as.data.frame(setDT(extended)[, .SD[which.max(pct_query_len)], by=X.query_acc])  
  
  original=original %>% select(X.query_acc,target_acc,min_Evalue,pct_query_len,pct_target_len)
  colnames(original)=paste("orig",colnames(original),sep="_")
  
  
  extended=extended %>% select(X.query_acc,target_acc,min_Evalue,pct_query_len,pct_target_len)
  colnames(extended)=paste("extn",colnames(extended),sep="_")
  
  
  
  merged_orig_extn_with_evalue_target_len=merge(original,extended,by.x=c("orig_X.query_acc","orig_target_acc"),by.y = c("extn_X.query_acc","extn_target_acc"))
  
  orig_extn_with_evalue_target_len=as.data.frame(cbind(merged_orig_extn_with_evalue_target_len,
                                                       Extn.cov_greater_than_original.cov= ifelse((merged_orig_extn_with_evalue_target_len$extn_pct_target_len > 
                                                                                                     merged_orig_extn_with_evalue_target_len$orig_pct_target_len), T, F),
                                                       Extn.cov_equals_original.cov= ifelse((merged_orig_extn_with_evalue_target_len$extn_pct_target_len == 
                                                                                                    merged_orig_extn_with_evalue_target_len$orig_pct_target_len), T, F),
                                                       Extn.cov_less_than_original.cov= ifelse((merged_orig_extn_with_evalue_target_len$extn_pct_target_len < 
                                                                                                  merged_orig_extn_with_evalue_target_len$orig_pct_target_len), T, F),
                                                       Extn.cov_greater_than_original.cov_and_Extn.eval_less_than_Orig.eval= ifelse((merged_orig_extn_with_evalue_target_len$extn_pct_target_len > merged_orig_extn_with_evalue_target_len$orig_pct_target_len & 
                                                                                                                                       merged_orig_extn_with_evalue_target_len$extn_min_Evalue < merged_orig_extn_with_evalue_target_len$orig_min_Evalue), T, F)))
  
  write.table(orig_extn_with_evalue_target_len,"Comparisions_of_extended.and.originals_vs_Uniprot_with_pct_target_eval.txt",sep="\t",row.names = F,quote=F)
  
  
  number_of_uniq_trans_Extn.cov_greater_than_Orig.cov=orig_extn_with_evalue_target_len %>% 
                                                      select(orig_X.query_acc,Extn.cov_greater_than_original.cov) %>% 
                                                      filter(Extn.cov_greater_than_original.cov==T) %>% n_distinct()
  
    
  
  number_of_uniq_trans_Extn.cov_equalto_Orig.cov=orig_extn_with_evalue_target_len %>% 
                                                  select(orig_X.query_acc,Extn.cov_equals_original.cov) %>% 
                                                  filter(Extn.cov_equals_original.cov == T) %>% n_distinct()
  
  
  number_of_uniq_trans_Extn.cov_greater_than_Orig.cov_eval= orig_extn_with_evalue_target_len %>%
                                                            select(orig_X.query_acc,Extn.cov_greater_than_original.cov_and_Extn.eval_less_than_Orig.eval) %>%
                                                            filter(Extn.cov_greater_than_original.cov_and_Extn.eval_less_than_Orig.eval == T) %>% n_distinct()
  
  message("Extn.cov > Orig.cov:",number_of_uniq_trans_Extn.cov_greater_than_Orig.cov)
  message("Extn.cov > Orig.cov and Extn.eval < Orig.eval:",number_of_uniq_trans_Extn.cov_greater_than_Orig.cov_eval)
  message("Extn.cov = Orig.cov:",number_of_uniq_trans_Extn.cov_equalto_Orig.cov)

}

