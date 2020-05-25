# Load the Groped blast hits into R and provide the respective objects into following function

orig_extn_comp=function(original, extended)
{
original$avg_per_id=NULL
original$max_pct_len=NULL
original$query_match_range=NULL
original$target_match_range=NULL
original$pct_query_len=NULL

extended$avg_per_id=NULL
extended$query_match_range=NULL
extended$target_match_range=NULL
extended$pct_query_len=NULL
extended$max_pct_len=NULL
extended$max_pct_len=NULL

orig_extn_with_evalue_target_len=merge(original,extended,by=c("X.query_acc","target_acc"))

orig_extn_with_evalue_target_len=as.data.frame(cbind(orig_extn_with_evalue_target_len,
Extn.cov_greater_than_original.cov= ifelse((orig_extn_with_evalue_target_len$pct_target_len.y > 
orig_extn_with_evalue_target_len$pct_target_len.x), T, F),
Extn.cov_greater_than_original.cov_and_Extn.eval_less_than_Orig.eval= ifelse((orig_extn_with_evalue_target_len$pct_target_len.y > orig_extn_with_evalue_target_len$pct_target_len.x & 
orig_extn_with_evalue_target_len$min_Evalue.y < orig_extn_with_evalue_target_len$min_Evalue.x), T, F)))

write.table(orig_extn_with_evalue_target_len,"Comparisions_of_extended.and.originals_vs_Pooled_planmine_transcriptomes_with_pct_target_eval.txt",sep="\t",row.names = F,quote=F)

Extn.cov_greater_than_Orig.cov=orig_extn_with_evalue_target_len[grep("TRUE",orig_extn_with_evalue_target_len$Extn.cov_greater_than_original.cov),]
number_of_uniq_trans_Extn.cov_greater_than_Orig.cov= length(unique(Extn.cov_greater_than_Orig.cov$X.query_acc))


Extn.cov_and_eval_greater_than_Orig.cov_eval=orig_extn_with_evalue_target_len[grep("TRUE",orig_extn_with_evalue_target_len$Extn.cov_greater_than_original.cov_and_Extn.eval_less_than_Orig.eval),]
number_of_uniq_trans_Extn.cov_greater_than_Orig.cov_eval=length(unique(Extn.cov_and_eval_greater_than_Orig.cov_eval$X.query_acc))


message("Extn.cov > Orig.cov:",number_of_uniq_trans_Extn.cov_greater_than_Orig.cov)
message("Extn.cov > Orig.cov and Extn.eval < Orig.eval:",number_of_uniq_trans_Extn.cov_greater_than_Orig.cov_eval)

}


