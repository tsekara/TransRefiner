count.df <- read.csv("target_counts_tr5.csv");
count.mat <- xtabs(Reads ~ Target + Experiment, data = count.df);
count.mat <- count.mat[order(rowSums(count.mat), decreasing = TRUE),];
write.csv(count.mat,"tr5_counts_sam2profile.csv");
