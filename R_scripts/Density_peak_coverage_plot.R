coverage.df <-
  read.table("target_coverage_exp8_tr5_wig.wig",skip = 1, col.names = c("target","start","end","reads"));
t.name = "tr5_comp71_c0_seq1";
mygene.coverage.df <- subset(coverage.df, target == t.name);
par(mar = rep(2, 4))
plot(x = 0, type = "n", xlim = c(0,max(mygene.coverage.df$end)),
     ylim = c(0,max(mygene.coverage.df$reads)),
     xlab = sprintf("Transcript location (%s)",t.name), ylab = "Read count");
segments(x0 = mygene.coverage.df$start, x1 = mygene.coverage.df$end,
         y0 = mygene.coverage.df$reads, lwd = 2);
segments(x0 = mygene.coverage.df$end,
         x1 = c(mygene.coverage.df$start[-1], max(mygene.coverage.df$end)),
         y0 = mygene.coverage.df$reads,
         y1 = c(mygene.coverage.df$reads[-1],0), lwd = 2);