library(superheat)  # super heatmap function
source("functions_helper.R")

# # # 

# # # statistical correlation

my_table_app_main <- NULL
my_table_sim_train <- NULL
for (pop in pop_list){
  # train dataset
  my_table_sim_train[[pop]] <- load_data(my_files_sim_train, pop)
  my_table_sim_train[[pop]] <- preprocess_data(my_table_sim_train[[pop]])
  
  # application dataset
  my_table_app_main[[pop]] <- load_data(my_files_app_main, pop)
  my_table_app_main[[pop]] <- preprocess_data(my_table_app_main[[pop]])
}

for (pop in pop_list){
  stat_temp <- c(stat_misc, stat_list)
  
  # neutral simulation correlation matrix
  cor_neutral <- filter(my_table_sim_train[[pop]], site_class=="neutral")[stat_temp] 
  cor_neutral <- cor(cor_neutral, use="pairwise.complete.obs", method="spearman")
  colnames(cor_neutral) <- stat_code[stat_temp]; rownames(cor_neutral) <- stat_code[stat_temp]
  # pdf(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-neutral", plout, sep=""))
  tiff(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-neutral", ".tiff", sep=""),
       height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # heatmap.2(abs(cor_neutral), scale="none", col=brewer.pal(9, "Blues"), margins=c(15, 15), 
  #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  superheat(abs(cor_neutral), row.dendrogram=T, col.dendrogram=T, 
            grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
            left.label.col="white", bottom.label.col="white", left.label.text.size=3, left.label.text.alignment="right", 
            bottom.label.text.size=3, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  dev.off()
  
  tab_neutral <- data.frame(row=rownames(cor_neutral)[row(cor_neutral)[upper.tri(cor_neutral)]], 
    col= colnames(cor_neutral)[col(cor_neutral)[upper.tri(cor_neutral)]], 
    corr=cor_neutral[upper.tri(cor_neutral)]) %>% arrange(-corr)
  write.table(tab_neutral, paste("../../results/feature_plots_cor_txt/", label_sim_train, "-", 
                                 pop, "-spearman-cor-neutral.txt", sep=""), row.names=F)
  
  # hardsweep simulation correlation matrix
  cor_hardsweep <- filter(my_table_sim_train[[pop]], site_class=="hardsweep")[stat_temp] 
  cor_hardsweep <- cor(cor_hardsweep, use="pairwise.complete.obs", method="spearman")
  colnames(cor_hardsweep) <- stat_code[stat_temp]; rownames(cor_hardsweep) <- stat_code[stat_temp]
  # pdf(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-hardsweep", plout, sep=""))
  tiff(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-hardsweep", ".tiff", sep=""),
       height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # heatmap.2(abs(cor_hardsweep), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  superheat(abs(cor_hardsweep), row.dendrogram=T, col.dendrogram=T, 
            grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
            left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
            bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  dev.off()
  
  tab_hardsweep <- data.frame(row=rownames(cor_hardsweep)[row(cor_hardsweep)[upper.tri(cor_hardsweep)]], 
    col= colnames(cor_hardsweep)[col(cor_hardsweep)[upper.tri(cor_hardsweep)]], 
    corr=cor_hardsweep[upper.tri(cor_hardsweep)]) %>% arrange(-corr)
  write.table(tab_hardsweep, paste("../../results/feature_plots_cor_txt/", label_sim_train, "-", 
                                   pop, "-spearman-cor-hardsweep.txt", sep=""), row.names=F)
  
  # softsweep simulation correlation matrix
  cor_softsweep <- filter(my_table_sim_train[[pop]], site_class=="softsweep")[stat_temp] 
  cor_softsweep <- cor(cor_softsweep, use="pairwise.complete.obs", method="spearman")
  colnames(cor_softsweep) <- stat_code[stat_temp]; rownames(cor_softsweep) <- stat_code[stat_temp]
  # pdf(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-softsweep", plout, sep=""))
  tiff(file=paste("../../results/feature_plots_cor_fig/", label_sim_train, "-", pop, "-spearman-cor-softsweep", ".tiff", sep=""),
       height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # heatmap.2(abs(cor_softsweep), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  superheat(abs(cor_softsweep), row.dendrogram=T, col.dendrogram=T, 
            grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
            left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
            bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  dev.off()
  
  tab_softsweep <- data.frame(row=rownames(cor_softsweep)[row(cor_softsweep)[upper.tri(cor_softsweep)]], 
    col= colnames(cor_softsweep)[col(cor_softsweep)[upper.tri(cor_softsweep)]], 
    corr=cor_softsweep[upper.tri(cor_softsweep)]) %>% arrange(-corr)
  write.table(tab_softsweep, paste("../../results/feature_plots_cor_txt/", label_sim_train, "-", 
                                   pop, "-spearman-cor-softsweep.txt", sep=""), row.names=F)
  
  # application correlation matrix
  cor_app <- filter(my_table_app_main[[pop]], T)[stat_temp]
  cor_app <- cor(cor_app, use="pairwise.complete.obs", method="spearman")
  colnames(cor_app) <- stat_code[stat_temp]; rownames(cor_app) <- stat_code[stat_temp]
  # pdf(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app", plout, sep=""))
  tiff(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app", ".tiff", sep=""),
       height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # heatmap.2(abs(cor_app), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  superheat(abs(cor_app), row.dendrogram=T, col.dendrogram=T, 
            grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
            left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
            bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  dev.off()
  
  tab_app <- data.frame(row=rownames(cor_app)[row(cor_app)[upper.tri(cor_app)]], 
    col= colnames(cor_app)[col(cor_app)[upper.tri(cor_app)]], 
    corr=cor_app[upper.tri(cor_app)]) %>% arrange(-corr)
  write.table(tab_app, paste("../../results/feature_plots_cor_txt/", label_app_main, "-", 
                             pop, "-spearman-cor-app.txt", sep=""), row.names=F)
  
  # # difference between app and neutral
  # # pdf(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-neutral", plout, sep=""))
  # tiff(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-neutral", ".tiff", sep=""),
  #      height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # # heatmap.2(abs(cor_app-cor_neutral), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  # #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  # superheat(abs(cor_app-cor_neutral), row.dendrogram=T, col.dendrogram=T, 
  #           grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
  #           left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
  #           bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  # dev.off()
  # 
  # # difference between app and softsweep
  # # pdf(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-softsweep", plout, sep=""))
  # tiff(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-softsweep", ".tiff", sep=""),
  #      height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # # heatmap.2(abs(cor_app-cor_softsweep), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  # #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  # superheat(abs(cor_app-cor_softsweep), row.dendrogram=T, col.dendrogram=T, 
  #           grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
  #           left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
  #           bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  # dev.off()
  # 
  # # difference between app and hardsweep
  # # pdf(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-hardsweep", plout, sep=""))
  # tiff(file=paste("../../results/feature_plots_cor_fig/", label_app_main, "-", pop, "-spearman-cor-app-hardsweep", ".tiff", sep=""),
  #      height=7.5, width=6.5, units="in", compression="lzw", res=300)
  # # heatmap.2(abs(cor_app-cor_hardsweep), scale="none", col=brewer.pal(9, "Blues"), margins=c(15 ,15), 
  # #           tracecol=NA, density.info="none", breaks=seq(0,1,length.out=10), main=pop)
  # superheat(abs(cor_app-cor_hardsweep), row.dendrogram=T, col.dendrogram=T, 
  #           grid.hline.col="white", grid.vline.col="white", heat.pal=brewer.pal(9, "Blues"), heat.lim=c(0,1), 
  #           left.label.col="white", bottom.label.col="white", left.label.text.size=3.5, left.label.text.alignment="right", 
  #           bottom.label.text.size=3.5, bottom.label.text.alignment="right", bottom.label.text.angle=90, title=pop)
  # dev.off()
}
my_table_app_main <- NULL
my_table_sim_train <- NULL

# # # statistical power analysis
# # # & kolmogovor-smirnov analysis

pow <- 0.99

stat_pow_helper <- function(tmp, temp, level) {
  resl <- NULL
  resl$right <- sum(temp > quantile(tmp$ne, probs=level))/length(temp)  # right-tail
  resl$left <- sum(temp < quantile(tmp$ne, probs=1-level))/length(temp)  # left-tail
  resl$two <- (sum(temp > quantile(tmp$ne, probs=(1+level)/2)) + 
                 sum(temp < quantile(tmp$ne, probs=(1-level)/2)))/length(temp)  # two-tail
  resl$maxp <- max(c(resl$right, resl$left, resl$two))
  resl$argmax <- c("right-tailed", "left-tailed", "two-tailed")[which.max(c(resl$right, resl$left, resl$two))]
  return(resl)
}

statistical_power <- function(my_table, my_stat, power) {
  tmp <- NULL
  tmp$neutral <- na.omit(filter(my_table, site_class=="neutral")[[my_stat]])
  tmp$hardsweep <- na.omit(filter(my_table, site_class=="hardsweep")[[my_stat]])
  tmp$softsweep <- na.omit(filter(my_table, site_class=="softsweep")[[my_stat]])
  res <- NULL
  res$neutral <- stat_pow_helper(tmp, tmp$neutral, power)
  res$hardsweep <- stat_pow_helper(tmp, tmp$hardsweep, power)
  res$softsweep <- stat_pow_helper(tmp, tmp$softsweep, power)
  return(res)
}

ks_test_helper <- function(my_table, my_stat) {
  tmp <- NULL
  tmp$neutral <- na.omit(filter(my_table, site_class=="neutral")[[my_stat]])
  tmp$hardsweep <- na.omit(filter(my_table, site_class=="hardsweep")[[my_stat]])
  tmp$softsweep <- na.omit(filter(my_table, site_class=="softsweep")[[my_stat]])
  res <- NULL
  ks_temp <- ks.test(tmp$neutral, tmp$hardsweep)
  res$neutral_hardsweep$statistic <- ks_temp$statistic
  res$neutral_hardsweep$pvalue <- ks_temp$p.value
  ks_temp <- ks.test(tmp$neutral, tmp$softsweep)
  res$neutral_softsweep$statistic <- ks_temp$statistic
  res$neutral_softsweep$pvalue <- ks_temp$p.value
  ks_temp <- ks.test(tmp$hardsweep, tmp$softsweep)
  res$hardsweep_softsweep$statistic <- ks_temp$statistic
  res$hardsweep_softsweep$pvalue <- ks_temp$p.value
  return(res)
}

ks_test_app <- function(my_app, my_sim, my_stat) {
  tmp <- NULL
  tmp$app <- na.omit(my_app[[my_stat]])
  tmp$neutral <- na.omit(filter(my_sim, site_class=="neutral")[[my_stat]])
  tmp$hardsweep <- na.omit(filter(my_sim, site_class=="hardsweep")[[my_stat]])
  tmp$softsweep <- na.omit(filter(my_sim, site_class=="softsweep")[[my_stat]])
  res <- NULL
  ks_temp <- ks.test(tmp$app, tmp$neutral)
  res$neutral$statistic <- ks_temp$statistic
  res$neutral$pvalue <- ks_temp$p.value
  ks_temp <- ks.test(tmp$app, tmp$hardsweep)
  res$hardsweep$statistic <- ks_temp$statistic
  res$hardsweep$pvalue <- ks_temp$p.value
  ks_temp <- ks.test(tmp$app, tmp$softsweep)
  res$softsweep$statistic <- ks_temp$statistic
  res$softsweep$pvalue <- ks_temp$p.value
  return(res)
}

stat_power <- NULL
stat_kstest <- NULL
stat_ksapp <- NULL
my_table_app_main <- NULL
my_table_sim_train <- NULL
for (pop in pop_list){
  # app load
  my_table_app_main[[pop]] <- load_data(my_files_app_main, pop)
  my_table_app_main[[pop]] <- preprocess_data(my_table_app_main[[pop]])
  # train dataset
  my_table_sim_train[[pop]] <- load_data(my_files_sim_train, pop)
  my_table_sim_train[[pop]] <- preprocess_data(my_table_sim_train[[pop]])
  
  # statistical power
  stat_temp <- c(stat_misc, stat_list)
  for (stat in stat_temp) {
    stat_power[[pop]][[stat]][["all"]] <- statistical_power(my_table_sim_train[[pop]], stat, pow)
    
    stat_power[[pop]][[stat]][["complete"]] <- 
      statistical_power(filter(my_table_sim_train[[pop]], complete_bin%in%c("0","2")), stat, pow)
    stat_power[[pop]][[stat]][["incomplete"]] <- 
      statistical_power(filter(my_table_sim_train[[pop]], complete_bin%in%c("0","1")), stat, pow)
    
    # stat_power[[pop]][[stat]][["scoeff_high"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","3")), stat, pow)
    # stat_power[[pop]][[stat]][["scoeff_mid"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","2")), stat, pow)
    # stat_power[[pop]][[stat]][["scoeff_low"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","1")), stat, pow)
    # stat_power[[pop]][[stat]][["stime_high"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","3")), stat, pow)
    # stat_power[[pop]][[stat]][["stime_mid"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","2")), stat, pow)
    # stat_power[[pop]][[stat]][["stime_low"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","1")), stat, pow)
    # stat_power[[pop]][[stat]][["sfreq_high"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","3")), stat, pow)
    # stat_power[[pop]][[stat]][["sfreq_mid"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","2")), stat, pow)
    # stat_power[[pop]][[stat]][["sfreq_low"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","1")), stat, pow)
    
    # stat_power[[pop]][[stat]][["sprefix_low"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==1), stat, pow)
    # stat_power[[pop]][[stat]][["sprefix_mid"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==2), stat, pow)
    # stat_power[[pop]][[stat]][["sprefix_high"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==3), stat, pow)
    # stat_power[[pop]][[stat]][["spostfix_low"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==1), stat, pow)
    # stat_power[[pop]][[stat]][["spostfix_mid"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==2), stat, pow)
    # stat_power[[pop]][[stat]][["spostfix_high"]] <- 
    #   statistical_power(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==3), stat, pow)
    
  }
  
  # kolmogorov-smirnov test
  stat_temp <- c(stat_misc, stat_list)
  for (stat in stat_temp) {
    # compare simulation
    stat_kstest[[pop]][[stat]][["all"]] <- ks_test_helper(my_table_sim_train[[pop]], stat)
    
    stat_kstest[[pop]][[stat]][["complete"]] <- 
      ks_test_helper(filter(my_table_sim_train[[pop]], complete_bin%in%c("0","2")), stat)
    stat_kstest[[pop]][[stat]][["incomplete"]] <- 
      ks_test_helper(filter(my_table_sim_train[[pop]], complete_bin%in%c("0","1")), stat)
    
    # stat_kstest[[pop]][[stat]][["scoeff_high"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","3")), stat)
    # stat_kstest[[pop]][[stat]][["scoeff_mid"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","2")), stat)
    # stat_kstest[[pop]][[stat]][["scoeff_low"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], scoeff_bin%in%c("0","1")), stat)
    # stat_kstest[[pop]][[stat]][["stime_high"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","3")), stat)
    # stat_kstest[[pop]][[stat]][["stime_mid"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","2")), stat)
    # stat_kstest[[pop]][[stat]][["stime_low"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], stime_bin%in%c("0","1")), stat)
    # stat_kstest[[pop]][[stat]][["sfreq_high"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","3")), stat)
    # stat_kstest[[pop]][[stat]][["sfreq_mid"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","2")), stat)
    # stat_kstest[[pop]][[stat]][["sfreq_low"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], sfreq_bin%in%c("0","1")), stat)
    
    # stat_kstest[[pop]][[stat]][["sprefix_low"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==1), stat)
    # stat_kstest[[pop]][[stat]][["sprefix_mid"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==2), stat)
    # stat_kstest[[pop]][[stat]][["sprefix_high"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|sprefix_freq==3), stat)
    # stat_kstest[[pop]][[stat]][["spostfix_low"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==1), stat)
    # stat_kstest[[pop]][[stat]][["spostfix_mid"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==2), stat)
    # stat_kstest[[pop]][[stat]][["spostfix_high"]] <- 
    #   ks_test_helper(filter(my_table_sim_train[[pop]], site_class=="neutral"|spostfix_time==3), stat)
    
    # compare application
    stat_ksapp[[pop]][[stat]] <- ks_test_app(my_table_app_main[[pop]], my_table_sim_train[[pop]], stat)
  }
}

# convert to tibbles

subsets <- c("all", "complete", "incomplete")
# subsets <- c("all", "scoeff_high", "scoeff_mid", "scoeff_low", 
#   "stime_high", "stime_mid", "stime_low", "sfreq_high", 
#   "sfreq_mid", "sfreq_low", "complete", "incomplete",
#   "sprefix_high", "sprefix_mid", "sprefix_low",
#   "spostfix_high", "spostfix_mid", "spostfix_low")

stat_temp <- c(stat_list)
stat_power_df <- data.frame(matrix(ncol=6, nrow=0))
colnames(stat_power_df) <- c("set", "pop", "stat", "site_class", "maxp", "argmax")
for (set in subsets) {
  for (pop in pop_list) {
    for (stat in stat_temp) {
      for (site_class in c("neutral", "hardsweep", "softsweep")) {
        stat_power_df[nrow(stat_power_df)+1,] = list(set, pop, stat, site_class, 
         stat_power[[pop]][[stat]][[set]][[site_class]]$maxp, stat_power[[pop]][[stat]][[set]][[site_class]]$argmax)
      }
    }
  }
}
stat_power_df <- tbl_df(stat_power_df) %>% mutate_at(vars(pop, stat, site_class, argmax), funs(as.factor))

stat_temp <- c(stat_list)
stat_kstest_df <- data.frame(matrix(ncol=6, nrow=0))
colnames(stat_kstest_df) <- c("set", "pop", "stat", "site_class", "D", "nlog10p")
for (set in subsets) {
  for (pop in pop_list) {
    for (stat in stat_temp) {
      for (site_class in c("neutral_hardsweep", "neutral_softsweep", "hardsweep_softsweep")) {
        stat_kstest_df[nrow(stat_kstest_df)+1,] = list(set, pop, stat, site_class,
          stat_kstest[[pop]][[stat]][[set]][[site_class]]$statistic, stat_kstest[[pop]][[stat]][[set]][[site_class]]$pvalue)
      }
    }
  }
}
stat_kstest_df <- tbl_df(stat_kstest_df) %>% mutate_at(vars(pop, stat, site_class), funs(as.factor))

stat_temp <- c(stat_misc, stat_list)
stat_ksapp_df <- data.frame(matrix(ncol=5, nrow=0))
colnames(stat_ksapp_df) <- c("pop", "stat", "site_class", "D", "nlog10p")
for (pop in pop_list) {
  for (stat in stat_temp) {
    for (site_class in c("neutral", "hardsweep", "softsweep")) {
      stat_ksapp_df[nrow(stat_ksapp_df)+1,] = list(pop, stat, site_class,
       stat_ksapp[[pop]][[stat]][[site_class]]$statistic, 
       stat_ksapp[[pop]][[stat]][[site_class]]$pvalue)
    }
  }
}
stat_ksapp_df <- tbl_df(stat_ksapp_df) %>% mutate_at(vars(pop, stat), funs(as.factor))

saveRDS(stat_power_df, paste("../../results/feature_plots_data/", "stat_power_df.rds", sep=""))
saveRDS(stat_kstest_df, paste("../../results/feature_plots_data/", "stat_kstest_df.rds", sep=""))
saveRDS(stat_ksapp_df, paste("../../results/feature_plots_data/", "stat_ksapp_df.rds", sep=""))

# # # plot statistical power ranks
# # # & plot kolmogovor-smirnov ranks

stat_temp <- c(stat_list)
for (c in c("hardsweep", "softsweep")) {
  for (s in subsets) {
    ggplot(data=filter(stat_power_df, site_class==c & set==s)) + 
      geom_bar(mapping=aes(x=reorder(stat, maxp), y=maxp), fill=color_code[c], stat="identity", position="dodge") + 
      facet_wrap(~pop, nrow=1) + theme(aspect.ratio=1) + 
      geom_text(mapping=aes(x=stat, y=maxp, label=round(maxp, 2), hjust=-0.1), angle=90) + 
      xlab("") + ylab("Power") + ylim(c(0, 1.2)) + # ggtitle("?!") + 
      scale_color_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) + 
      scale_fill_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) + 
      scale_x_discrete(labels=stat_code[stat_temp]) + 
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + theme(panel.grid.major = element_blank())
    ggsave(file=paste("../../results/feature_plots_power/", c, "_", s, "_power", plout, sep=""))
  }
}

stat_temp <- c(stat_list)
color_temp <- c("hardsweep_softsweep"="sweep", "neutral_hardsweep"="hardsweep", "neutral_softsweep"="softsweep")
for (c in c("hardsweep_softsweep", "neutral_hardsweep", "neutral_softsweep")) {
  for (s in subsets) {
    ggplot(data=filter(stat_kstest_df, site_class==c & set==s)) + 
      geom_bar(mapping=aes(x=reorder(stat, D), y=D), fill=color_code[color_temp[c]], stat="identity", position="dodge") + 
      facet_wrap(~pop, nrow=1) + theme(aspect.ratio=1) + 
      geom_text(mapping=aes(x=stat, y=D, label=round(D, 2), hjust=-0.1), angle=90) + 
      xlab("") + ylab("K-S test statistic") + ylim(c(0, 1.2)) + # ggtitle("?!") + 
      scale_color_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) +
      scale_fill_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) +
      scale_x_discrete(labels=stat_code[stat_temp]) + 
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + theme(panel.grid.major = element_blank())
    ggsave(file=paste("../../results/feature_plots_power/", c, "_", s, "_kstest", plout, sep=""))
  }
}

stat_temp <- c(stat_misc, stat_list)
for (c in c("neutral")) {
  ggplot(data=filter(stat_ksapp_df, site_class==c)) + 
    geom_bar(mapping=aes(x=reorder(stat, D), y=D), fill=color_code[c], stat="identity", position="dodge") + 
    facet_wrap(~pop, nrow=1) + theme(aspect.ratio=1) + 
    geom_text(mapping=aes(x=stat, y=D, label=round(D, 2), hjust=-0.1), angle=90) + 
    xlab("") + ylab("K-S test statistic") + ylim(c(0, 1.2)) + # ggtitle("?!") + 
    scale_color_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) +
    scale_fill_manual(name="sweep class", breaks=names(color_code), labels=class_code, values=color_code) +
    scale_x_discrete(labels=stat_code[stat_temp]) + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) + theme(panel.grid.major = element_blank())
  ggsave(file=paste("../../results/feature_plots_power/", c, "_ksapp", plout, sep=""))
}
