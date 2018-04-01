source("functions_helper.R")
set.seed(2017070901)

# # # 

# # # statistics

# plots_uni_dens <- NULL
# plots_bi_dens <- NULL
# plots_taj_dens <- NULL

my_table_app_main <- NULL
my_table_sim_train <- NULL
for (pop in pop_list){
  # app load
  my_table_app_main[[pop]] <- load_data(my_files_app_main, pop)
  my_table_app_main[[pop]] <- preprocess_data(my_table_app_main[[pop]])
  
  # train load
  my_table_sim_train[[pop]] <- load_data(my_files_sim_train, pop)
  my_table_sim_train[[pop]] <- preprocess_data(my_table_sim_train[[pop]])
}
my_table_app_main <- preprocess_data(bind_rows(my_table_app_main[pop_list]))
my_table_sim_train <- preprocess_data(bind_rows(my_table_sim_train[pop_list]))

# # # metadata statistics

# metadata univariate
meta_temp <- c("stime", "scoeff", "sfreq", "complete", "sweep_freq_bin",
               "sweep_freq_pop", "sweep_fix_time", "sweep_dur_time", "origins")
for (meta in meta_temp) {
  if (!is.factor(my_table_sim_train[[meta]])) {
    let(list(COL=meta), ggplot() + 
      # geom_rug(data=my_table_sim_train, mapping=aes(x=COL, group=site_class, color=site_class), alpha=0.5) + 
      # guides(color=F, fill=guide_legend(override.aes=list(alpha=1))) + 
      geom_density(data=my_table_sim_train, mapping=aes(x=COL, group=site_class, color=site_class, fill=site_class), alpha=0.2) + 
      xlab(stat_code[[meta]])) + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) +
      scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank())
  }
  ggsave(file=paste("../../results/summary_plots_uni_meta/", meta, "_density", plout, sep=""))
}

# metadata bivariate
meta_temp <- c("stime", "scoeff", "sfreq", "complete", "sweep_freq_bin",
               "sweep_freq_pop", "sweep_fix_time", "sweep_dur_time", "origins")
for (i in 1:(length(meta_temp)-1)) {
  for (j in (i+1):length(meta_temp)) {
    meta1 = meta_temp[i]; meta2 = meta_temp[j]
    if (!is.factor(my_table_sim_train[[meta1]]) & !is.factor(my_table_sim_train[[meta2]])) {
      let(list(COL1=meta1, COL2=meta2), ggplot() + 
        # geom_density2d(data=my_table_sim_train, mapping=aes(x=COL1, y=COL2, group=site_class, color=site_class)) + 
        geom_point(data=my_table_sim_train, mapping=aes(x=COL1, y=COL2, group=site_class, color=site_class), alpha=0.5, size=0.8) + 
        xlab(stat_code[[meta1]]) + ylab(stat_code[[meta2]]) + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank()))
    } else if (is.factor(my_table_sim_train[[meta1]])) {
      let(list(COL1=meta1, COL2=meta2), ggplot() + 
        geom_boxplot(data=my_table_sim_train, mapping=aes(x=COL1, y=COL2, fill=site_class, color=site_class), alpha=0.5, outlier.size=0.8) + 
        xlab(stat_code[[meta1]]) + ylab(stat_code[[meta2]]) + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank()) + 
        theme(panel.grid.major = element_blank()))
    } else if (is.factor(my_table_sim_train[[meta2]])) {
      let(list(COL1=meta1, COL2=meta2), ggplot() + 
        geom_boxplot(data=my_table_sim_train, mapping=aes(x=COL2, y=COL1, fill=site_class, color=site_class), alpha=0.5, outlier.size=0.8) + 
        xlab(stat_code[[meta1]]) + ylab(stat_code[[meta2]]) + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank()) + 
        theme(panel.grid.major = element_blank()))
    }
    ggsave(file=paste("../../results/summary_plots_bi_meta/", meta1, "_", meta2, "_density", plout, sep=""))
  }
}

# # # summary statistics

# univariate density
stat_temp <- c(stat_misc, stat_list)
for (stat in stat_temp) {
  print(c(stat))
  # # reduce plot size cheat
  # #  by subsampling dense region
  # rug_temp <- let(list(COL=stat), filter(my_table_app_main, !is.na(my_table_app_main$COL) & 
  #   (COL<quantile(my_table_app_main$COL, 0.05)|COL>quantile(my_table_app_main$COL, 0.95))))
  # rug_samp <- sample_frac(let(list(COL=stat), filter(my_table_app_main, !is.na(my_table_app_main$COL) & 
  #   (!(COL<quantile(my_table_app_main$COL, 0.05)|COL>quantile(my_table_app_main$COL, 0.95))))), 0.2)
  # rug_temp <- bind_rows(rug_temp, rug_samp)
  let(list(COL=stat), ggplot() + 
    # geom_rug(data=rug_temp, mapping=aes(x=COL), sides="b", alpha=1, size=0.5) + 
    geom_density(data=sample_n(my_table_sim_train, min(10000, nrow(my_table_sim_train))), 
      mapping=aes(x=COL, group=site_class, color=site_class, fill=site_class), alpha=0.2) +
    geom_density(data=sample_n(my_table_app_main, min(10000, nrow(my_table_app_main))), 
      mapping=aes(x=COL), color="black", fill=NA, linetype="dashed") +
    # geom_bkde(data=my_table_sim_train, mapping=aes(x=COL, group=site_class, color=site_class, fill=site_class), alpha=0.2, gridsize=100) +
    # geom_bkde(data=my_table_app_main, mapping=aes(x=COL), color="black", fill=NA, linetype="dashed", gridsize=100) + 
    xlab(stat_code[[stat]])) + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) + # ggtitle(pop) + 
    scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
    scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
    theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank())
  ggsave(file=paste("../../results/summary_plots_uni_summ/", stat, "_density", plout, sep=""))
}

# bivariate density
stat_temp <- c(stat_misc, stat_list)
for (i in 1:(length(stat_temp)-1)) {
  for (j in (i+1):length(stat_temp)) {
    stat1 = stat_temp[i]; stat2 = stat_temp[j]
    print(c(stat1, stat2))
    let(list(COL1=stat1, COL2=stat2), ggplot() + 
      # geom_hex(data=my_table_sim_train, mapping=aes(x=COL1, y=COL2, fill=site_class), bins=100, show.legend=F, alpha=0.2) + 
      geom_hex(data=my_table_app_main, mapping=aes(x=COL1, y=COL2), color="grey", fill="grey", bins=75, show.legend=F, alpha=1) + 
      # geom_bkde2d(data=my_table_sim_train, mapping=aes(x=COL1, y=COL2, group=site_class, color=site_class), grid_size=c(30,30)) +
      # geom_bkde2d(data=my_table_app_main, mapping=aes(x=COL1, y=COL2), color="black", linetype="dashed", grid_size=c(30,30)) + 
      geom_density2d(data=sample_n(my_table_sim_train, min(10000, nrow(my_table_sim_train))), 
        mapping=aes(x=COL1, y=COL2, group=site_class, color=site_class)) +
      geom_density2d(data=sample_n(my_table_app_main, min(10000, nrow(my_table_app_main))), 
        mapping=aes(x=COL1, y=COL2), color="black", linetype="dashed") +
      xlab(stat_code[[stat1]]) + ylab(stat_code[[stat2]]) + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      # scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.title=element_blank()))
    ggsave(file=paste("../../results/summary_plots_bi_summ/", stat1, "_", stat2, "_density", plout, sep=""))
  }
}

my_table_app_main <- NULL
my_table_sim_train <- NULL

# # # summary at linked sites

my_table_sim_link <- NULL
for (pop in pop_list){
  # link load
  my_table_sim_link[[pop]] <- load_data(my_files_sim_link, pop)
  my_table_sim_link[[pop]] <- preprocess_data(my_table_sim_link[[pop]])
}
my_table_sim_link <- preprocess_data(bind_rows(my_table_sim_link[pop_list]))

# univariate density
stat_temp <- c(stat_misc, stat_list)
for (stat in stat_temp) {
  print(c(stat))
  options(scipen=999)
  let(list(COL=stat), ggplot() + 
    geom_smooth(data=my_table_sim_link, 
       mapping=aes(x=as.numeric(dist_physpos/1000), y=COL, color=site_class, fill=site_class), level=0) + 
    # geom_smooth(data=filter(my_table_sim_link, complete==1 | site_class=="neutral"), 
    #   mapping=aes(x=as.numeric(dist_physpos/1000), y=COL, color=site_class, fill=site_class), level=0) + 
    # geom_smooth(data=filter(my_table_sim_link, complete==0 | site_class=="neutral"), 
    #   mapping=aes(x=as.numeric(dist_physpos/1000), y=COL, color=site_class, fill=site_class), level=0, linetype="dashed") + 
    xlab(stat_code[["dist_physpos"]]) + ylab(stat_code[[stat]]) + facet_wrap(~pop) + theme(aspect.ratio=1) +  # ggtitle(pop) + 
    # geom_boxplot(data=my_table_sim_link, mapping=aes(x=as.factor(dist_physpos/1000), y=COL, fill=site_class), 
    # outlier.size=0.1, outlier.colour=NULL) + 
    # xlab(stat_code[["dist_physpos"]]) + ylab(stat_code[[stat]]) + facet_wrap(~pop, nrow=3) + theme(aspect.ratio=1/3) + ggtitle(pop) + 
    scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
    scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
    theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + 
    theme(legend.title=element_blank()))  #  + theme(panel.grid.major = element_blank()))
  ggsave(file=paste("../../results/summary_plots_uni_linked/", stat, "_linked", plout, sep=""))
  options(scipen=0)
}
my_table_sim_link <- NULL

# # # 

# # count iHS and XPEHH
# for (pop in pop_list) {
#   # count -inf values for iHS and XPEHH
#   count_app_main <- NULL
#   count_app_main$total <- count(my_table_app_main[[pop]])
#   count_app_main$iHS <- count(filter(my_table_app_main[[pop]], is.infinite(iHS_win)))
#   count_app_main$XPEHH <- count(filter(my_table_app_main[[pop]], is.infinite(XPEHH_win)))
#   count_app_main$iHS_XPEHH <- count(filter(my_table_app_main[[pop]], is.infinite(iHS_win) & is.infinite(XPEHH_win)))
#   
#   count_sim_train <- NULL
#   count_sim_train$total <- count(my_table_sim_train[[pop]])
#   count_sim_train$site_class <- count(my_table_sim_train[[pop]], site_class)
#   count_sim_train$iHS <- count(filter(my_table_sim_train[[pop]], is.infinite(iHS_win)), site_class)
#   count_sim_train$XPEHH <- count(filter(my_table_sim_train[[pop]], is.infinite(XPEHH_win)), site_class)
#   count_sim_train$iHS_XPEHH <- count(filter(my_table_sim_train[[pop]], is.infinite(iHS_win) & is.infinite(XPEHH_win)), site_class)
#   
#   saveRDS(count_app_main, paste("../../results/summary_plots_count/", pop, "-", "iHS_win", "_", "XPEHH_win", "_count_app_main.rds", sep=""))
#   saveRDS(count_sim_train, paste("../../results/summary_plots_count/", pop, "-", "iHS_win", "_", "XPEHH_win", "_count_sim_train.rds", sep=""))
# }
