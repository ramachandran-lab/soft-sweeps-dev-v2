library(ggbio)
library(caret)
library(pROC)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

source("functions_helper.R")
source("functions_classifier.R")
set.seed(2017070901)

my_table_app_main <- NULL
my_table_sim_train <- NULL
my_table_sim_test <- NULL
my_table_sim_link <- NULL
my_table_app_main_train <- NULL
for (pop in pop_list){
  # train dataset
  my_table_sim_train[[pop]] <- load_data(my_files_sim_train, pop)
  my_table_sim_train[[pop]] <- preprocess_data(my_table_sim_train[[pop]])
  # my_table_sim_train[[pop]] <- filter(my_table_sim_train[[pop]], is.na(complete) | complete==0)
  
  # test dataset
  my_table_sim_test[[pop]] <- load_data(my_files_sim_test, pop)
  my_table_sim_test[[pop]] <- preprocess_data(my_table_sim_test[[pop]])
  # my_table_sim_test[[pop]] <- filter(my_table_sim_test[[pop]], is.na(complete) | complete==0)
  
  # linked dataset
  my_table_sim_link[[pop]] <- load_data(my_files_sim_link, pop)
  my_table_sim_link[[pop]] <- preprocess_data(my_table_sim_link[[pop]])
  # my_table_sim_link[[pop]] <- filter(my_table_sim_link[[pop]], is.na(complete) | complete==0)
}

# classifier
simulation_helper <- function(prob_prior, stat_temp, stat_sweep_final, label) {
  # # # classifier
  for (pop in pop_list){
    # datasets
    train_x <- my_table_sim_train[[pop]][, stat_temp]
    train_y <- my_table_sim_train[[pop]][, "site_class"]
    test_x <- my_table_sim_test[[pop]][, stat_temp]
    test_y <- my_table_sim_test[[pop]][, "site_class"]
    link_x <- my_table_sim_link[[pop]][, stat_temp]
    link_y <- my_table_sim_link[[pop]][, "stype"]
    
    # # # simulation
    prob_postodds <- 10^1
    fit_classifier_uni <- fit_mix_uni(train_x, train_y)  # , 0.995)
    test_y_nb <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_uni, test_x), stat_sweep_final), prob_prior, prob_postodds)
    link_y_nb <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_uni, link_x), stat_sweep_final), prob_prior, prob_postodds)
    
    # posterior biplots
    temp <- as_tibble(cbind(my_table_sim_test[[pop]], test_y_nb$main))
    ggplot() + 
      geom_point(data=temp, mapping=aes(y=logit((1-neutral)), x=logit(softsweep/(softsweep+hardsweep)), color=pred_class_alt), alpha=0.1) +
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) +
      geom_hline(yintercept=logit(1/2), color="grey") + geom_vline(xintercept=log10(prob_postodds), color="grey") +
      geom_vline(xintercept=logit(1/2), color="grey") + geom_vline(xintercept=-log10(prob_postodds), color="grey") +
      facet_wrap(~site_class, nrow=1) + theme(aspect.ratio=1) + guides(color=guide_legend(override.aes=list(alpha=1)))
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_logit_posterior/", 
      pop, "-biplot_logit_posterior_class-", label, plout, sep=""))
    ggplot() + 
      geom_density(data=filter(temp, T), mapping=aes(x=logit(softsweep/(softsweep+hardsweep)), color=site_class)) +
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) +
      geom_vline(xintercept=logit(1/2), color="grey") + geom_vline(xintercept=log10(1/prob_postodds), color="grey") +
      geom_vline(xintercept=log10(prob_postodds), color="grey") + theme(aspect.ratio=1) +
      facet_wrap(~site_class, nrow=1) + xlab("logit(softsweep_hardsweep.posterior)")
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_logit_posterior/", 
      pop, "-uniplot_logit_posterior_class-", label, plout, sep=""))
    ggplot() + 
      geom_point(data=temp, mapping=aes(y=logit(1-neutral), x=logit(softsweep/(softsweep+hardsweep)), color=site_class), alpha=0.05) +
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) +
      geom_hline(yintercept=logit(1/2), color="grey") + geom_vline(xintercept=log10(prob_postodds), color="grey") + 
      geom_vline(xintercept=logit(1/2), color="grey") + geom_vline(xintercept=-log10(prob_postodds), color="grey") +
      xlab("logit(softsweep_hardsweep_posterior") + ylab("logit(sweep_neutral_posterior") +
      theme(aspect.ratio=1) + guides(color=guide_legend(override.aes=list(alpha=1))) +
      xlab("logit(softsweep_hardsweep.posterior)") + ylab("logit(sweep_neutral.posterior)")
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_logit_posterior/", 
      pop, "-biplot_logit_posterior_all-", label, plout, sep=""))
    
    # confusion matrix linked
    temp <- as_tibble(cbind(my_table_sim_link[[pop]], link_y_nb$main))
    conf_mat <- as_tibble(table(temp[c("stype", "pred_class_alt", "dist_physpos")])) %>%
      mutate(pred_class_alt=as.factor(pred_class_alt)) %>% mutate(stype=as.factor(stype))
    conf_mat <- conf_mat %>% group_by(stype, dist_physpos) %>% mutate(m = n/sum(n))
    ggplot() + geom_bar(data=conf_mat, mapping=aes(x=as.factor(as.numeric(dist_physpos)), y=m, fill=pred_class_alt), stat="identity") +
      scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) +
      theme(panel.grid.major = element_blank()) + facet_wrap(~stype, nrow=3) +
      theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5))
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
      pop, "-confusion_prior_class_pop_linked-", label, plout, sep=""))
    write_tsv(conf_mat, paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
      pop, "-confusion_prior_class_pop_linked-", label, ".txt", sep=""))
    
    # confusion matrix window
    temp <- as_tibble(cbind(my_table_sim_test[[pop]], test_y_nb$main))
    z <- temp %>% mutate(
        # sweep_neutral_value = logit(neutral/(neutral+softsweep+hardsweep)),
        sweep_neutral_value = -log10(neutral),  # -log neutral posterior

        # 
        # # 
        # 

        # softsweep_hardsweep_value = logit(softsweep/(softsweep+hardsweep)))
        softsweep_hardsweep_value = log10(softsweep/hardsweep))  # log posterior odds of soft vs hard

        # # 
        # # 
        # # 

    w <- as_tibble(table(z[c("stype", "pred_class_alt")])); w <- w %>% group_by(stype) %>% mutate(m = n/sum(n))
    ggplot() + geom_bar(data=w, mapping=aes(x=stype, y=m, fill=pred_class_alt), stat="identity") + 
      scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      theme(panel.grid.major = element_blank()) + theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5))
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
                 pop, "-confusion_prior_class_pop_window-", label, plout, sep=""))
    temp_window <- w
    write_tsv(w, paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
                       pop, "-confusion_prior_class_pop_window-", label, ".txt", sep=""))
    
    # confusion matrix region
    temp <- as_tibble(cbind(my_table_sim_link[[pop]], link_y_nb$main))
    x <- temp %>% filter(pred_class!="neutral") %>% group_by(key) %>% filter(n() >= 3) %>% 
      dplyr::summarise(stype = stype[1], 
        # sweep_neutral_value = mean(logit(neutral/(neutral+softsweep+hardsweep))),
        sweep_neutral_value = mean(-log10(neutral)),  # -log neutral posterior

        # 
        # # 
        # 

        # softsweep_hardsweep_value = mean((-log10(neutral))*logit(softsweep/(softsweep+hardsweep)))/mean((-log10(neutral)))) %>% 
        softsweep_hardsweep_value = mean((-log10(neutral))*(log10(softsweep/(hardsweep))))/mean((-log10(neutral)))) %>%   # log posterior odds of soft vs hard

        # 
        # # 
        # 

        # # softsweep_hardsweep_value = mean(logit(softsweep/(softsweep+hardsweep)))) %>% 
      # mutate(pred_class=ifelse(softsweep_hardsweep_value>logit(1/(1+prob_postodds)), ifelse(softsweep_hardsweep_value>logit(0.5), 
      #   ifelse(softsweep_hardsweep_value>logit(prob_postodds/(1+prob_postodds)), "softsweep", "softsweepweak"), "hardsweepweak"), "hardsweep"))
      mutate(pred_class=ifelse(softsweep_hardsweep_value>log10(1/prob_postodds), ifelse(softsweep_hardsweep_value>log10(1),  # log posterior odds of soft vs hard
        ifelse(softsweep_hardsweep_value>log10(prob_postodds), "softsweep", "softsweepweak"), "hardsweepweak"), "hardsweep"))

      # # 
      # # 
      # # 

    y <- temp %>% filter(!(key %in% unique(x$key))) %>% group_by(key) %>% dplyr::summarise(stype=stype[1],
        sweep_neutral_value=NA, softsweep_hardsweep_value=NA, pred_class="neutral")
    z <- rbind(x, y)
    w <- as_tibble(table(z[c("stype", "pred_class")])); w <- w %>% group_by(stype) %>% mutate(m = n/sum(n))
    w$pred_class <- factor(as.factor(w$pred_class), levels(as.factor(w$pred_class))[c(1,2,3,5,4)])
    table(z[c("stype", "pred_class")])/rowSums(table(z[c("stype", "pred_class")]))
    ggplot() + geom_bar(data=w, mapping=aes(x=stype, y=m, fill=pred_class), stat="identity") + 
      scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      theme(panel.grid.major = element_blank()) + theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5))
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
      pop, "-confusion_prior_class_pop_region-", label, plout, sep=""))
    write_tsv(w, paste("../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/",
      pop, "-confusion_prior_class_pop_region-", label, ".txt", sep=""))
  }
} 

# params
stat_temp <- stat_list
stat_sweep_final <- 
  c("XPEHH_win", "FST_win", "D_theta_Pi_win", "H12_win", "H2H1_win", "Taj_D_win", "H2_win", "ZA_win", "iHS_win")
  # c("H12_win", "XPEHH_win", "H2_win", "ZA_win", "iHS_win", "H2H1_win", "FST_win", "D_theta_Pi_win")

prob_neutral <- 0.98
prob_prior <-
  c("neutral"=prob_neutral,
    "hardsweep"=(1-prob_neutral)/2,
    "softsweep"=(1-prob_neutral)/2)
label <- "0.98-simulated"
simulation_helper(prob_prior, stat_temp, stat_sweep_final, label)

prob_neutral <- 0.90
prob_prior <-
  c("neutral"=prob_neutral,
    "hardsweep"=(1-prob_neutral)/2,
    "softsweep"=(1-prob_neutral)/2)
label <- "0.90-simulated"
simulation_helper(prob_prior, stat_temp, stat_sweep_final, label)

prob_neutral <- 0.34
prob_prior <-
  c("neutral"=prob_neutral,
    "hardsweep"=(1-prob_neutral)/2,
    "softsweep"=(1-prob_neutral)/2)
label <- "0.34-simulated"
simulation_helper(prob_prior, stat_temp, stat_sweep_final, label)

# # # 
confusion_window <- NULL
file_dir <- "../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/"
confusion_window$lowest_window$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_window-0.98-simulated.txt"))
confusion_window$lower_window$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_window-0.90-simulated.txt"))
confusion_window$low_window$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_window-0.34-simulated.txt"))
confusion_window$lowest_window$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_window-0.98-simulated.txt"))
confusion_window$lower_window$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_window-0.90-simulated.txt"))
confusion_window$low_window$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_window-0.34-simulated.txt"))
confusion_window$lowest_window$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_window-0.98-simulated.txt"))
confusion_window$lower_window$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_window-0.90-simulated.txt"))
confusion_window$low_window$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_window-0.34-simulated.txt"))

prior_map <- c("low_window"=0.34, "lower_window"=0.9, "lowest_window"=0.98)
for (i in names(prior_map)) {
  for (j in pop_list) {
    confusion_window[[i]][[j]] <- confusion_window[[i]][[j]] %>% mutate(prior=prior_map[[i]], population=j)
  }
}
confusion_window_final <- NULL
for (i in names(prior_map)) {
  for (j in pop_list) {
    if (is.null(confusion_window_final)) {
      confusion_window_final <- confusion_window[[i]][[j]]
    } else {
      confusion_window_final <- rbind(confusion_window_final, confusion_window[[i]][[j]])
    }
  }
}
confusion_window_final$population <- factor(confusion_window_final$population, pop_list)
confusion_window_final$pred_class_alt <-
  factor(confusion_window_final$pred_class_alt, c("softsweep", "softsweepweak", "neutral", "hardsweepweak", "hardsweep"))
confusion_window_final <- dplyr::rename(confusion_window_final, pred_class=pred_class_alt)

ggplot(data=confusion_window_final) + geom_bar(mapping=aes(x=stype, y=m, fill=pred_class), stat="identity") + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  facet_wrap(~population+prior) + theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + 
  xlab("Sweep simulation") + ylab("Proportion")
ggsave(paste("../../results/classifier_plots_miscellaneous/confusion_prior_class_pop_window", plout, sep=""))
write_tsv(confusion_window_final, paste("../../results/classifier_plots_miscellaneous/confusion_prior_class_pop_window", ".tsv", sep=""))

# # # 
confusion_region <- NULL
file_dir <- "../../results/classifier_plots_miscellaneous/classifier_plots_confusion_matrix/"
confusion_region$lowest_region$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_region-0.98-simulated.txt"))
confusion_region$lower_region$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_region-0.90-simulated.txt"))
confusion_region$low_region$YRI <- read_tsv(paste(file_dir, "", sep="YRI-confusion_prior_class_pop_region-0.34-simulated.txt"))
confusion_region$lowest_region$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_region-0.98-simulated.txt"))
confusion_region$lower_region$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_region-0.90-simulated.txt"))
confusion_region$low_region$CEU <- read_tsv(paste(file_dir, "", sep="CEU-confusion_prior_class_pop_region-0.34-simulated.txt"))
confusion_region$lowest_region$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_region-0.98-simulated.txt"))
confusion_region$lower_region$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_region-0.90-simulated.txt"))
confusion_region$low_region$CHB <- read_tsv(paste(file_dir, "", sep="CHB-confusion_prior_class_pop_region-0.34-simulated.txt"))

prior_map <- c("low_region"=0.34, "lower_region"=0.9, "lowest_region"=0.98)
for (i in names(prior_map)) {
  for (j in pop_list) {
    confusion_region[[i]][[j]] <- confusion_region[[i]][[j]] %>% mutate(prior=prior_map[[i]], population=j)
  }
}
confusion_region_final <- NULL
for (i in names(prior_map)) {
  for (j in pop_list) {
    if (is.null(confusion_region_final)) {
      confusion_region_final <- confusion_region[[i]][[j]]
    } else {
      confusion_region_final <- rbind(confusion_region_final, confusion_region[[i]][[j]])
    }
  }
}
confusion_region_final$population <- factor(confusion_region_final$population, pop_list)
confusion_region_final$pred_class <-
  factor(confusion_region_final$pred_class, c("softsweep", "softsweepweak", "neutral", "hardsweepweak", "hardsweep"))

ggplot(data=confusion_region_final) + geom_bar(mapping=aes(x=stype, y=m, fill=pred_class), stat="identity") + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  facet_wrap(~population+prior) + theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + 
  xlab("Sweep simulation") + ylab("Proportion")
ggsave(paste("../../results/classifier_plots_miscellaneous/confusion_prior_class_pop_region", plout, sep=""))
write_tsv(confusion_region_final, paste("../../results/classifier_plots_miscellaneous/confusion_prior_class_pop_region", ".tsv", sep=""))


confusion_window_final$type <- "window"
confusion_region_final$type <- "region"
confusion_final <- rbind(confusion_window_final, confusion_region_final)

ggplot(data=confusion_final) + geom_bar(mapping=aes(x=stype, y=m, fill=pred_class), stat="identity") + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  facet_wrap(~population+prior) + theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + 
  xlab("Sweep simulation") + ylab("Proportion")

# # # 

# # # simulation, condition on sweep
# prob_postodds <- 10^1.5
# fit_classifier_uni <- fit_mix_uni(train_x, train_y)  # , 0.995)
# train_y_a <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_uni, train_x), stat_sweep_final), prob_prior, prob_postodds)
# test_y_a <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_uni, test_x), stat_sweep_final), prob_prior, prob_postodds)
# link_y_a <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_uni, link_x), stat_sweep_final), prob_prior, prob_postodds)
# train_x_b <- as_tibble(cbind(my_table_sim_train[[pop]], train_y_a$main)) %>% filter(pred_class!="neutral") %>% .[, stat_temp]
# train_y_b <- as_tibble(cbind(my_table_sim_train[[pop]], train_y_a$main)) %>% filter(pred_class!="neutral") %>% .[, "site_class"]
# fit_classifier_temp <- fit_mix_uni(train_x_b, train_y_b)  # , 0.995)
# test_y_b <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_temp, test_x), stat_sweep_final), 
#   c("neutral"=1/3, "hardsweep"=1/3, "softsweep"=1/3), prob_postodds)
# link_y_b <- pred_mix_nb(logd_mix_nb(logd_mix_uni(fit_classifier_temp, link_x), stat_sweep_final), 
#   c("neutral"=1/3, "hardsweep"=1/3, "softsweep"=1/3), prob_postodds)
# 
# test_y_nb <- test_y_a$main
# test_y_nb <- test_y_b$main
# test_y_nb$neutral <- test_y_a$main$neutral
# test_y_nb$pred_class_alt[which(test_y_a$main$pred_class_alt=="neutral")] <- "neutral"
# test_y_nb$pred_class_alt[which(test_y_a$main$pred_class_alt!="neutral")] <- 
#   test_y_b$main$pred_soft_hard_alt[which(test_y_a$main$pred_class_alt!="neutral")]
# 
# link_y_nb <- link_y_a$main
# link_y_nb <- link_y_b$main
# link_y_nb$neutral <- link_y_a$main$neutral
# link_y_nb$hardsweep <- link_y_a$main$hardsweep
# link_y_nb$softsweep <- link_y_a$main$softsweep
# link_y_nb$pred_class_alt[which(link_y_a$main$pred_class_alt=="neutral")] <- "neutral"
# link_y_nb$pred_class_alt[which(link_y_a$main$pred_class_alt!="neutral")] <- 
#   link_y_b$main$pred_soft_hard_alt[which(link_y_a$main$pred_class_alt!="neutral")]

# # # # 
# temp <- filter(as_tibble(cbind(my_table_sim_link[[pop]], link_y_nb$main)), is.na(complete)|complete==1)
# temp <- filter(as_tibble(cbind(my_table_sim_link[[pop]], link_y_nb$main)), is.na(sweep_freq_pop)|sweep_freq_pop<=0.8)
