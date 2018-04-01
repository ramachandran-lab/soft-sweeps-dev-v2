library(ggbio)
# library(caret)
# library(pROC)
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
  
  # application dataset
  my_table_app_main[[pop]] <- load_data(my_files_app_main, pop)
  my_table_app_main[[pop]] <- preprocess_data(my_table_app_main[[pop]])
  # my_table_app_main[[pop]] <- filter(my_table_app_main[[pop]], is.na(complete) | complete==0)
  
  # application dataset
  # (empirical neutral)
  my_table_app_main_train[[pop]] <- my_table_sim_train[[pop]]
  neutral_train_size <- nrow(filter(my_table_sim_train[[pop]], site_class=="neutral"))
  my_table_new_neutral <- my_table_app_main[[pop]]  # [sample(nrow(my_table_app_main[[pop]]), neutral_train_size),]
  my_table_new_neutral <- my_table_new_neutral %>% mutate(stype="neutral") %>% mutate(site_class="neutral")
  my_table_app_main_train[[pop]] <- filter(my_table_app_main_train[[pop]], site_class!="neutral")
  my_table_app_main_train[[pop]] <- rbind(my_table_app_main_train[[pop]], my_table_new_neutral)
}

# # # application
application_helper <- function(prob_prior, stat_temp, stat_sweep_final, label, simulated) {
  for (pop in pop_list){
    # datasets
    train_x <- my_table_sim_train[[pop]][, stat_temp]
    train_y <- my_table_sim_train[[pop]][, "site_class"]
    test_x <- my_table_sim_test[[pop]][, stat_temp]
    test_y <- my_table_sim_test[[pop]][, "site_class"]
    link_x <- my_table_sim_link[[pop]][, stat_temp]
    link_y <- my_table_sim_link[[pop]][, "stype"]
    train_app_c <- my_table_app_main_train[[pop]][, "chr"]
    train_app_x <- my_table_app_main_train[[pop]][, stat_temp]
    train_app_y <- my_table_app_main_train[[pop]][, "site_class"]
    appl_c <- my_table_app_main[[pop]][, "chr"]
    appl_x <- my_table_app_main[[pop]][, stat_temp]
    
    # # # application
    prob_postodds <- 10^1
    if (simulated) {
      fit_classifier_appl <- fit_mix_uni(train_x, train_y)  # , 0.995)
      appl_y_nb <- pred_mix_nb(logd_mix_nb(logd_mix_uni(
        fit_classifier_appl, appl_x), stat_sweep_final), prob_prior, prob_postodds)
    } else {
      # by genome
      fit_classifier_appl <- fit_mix_uni(train_app_x, train_app_y)  # , 0.995)
      appl_y_nb <- pred_mix_nb(logd_mix_nb(logd_mix_uni(
        fit_classifier_appl, appl_x), stat_sweep_final), prob_prior, prob_postodds)
      # # by chromosome
      # fit_classifier_appl_chr <- NULL
      # appl_y_nb_chr <- NULL
      # for (c in chr_list) {
      #   fit_classifier_appl_chr[[c]] <-
      #     fit_mix_uni(train_app_x[which(is.na(train_app_c) | train_app_c==c),],
      #                 train_app_y[which(is.na(train_app_c) | train_app_c==c),])  # , 0.995)
      #   appl_y_nb_chr[[c]] <- pred_mix_nb(logd_mix_nb(logd_mix_uni(
      #     fit_classifier_appl_chr[[c]], appl_x[which(is.na(appl_c) | appl_c==c),]),
      #     stat_sweep_final), prob_prior, prob_postodds)
      # }
      # appl_y_nb <- NULL
      # appl_y_nb$main <- appl_y_nb_chr[[1]]$main
      # for (c in chr_list[-1]) {
      #   appl_y_nb$main <- rbind(appl_y_nb$main, appl_y_nb_chr[[c]]$main)
      # }
    }
    temp <- as_tibble(cbind(my_table_app_main[[pop]], appl_y_nb$main))
    prediction_table_windows <- temp
    write_tsv(prediction_table_windows, paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/", 
      pop, "-prediction_table_windows-", label, ".tsv", sep=""))
    
    # karyogram for each population
    pred_temp <- 
      c("hardsweep", "neutral", "softsweep", "pred_class", "pred_class_alt", 
        "pred_sweep_neutral", "pred_sweep_neutral_alt", "pred_soft_hard", "pred_soft_hard_alt")
    # recover 50 Kb windows
    a <- GRanges(seqnames=paste("chr", temp$chr, sep=""), 
                 ranges=IRanges(start=temp$physpos-25000, end=temp$physpos+24999), mcols=temp[,pred_temp])
    a <- as_tibble(a); names(a) <- c("seqnames", "start", "end", "width", "strand", pred_temp); a <- GRanges(a)
    seqinfo(a) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(a))]
    
    # filter sweep windows
    b <- a[a$pred_class_alt != "neutral"]
    gr_karyogram_windows <- b; tb_karyogram_windows <- as_tibble(gr_karyogram_windows)
    autoplot(seqinfo(gr_karyogram_windows)) + 
      layout_karyogram(gr_karyogram_windows, aes(color=pred_class_alt, fill=pred_class_alt)) + 
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
    ggsave(paste("../../results/classifier_plots_karyograms/", pop, "-karyogram_sweep_windows-", label, plout, sep=""))
    write_tsv(tb_karyogram_windows, paste("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/", 
                                          pop, "-karyogram_sweep_windows-", label, ".tsv", sep=""))
    
    # merge windows within 0.5 Mb for empirical
    if (simulated==F) {
      max_distance <- 500000
    } else {  # within 0.1 Mb for simulated
      max_distance <- 500000
    }
    
    d <- GenomicRanges::reduce(extend(b, upstream=max_distance, downstream=max_distance))
    e <- subjectHits(findOverlaps(b, d)); b_df <- as_tibble(b); b_df$group <- e
    start(d) <- (b_df %>% group_by(group) %>% summarise(start=min(start)))$start
    end(d) <- (b_df %>% group_by(group) %>% summarise(end=max(end)))$end
    seqinfo(d) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(d))]
    
    # identify chrm ends and large gaps (5 Mb)
    j <- GRanges(seqnames=paste("chr", my_table_app_main[[pop]]$chr, sep=""), 
                 ranges=IRanges(start=my_table_app_main[[pop]]$physpos-25000, end=my_table_app_main[[pop]]$physpos+24999))
    seqinfo(j) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(j))]
    j <- GenomicRanges::reduce(j); j <- gaps(j) %>% .[strand(.)=="*"]; j <- j[width(j)>5000000]  # remove near large gaps (5Mb)
    m <- GRanges(seqnames=paste("chr", c(1:22), sep=""), ranges=IRanges(start=rep(1,22), end=rep(0,22)))
    seqinfo(m) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(m))]
    n <- GRanges(seqnames=paste("chr", c(1:22), sep=""), ranges=IRanges(start=seqlengths(j), end=seqlengths(j)-1))
    seqinfo(n) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(n))]
    
    # remove regions within 1 window of above
    j <- GenomicRanges::reduce(extend(c(j,m,n), upstream=0, downstream=0))
    # j <- GenomicRanges::reduce(extend(c(j,m,n), upstream=50000, downstream=50000))
    autoplot(seqinfo(j)) + layout_karyogram(j, color="blue", fill="blue")
    k <- queryHits(findOverlaps(d, j)); if (length(k) > 0) {d <- d[-k]}
    
    # calculate weighted mean logit scores
    f <- b_df %>% group_by(group) %>%
      # summarise(softsweep_hardsweep_value=mean(logit(softsweep/(softsweep+hardsweep))))
      summarise(softsweep_hardsweep_value=mean((-log10(neutral))*log10(softsweep/(hardsweep)))/mean((-log10(neutral))))  # log posterior odds of soft vs hard

      # 
      # # 
      # 

    f$sweep_neutral_value <- b_df %>% group_by(group) %>% # mean logit sweep_neutral
      # summarise(sweep_neutral_value=mean(logit(neutral/(neutral+softsweep+hardsweep)))) %>% .$sweep_neutral_value
      summarise(sweep_neutral_value=mean(-log10(neutral))) %>% .$sweep_neutral_value  # -log neutral posterior

      # 
      # # 
      # 

    # f <- f %>% mutate(pred_class=ifelse(softsweep_hardsweep_value>logit(1/(1+prob_postodds)), ifelse(softsweep_hardsweep_value>logit(0.5), 
    #   ifelse(softsweep_hardsweep_value>logit(prob_postodds/(1+prob_postodds)), "softsweep", "softsweepweak"), "hardsweepweak"), "hardsweep"))
    f <- f %>% mutate(pred_class=ifelse(softsweep_hardsweep_value>log10(1/(prob_postodds)), ifelse(softsweep_hardsweep_value>log10(1),  # log posterior odds of soft vs hard
      ifelse(softsweep_hardsweep_value>log10(prob_postodds), "softsweep", "softsweepweak"), "hardsweepweak"), "hardsweep"))

    # # 
    # # 
    # # 

    if (length(k) > 0) {f <- f[-k,]}
    d$pred_class <- f$pred_class
    d$softsweep_hardsweep_value <- f$softsweep_hardsweep_value
    d$sweep_neutral_value <- f$sweep_neutral_value
    d$group <- c(1:length(d))
    
    # filter sweep regions
    # drop regions with >= range width
    # d <- d[width(d) > 50000]
    # drop regions with >= sweep windows
    d <- as_tibble(d[queryHits(findOverlaps(d, b))]) %>% group_by(group) %>% filter(n() >= 3) %>% 
      ungroup() %>% select(-c(group)) %>% unique() %>% GRanges()
    seqinfo(d) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[levels(seqnames(d))]
    gr_karyogram_regions <- d; tb_karyogram_regions <- as_tibble(gr_karyogram_regions)
    if (length(d) > 0){
      autoplot(seqinfo(gr_karyogram_regions)) + layout_karyogram(gr_karyogram_regions, aes(color=pred_class, fill=pred_class)) + 
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
      ggsave(paste("../../results/classifier_plots_karyograms/", pop, "-karyogram_sweep_regions-", label, plout, sep=""))
      # autoplot(seqinfo(gr_karyogram_regions)) + 
      #   layout_karyogram(gr_karyogram_regions, aes(color=-softsweep_hardsweep_value, fill=-softsweep_hardsweep_value)) + 
      #   scale_fill_distiller(palette = "RdBu") + scale_color_distiller(palette = "RdBu")
    }
    write_tsv(tb_karyogram_regions, paste("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/", 
      pop, "-karyogram_sweep_regions-", label, ".tsv", sep=""))
    table(tb_karyogram_regions$pred_class)
    
    # plot manhattan for each chromosome
    for (chrm in chr_list) {
      temp_range <- gr_karyogram_regions[seqnames(gr_karyogram_regions)==paste("chr", chrm, sep="")]
      p <- ggplot() + 
        geom_line(data=filter(temp, chr==chrm), mapping=aes(x=physpos, y=-log10(neutral)), alpha=1.0, color="grey") + 
        geom_point(data=filter(temp, chr==chrm), mapping=aes(x=physpos, y=-log10(neutral), color=pred_class_alt), alpha=1.0) +
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
      if (length(temp_range) > 0) {
        p <- p + geom_rect(data=as_tibble(temp_range), 
          mapping=aes(xmin=start(temp_range), xmax=end(temp_range), ymin=-2, ymax=-1, fill=temp_range$pred_class)) + 
          scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
      }
      p; ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_chrm/", 
        pop, "-chr", chrm, "-posterior-", label, plout, sep=""))
    }
    
    # plot manhattan for known regions
    get_known_regions <- function(temp_data, temp_sweep, pop_name, gene_name, chr_name, start, end, buffer) {
      temp_data <- filter(temp, chr==chr_name & (physpos>=start-buffer) & (physpos<=end+buffer))
      temp_range <- (GRanges(seqnames=paste("chr", chr_name, sep=""), ranges=IRanges(start=start-buffer, end=end+buffer)))
      g <- temp_sweep[queryHits(findOverlaps(temp_sweep, temp_range))]
      p <- ggplot() + 
        geom_line(data=temp_data, mapping=aes(x=physpos, y=-log10(neutral)), alpha=1.0, color="grey", size=2) + 
        geom_vline(xintercept=start, size=1, alpha=0.8, linetype="dashed", color="grey") + 
        geom_vline(xintercept=end, size=1, alpha=0.8, linetype="dashed", color="grey") +
        geom_point(data=temp_data, mapping=aes(x=physpos, y=-log10(neutral), color=pred_class_alt), alpha=1.0, size=7) +
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        theme(panel.grid.major = element_blank()) + ggtitle(paste(gene_name, " (", pop_name, ")", sep="")) + 
        xlim(c(start-buffer, end+buffer)) + theme(axis.text.x = element_text(angle=45, hjust=0.4, vjust=0.4)) + 
        # ylab("-log10(1-sweep_neutral.posterior)") + xlab("Position") + theme(legend.position="none") + 
        ylab("-log10(neutral posterior)") + xlab(paste("Position on chr", chr_name, sep="")) + theme(legend.position="none") + 

        # 
        # # 
        # 

        theme(text=element_text(size=24)) + ylim(c(-2,-log10(.Machine$double.eps)))
      if (length(g) > 0) {
        p <- p + geom_rect(data=as_tibble(g), mapping=aes(xmin=start(g), xmax=end(g), ymin=-2, ymax=-1, fill=g$pred_class)) + 
          scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
          xlim(c(start-buffer, end+buffer))
      }
      return(p)
    }
    
    # if (TRUE) {
    if (pop == "CEU") {
      get_known_regions(temp, gr_karyogram_regions, pop, "ASPM", 1, 197053257, 197115824, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ASPM", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "LCT-MCM6", 2, 136545415, 136634047, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-LCT.MCM6", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TLR1", 4, 38797876, 38806814, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TLR1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC24A5" , 15, 48413169, 48434926, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC24A5", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "KITLG", 12, 88886570, 88974250, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-KITLG", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TRPV6", 7, 142568956, 142583490, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TRPV6", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "OCA2-HERC2", 15, 28000021, 28567313, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-OCA2.HERC2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC45A2", 5, 33944721, 33984835, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC45A2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "BNC2", 9, 16409501, 16870841, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-BNC2", "-posterior-", label, plout, sep=""))
      
      # voight et al. 2006
      get_known_regions(temp, gr_karyogram_regions, pop, "SPAG4", 20, 34203814, 34208971, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SPAG4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "ODF2", 9, 131217465, 131263571, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ODF2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "ACVR1", 2, 158592958, 158732374, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ACVR1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TGM4", 3, 44916100, 44956482, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TGM4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "MYO5A", 15, 52599480, 52821247, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-MYO5A", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "DTNBP1", 6, 15523032, 15663289, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-DTNBP1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TYRP1", 9, 12685439, 12710290, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TYRP1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "BMP3", 4, 81952119, 81978685, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-BMP3", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "BMPR2", 2, 203241659, 203432474, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-BMPR2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "GDF5", 20, 34021145, 34042568, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-GDF5", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC27A4", 9, 131102925, 131123749, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC27A4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "PPARD", 6, 35310335, 35395968, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-PPARD", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "CENPJ", 13, 25457171, 25497018, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-CENPJ", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC6A4", 17, 28521337, 28563020, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC6A4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SNTG1", 8, 50822349, 51706678, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SNTG1", "-posterior-", label, plout, sep=""))
    }
    
    # if (TRUE) {
    if (pop == "CHB") {
      get_known_regions(temp, gr_karyogram_regions, pop, "ABCC11", 16, 48200782, 48281479, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ABCC11", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "ADH1B", 4, 100227527, 100242572, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ADH1B", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "EDAR", 2, 109510927, 109605828, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-EDAR", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "MC1R", 16, 89984287, 89987385, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-MC1R", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "PCDH15", 10, 55562531, 56561051, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-PCDH15", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "KITLG", 12, 88886570, 88974250, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-KITLG", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TRPV6", 7, 142568956, 142583490, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TRPV6", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "OCA2-HERC2", 15, 28000021, 28567313, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-OCA2.HERC2", "-posterior-", label, plout, sep=""))
      
      # voight et al. 2006
      get_known_regions(temp, gr_karyogram_regions, pop, "RSBN1", 1, 114304454, 114355098, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-RSBN1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SPAG4", 20, 34203814, 34208971, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SPAG4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "BMP5", 6, 55618443, 55740362, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-BMP5", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "GDF5", 20, 34021145, 34042568, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-GDF5", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "MAN2A1", 5, 109025067, 109205326, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-MAN2A1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SI", 3, 164696686, 164796283, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SI", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC25A20", 3, 48894369, 48936426, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC25A20", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "LEPR", 1, 65886248, 66107242, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-LEPR", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "CENPJ", 13, 25457171, 25497018, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-CENPJ", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC6A4", 17, 28521337, 28563020, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC6A4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SNTG1", 8, 50822349, 51706678, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SNTG1", "-posterior-", label, plout, sep=""))
    }
    
    # if (TRUE) {
    if (pop == "YRI") {
      get_known_regions(temp, gr_karyogram_regions, pop, "ACKR1 (DARC)", 1, 159173803, 159176290, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-ACKR1.DARC", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "HBB", 11, 5246696, 5248301, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-HBB", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "PSCA", 8, 143751726, 143764145, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-PSCA", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "RBFOX2", 22, 36134783, 36424679, 2.5e6)  # sugden 2018
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-RBFOX2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SHPK-CTNS-EMC6-TRPV1-ITGAE-GSG2-TAX1BP3", 17, 3468740, 3704537, 2.5e6)  # Sugden 2018
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SHPK.CTNS.EMC6.TRPV1.ITGAE.GSG2.TAX1BP3", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "HS3ST2-USP31", 16, 22825860, 23160591, 2.5e6)  # added, see Pickrell 2009, Granka 2012, Sugden 2018
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-HS3ST2.USP31", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "KITLG ", 12, 88886570, 88974250, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-KITLG", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "TRPV6", 7, 142568956, 142583490, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-TRPV6", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "OCA2-HERC2", 15, 28000021, 28567313, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-OCA2.HERC2", "-posterior-", label, plout, sep=""))
      
      # ferrer-admetlla et al. 2014
      get_known_regions(temp, gr_karyogram_regions, pop, "IL34", 16, 70613798, 70694585, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-IL34", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "APOL1", 22, 36649056, 36663576, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-APOL1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "CD36", 7, 79998891, 80308593, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-CD36", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "OSBP2", 22, 31089769, 31303811, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-OSBP2", "-posterior-", label, plout, sep=""))
      
      # harris et al. 2018
      get_known_regions(temp, gr_karyogram_regions, pop, "KIAA0825", 5, 93488671, 93954309, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-KIAA0825", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "NNT", 5, 43602794, 43707507, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-NNT", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "HEMGN", 9, 100689073, 100707138, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-HEMGN", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "RGS18", 1, 192127587, 192154945, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-RGS18", "-posterior-", label, plout, sep=""))
      
      # voight et al. 2006
      get_known_regions(temp, gr_karyogram_regions, pop, "RSBN1", 1, 114304454, 114355098, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-RSBN1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "CPEB2", 4, 15004298, 15071777, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-CPEB2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "FZD6", 8, 104310661, 104345094, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-FZD6", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "MAN2A1", 5, 109025067, 109205326, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-MAN2A1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "NCOA1", 2, 24714783, 24993571, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-NCOA1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "CDK5RAP2", 9, 123151147, 123342448, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-CDK5RAP2", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "GABRA4", 4, 46920917, 46996424, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-GABRA4", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "PSEN1", 14, 73603126, 73690399, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-PSEN1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SYT1", 12, 79257773, 79845788, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SYT1", "-posterior-", label, plout, sep=""))
      get_known_regions(temp, gr_karyogram_regions, pop, "SNTG1", 8, 50822349, 51706678, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SNTG1", "-posterior-", label, plout, sep=""))
      
      # szpak et al. 2018
      get_known_regions(temp, gr_karyogram_regions, pop, "SLC39A4", 8, 145635126, 145642279, 2.5e6)
      ggsave(paste("../../results/classifier_plots_nb_region/", pop, "-SLC39A4", "-posterior-", label, plout, sep=""))
    }
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
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.98-simulated", T)
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.98-empirical", F)

prob_neutral <- 0.90
prob_prior <-
  c("neutral"=prob_neutral,
    "hardsweep"=(1-prob_neutral)/2,
    "softsweep"=(1-prob_neutral)/2)
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.90-simulated", T)
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.90-empirical", F)

prob_neutral <- 0.34
prob_prior <-
  c("neutral"=prob_neutral,
    "hardsweep"=(1-prob_neutral)/2,
    "softsweep"=(1-prob_neutral)/2)
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.34-simulated", T)
application_helper(prob_prior, stat_temp, stat_sweep_final, "0.34-empirical", F)

# # # sweep windows
sweeps_list_windows_empirical <- NULL
sweeps_list_windows_empirical$low_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.34-empirical.tsv")
sweeps_list_windows_empirical$low_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.34-empirical.tsv")
sweeps_list_windows_empirical$low_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.34-empirical.tsv")
sweeps_list_windows_empirical$lower_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.90-empirical.tsv")
sweeps_list_windows_empirical$lower_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.90-empirical.tsv")
sweeps_list_windows_empirical$lower_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.90-empirical.tsv")
sweeps_list_windows_empirical$lowest_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.98-empirical.tsv")
sweeps_list_windows_empirical$lowest_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.98-empirical.tsv")
sweeps_list_windows_empirical$lowest_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.98-empirical.tsv")

prior_map <- c("low_windows"=0.34, "lower_windows"=0.9, "lowest_windows"=0.98)
for (i in names(prior_map)) {
  for (j in pop_list) {
    sweeps_list_windows_empirical[[i]][[j]] <- sweeps_list_windows_empirical[[i]][[j]] %>% mutate(prior=prior_map[[i]], population=j)
  }
}

sweep_table_windows_empirical <- NULL
for (i in names(prior_map)) {
  for (j in pop_list) {
    if (is.null(sweep_table_windows_empirical)) {
      sweep_table_windows_empirical <- sweeps_list_windows_empirical[[i]][[j]]
    }
    sweep_table_windows_empirical <- rbind(sweep_table_windows_empirical, sweeps_list_windows_empirical[[i]][[j]])
  }
}

sweep_table_windows_empirical$population <- factor(sweep_table_windows_empirical$population, pop_list)
sweep_table_windows_empirical$pred_class_alt <- 
  factor(sweep_table_windows_empirical$pred_class_alt, c("softsweep", "softsweepweak", "hardsweepweak", "hardsweep"))

sweeps_list_windows_simulated <- NULL
sweeps_list_windows_simulated$low_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.34-simulated.tsv")
sweeps_list_windows_simulated$low_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.34-simulated.tsv")
sweeps_list_windows_simulated$low_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.34-simulated.tsv")
sweeps_list_windows_simulated$lower_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.90-simulated.tsv")
sweeps_list_windows_simulated$lower_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.90-simulated.tsv")
sweeps_list_windows_simulated$lower_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.90-simulated.tsv")
sweeps_list_windows_simulated$lowest_windows$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/YRI-karyogram_sweep_windows-0.98-simulated.tsv")
sweeps_list_windows_simulated$lowest_windows$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CEU-karyogram_sweep_windows-0.98-simulated.tsv")
sweeps_list_windows_simulated$lowest_windows$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_windows/CHB-karyogram_sweep_windows-0.98-simulated.tsv")

prior_map <- c("low_windows"=0.34, "lower_windows"=0.9, "lowest_windows"=0.98)
for (i in names(prior_map)) {
  for (j in pop_list) {
    sweeps_list_windows_simulated[[i]][[j]] <- sweeps_list_windows_simulated[[i]][[j]] %>% mutate(prior=prior_map[[i]], population=j)
  }
}

sweep_table_windows_simulated <- NULL
for (i in names(prior_map)) {
  for (j in pop_list) {
    if (is.null(sweep_table_windows_simulated)) {
      sweep_table_windows_simulated <- sweeps_list_windows_simulated[[i]][[j]]
    } else {
      sweep_table_windows_simulated <- rbind(sweep_table_windows_simulated, sweeps_list_windows_simulated[[i]][[j]])
    }
  }
}

sweep_table_windows_simulated$population <- factor(sweep_table_windows_simulated$population, pop_list)
sweep_table_windows_simulated$pred_class_alt <- 
  factor(sweep_table_windows_simulated$pred_class_alt, c("softsweep", "softsweepweak", "hardsweepweak", "hardsweep"))

sweep_table_windows <- rbind(
  sweep_table_windows_simulated %>% mutate(label="Simulated"),
  sweep_table_windows_empirical %>% mutate(label="Empirical"))
sweep_table_windows$label <- 
  factor(sweep_table_windows$label, c("Simulated", "Empirical"))
sweep_table_windows_final <- filter(sweep_table_windows, prior=="0.98" & label=="Empirical")

ggplot(data=sweep_table_windows) + geom_bar(mapping=aes(x=as.factor(prior), fill=pred_class_alt, alpha=label), position="stack") + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  scale_alpha_manual(name="Classifier", breaks=c("Simulated", "Empirical"), labels=c("Simulated", "Empirical"), values=c(1, 1), guide="none") + 
  xlab("Prior probability of neutral class") + ylab("Count") + facet_grid(~population + label) + 
  theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) + theme(legend.position = c(0.87, 0.8)) + theme(aspect.ratio=4)
ggsave(paste("../../results/classifier_plots_miscellaneous/sweep_windows_comparison_classified_prior_pop", plout, sep=""))
ggplot(data=sweep_table_windows) + geom_bar(mapping=aes(x=as.factor(prior), fill=pred_class_alt, alpha=label), position="stack") + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  scale_alpha_manual(name="Classifier", breaks=c("Simulated", "Empirical"), labels=c("Simulated", "Empirical"), values=c(1, 1), guide="none") + 
  xlab("Prior probability of neutral class") + ylab("Count") + facet_grid(~population + label) + 
  theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5)) 
ggsave(paste("../../results/classifier_plots_miscellaneous/sweep_windows_comparison_classified_prior_pop_alt", plout, sep=""))

# # # sweep regions
sweeps_list_regions_empirical <- NULL
sweeps_list_regions_empirical$low_regions$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/YRI-karyogram_sweep_regions-0.34-empirical.tsv")
sweeps_list_regions_empirical$low_regions$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CEU-karyogram_sweep_regions-0.34-empirical.tsv")
sweeps_list_regions_empirical$low_regions$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CHB-karyogram_sweep_regions-0.34-empirical.tsv")
sweeps_list_regions_empirical$lower_regions$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/YRI-karyogram_sweep_regions-0.90-empirical.tsv")
sweeps_list_regions_empirical$lower_regions$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CEU-karyogram_sweep_regions-0.90-empirical.tsv")
sweeps_list_regions_empirical$lower_regions$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CHB-karyogram_sweep_regions-0.90-empirical.tsv")
sweeps_list_regions_empirical$lowest_regions$YRI <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/YRI-karyogram_sweep_regions-0.98-empirical.tsv")
sweeps_list_regions_empirical$lowest_regions$CEU <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CEU-karyogram_sweep_regions-0.98-empirical.tsv")
sweeps_list_regions_empirical$lowest_regions$CHB <- 
  read_tsv("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_regions/CHB-karyogram_sweep_regions-0.98-empirical.tsv")

prior_map <- c("low_regions"=0.34, "lower_regions"=0.9, "lowest_regions"=0.98)
for (i in names(prior_map)) {
  for (j in pop_list) {
    sweeps_list_regions_empirical[[i]][[j]] <- sweeps_list_regions_empirical[[i]][[j]] %>% mutate(prior=prior_map[[i]], population=j)
  }
}

sweep_table_regions_empirical <- NULL
for (i in names(prior_map)) {
  for (j in pop_list) {
    if (is.null(sweep_table_regions_empirical)) {
      sweep_table_regions_empirical <- sweeps_list_regions_empirical[[i]][[j]]
    } else {
      sweep_table_regions_empirical <- rbind(sweep_table_regions_empirical, sweeps_list_regions_empirical[[i]][[j]])
    }
  }
}

sweep_table_regions_empirical$population <- factor(sweep_table_regions_empirical$population, pop_list)
sweep_table_regions_empirical$pred_class <- 
  factor(sweep_table_regions_empirical$pred_class, c("softsweep", "softsweepweak", "hardsweepweak", "hardsweep"))
sweep_table_regions_final <- filter(sweep_table_regions_empirical, prior==0.98)
sweep_regions_counts <- table(sweep_table_regions_final[c("population", "pred_class")])
write_tsv(as_tibble(sweep_regions_counts), 
          "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/sweep_regions_counts_final_pop_class.tsv")

ggplot(data=filter(sweep_table_regions_empirical, prior %in% c(0.90, 0.98))) + 
  geom_bar(mapping=aes(x=population, fill=pred_class)) + facet_grid(~as.factor(prior)) + theme(text=element_text(size=24)) + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code, guide="none") + 
  xlab("Population") + ylab("Count")
ggsave(paste("../../results/classifier_plots_miscellaneous/sweep_regions_final_classified_pop", plout, sep=""))
ggplot(data=filter(sweep_table_regions_empirical, prior %in% c(0.90, 0.98))) + 
  geom_bar(mapping=aes(x=population, fill=pred_class)) + facet_grid(~as.factor(prior)) + theme(text=element_text(size=24)) + 
  scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
  xlab("Population") + ylab("Count")
ggsave(paste("../../results/classifier_plots_miscellaneous/sweep_regions_final_classified_pop_alt", plout, sep=""))

tb_ucsc_genes <- read_tsv("../../data/applications/1000GP_Phase1_misc/hg19-ucsc_genes-known_canonical.txt")
names(tb_ucsc_genes) <- c("chrom", "chromStart", "chromEnd", "clusterId", "protein", "geneSymbol", "description")
gr_ucsc_genes <- GRanges(seqnames=tb_ucsc_genes$chrom,
  IRanges(start=tb_ucsc_genes$chromStart, end=tb_ucsc_genes$chromEnd), gene=tb_ucsc_genes$description)
tb_ucsc_genes <- as_tibble(gr_ucsc_genes)
overlaps <- as_tibble(findOverlaps(GRanges(sweep_table_regions_final), GRanges(tb_ucsc_genes)))
overlaps <- overlaps %>% mutate(gene=tb_ucsc_genes$gene[subjectHits])
nearest <- as_tibble(nearest(GRanges(sweep_table_regions_final), GRanges(tb_ucsc_genes)))
nearest <- nearest %>% mutate(gene=tb_ucsc_genes$gene[value])
sweep_table_regions_final$overlap_genes <- NA
for (i in 1:nrow(sweep_table_regions_final)) {
  query <- filter(overlaps, queryHits==i)
  if (nrow(query)==0) {
    sweep_table_regions_final$overlap_genes[i] <- paste("~", as.character(nearest$gene[i]), sep="")
  } else {
    sweep_table_regions_final$overlap_genes[i] <- paste(as.character(filter(overlaps, queryHits==i)$gene), collapse=", ")
  }
}
sweep_table_regions_final <- dplyr::select(sweep_table_regions_final, 
  c(population, seqnames, start, end, width, pred_class, sweep_neutral_value, softsweep_hardsweep_value, overlap_genes))
write_tsv(sweep_table_regions_final, "../../results/classifier_plots_miscellaneous/sweep_regions_final_genes.tsv")

for (p in pop_list) {
  for (c in chr_list) {
    if (nrow(filter(sweep_table_regions_final, seqnames==paste("chr", c, sep="") & population==p))>0) {
      ggplot(filter(sweep_table_regions_final, seqnames==paste("chr", c, sep="") & population==p)) + 
        geom_segment(mapping=aes(x=start, xend=end, color=pred_class), y=0, yend=0, size=2) + 
        geom_point(mapping=aes(x=(start+end)/2, y=sweep_neutral_value, color=pred_class)) + 
        geom_segment(mapping=aes(x=(start+end)/2, xend=(start+end)/2, y=sweep_neutral_value, color=pred_class), yend=0) + 
        geom_text(mapping=aes(x=(start+end)/2, y=0, label=overlap_genes), angle=90, size=3, hjust=0, vjust=-0.2) + 
        scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
        xlim(c(0, seqlengths(seqinfo(BSgenome.Hsapiens.UCSC.hg19))[paste("chr", c, sep="")]))
      ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_label/", p, "-chr", c, "-sweep_region_final_label", plout, sep=""))
    }
  }
}

for (c in chr_list) {
  if (nrow(filter(sweep_table_regions_final, seqnames==paste("chr", c, sep="")))>0) {
    ggplot(filter(sweep_table_regions_final, seqnames==paste("chr", c, sep=""))) + 
      geom_segment(mapping=aes(x=start, xend=end, color=pred_class), y=0, yend=0, size=2) + 
      geom_point(mapping=aes(x=(start+end)/2, y=sweep_neutral_value, color=pred_class)) + 
      geom_segment(mapping=aes(x=(start+end)/2, xend=(start+end)/2, y=sweep_neutral_value, color=pred_class), yend=0) + 
      geom_text(mapping=aes(x=(start+end)/2, y=0, label=overlap_genes), angle=90, size=1, hjust=0, vjust=-0.2) + 
      scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code) + 
      xlim(c(0, seqlengths(seqinfo(BSgenome.Hsapiens.UCSC.hg19))[paste("chr", c, sep="")])) + 
      facet_wrap(~population, nrow=3) + theme(aspect.ratio=1) + ggtitle(paste("Chrm ", c, sep=""))
    ggsave(paste("../../results/classifier_plots_miscellaneous/classifier_plots_sweep_label/chr", c, "-sweep_region_final_label", plout, sep="")) 
  }
}

# # # proportions plot
proportions_list <- NULL

prediction_list_windows <- NULL
prediction_list_windows$YRI <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/YRI-prediction_table_windows-0.98-empirical.tsv", sep=""))
prediction_list_windows$CEU <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CEU-prediction_table_windows-0.98-empirical.tsv", sep=""))
prediction_list_windows$CHB <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CHB-prediction_table_windows-0.98-empirical.tsv", sep=""))
prediction_table_windows_final <- rbind(prediction_list_windows$YRI, prediction_list_windows$CEU, prediction_list_windows$CHB)
prediction_list_windows$pop <- factor(prediction_list_windows$pop, pop_list)
prediction_list_windows$pred_class_alt <- 
  factor(prediction_list_windows$pred_class_alt, c("softsweep", "softsweepweak", "neutral", "hardsweepweak", "hardsweep"))
prediction_windows_counts <- table(prediction_table_windows_final[c("pop", "pred_class_alt")])
write_tsv(as_tibble(prediction_windows_counts), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_counts_0.98-empirical_pop_class.tsv")
prediction_windows_props <- prediction_windows_counts/rowSums(prediction_windows_counts)
write_tsv(as_tibble(prediction_windows_props), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_proportions_0.98-empirical_pop_class.tsv")
proportions_list$lowest <- prediction_windows_props %>% as_tibble() %>% mutate(prior=0.98)

prediction_list_windows <- NULL
prediction_list_windows$YRI <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/YRI-prediction_table_windows-0.90-empirical.tsv", sep=""))
prediction_list_windows$CEU <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CEU-prediction_table_windows-0.90-empirical.tsv", sep=""))
prediction_list_windows$CHB <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CHB-prediction_table_windows-0.90-empirical.tsv", sep=""))
prediction_table_windows_final <- rbind(prediction_list_windows$YRI, prediction_list_windows$CEU, prediction_list_windows$CHB)
prediction_list_windows$pop <- factor(prediction_list_windows$pop, pop_list)
prediction_list_windows$pred_class_alt <- 
  factor(prediction_list_windows$pred_class_alt, c("softsweep", "softsweepweak", "neutral", "hardsweepweak", "hardsweep"))
prediction_windows_counts <- table(prediction_table_windows_final[c("pop", "pred_class_alt")])
write_tsv(as_tibble(prediction_windows_counts), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_counts_0.90-empirical_pop_class.tsv")
prediction_windows_props <- prediction_windows_counts/rowSums(prediction_windows_counts)
write_tsv(as_tibble(prediction_windows_props), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_proportions_0.90-empirical_pop_class.tsv")
proportions_list$lower <- prediction_windows_props %>% as_tibble() %>% mutate(prior=0.90)

prediction_list_windows <- NULL
prediction_list_windows$YRI <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/YRI-prediction_table_windows-0.34-empirical.tsv", sep=""))
prediction_list_windows$CEU <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CEU-prediction_table_windows-0.34-empirical.tsv", sep=""))
prediction_list_windows$CHB <- 
  read_tsv(paste("../../results/classifier_plots_miscellaneous/classifier_plots_nb_pred/CHB-prediction_table_windows-0.34-empirical.tsv", sep=""))
prediction_table_windows_final <- rbind(prediction_list_windows$YRI, prediction_list_windows$CEU, prediction_list_windows$CHB)
prediction_list_windows$pop <- factor(prediction_list_windows$pop, pop_list)
prediction_list_windows$pred_class_alt <- 
  factor(prediction_list_windows$pred_class_alt, c("softsweep", "softsweepweak", "neutral", "hardsweepweak", "hardsweep"))
prediction_windows_counts <- table(prediction_table_windows_final[c("pop", "pred_class_alt")])
write_tsv(as_tibble(prediction_windows_counts), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_counts_0.34-empirical_pop_class.tsv")
prediction_windows_props <- prediction_windows_counts/rowSums(prediction_windows_counts)
write_tsv(as_tibble(prediction_windows_props), "../../results/classifier_plots_miscellaneous/classifier_plots_prediction_windows/prediction_windows_proportions_0.34-empirical_pop_class.tsv")
proportions_list$low <- prediction_windows_props %>% as_tibble() %>% mutate(prior=0.34)

proportions_list_final <- rbind(proportions_list$low, proportions_list$lower, proportions_list$lowest)
ggplot(data=filter(proportions_list_final, pred_class_alt=="neutral")) + 
  geom_point(mapping=aes(x=as.numeric(prior), y=n, group=pop, color=pop)) + 
  geom_line(mapping=aes(x=as.numeric(prior), y=n, group=pop, color=pop)) + 
  geom_abline(intercept=0, slope=1, color="darkgrey", linetype="dashed") + 
  scale_x_continuous(name="Prior probabilty of neutral class", breaks=c(0.34, 0.90, 0.98), labels=c(0.34, 0.90, 0.98)) + 
  ylab("Proportion predicted in neutral class") + guides(color=guide_legend(title="Population")) + 
  theme(panel.grid.major = element_blank()) + theme(text=element_text(size=24)) + theme(legend.position = c(0.2, 0.8))
ggsave(paste("../../results/classifier_plots_miscellaneous/proportion_windows_prior_predicted", plout, sep=""))
ggplot(data=filter(proportions_list_final, pred_class_alt=="neutral")) + 
  geom_point(mapping=aes(x=as.numeric(prior), y=n, group=pop, color=pop)) + 
  geom_line(mapping=aes(x=as.numeric(prior), y=n, group=pop, color=pop)) + 
  geom_abline(intercept=0, slope=1, color="darkgrey", linetype="dashed") + 
  scale_x_continuous(name="Prior probabilty of neutral class", breaks=c(0.34, 0.90, 0.98), labels=c(0.34, 0.90, 0.98)) + 
  ylab("Proportion predicted in neutral class") + guides(color=guide_legend(title="Population")) + 
  theme(panel.grid.major = element_blank()) + theme(text=element_text(size=24))
ggsave(paste("../../results/classifier_plots_miscellaneous/proportion_windows_prior_predicted_alt", plout, sep=""))

# # # 

# # summary statistics plots
# ggplot(temp) + geom_boxplot(mapping=aes(x=pred_class_alt, y=H2H1_win, fill=pred_class_alt), alpha=1) + 
#   scale_fill_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
# ggplot(temp) + geom_point(mapping=aes(x=H12_win, y=H2H1_win, color=pred_class_alt)) + 
#   scale_color_manual(name=stat_code["site_class"], breaks=names(color_code), labels=class_code, values=color_code)
