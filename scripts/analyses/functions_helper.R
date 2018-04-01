library(RColorBrewer)  # brewer.pal function
library(tidyverse)  # data frame tools
library(wrapr)  # variable names in ggplot
# library(ggalt)  # more ggplot options
# library(ggfortify)  # more ggplot options
library(ggthemr)  # more ggplot themes
# library(ggthemes)  # more ggplot themes
# library(ggpubr)  # pub quality ggplot

# # # 

options(warn=1)
pop_list <- c("YRI", "CEU", "CHB")
plout <- ".pdf"  # .tiff, .jpg, .png, .pdf

# # # custom theme

# ggthemr("fresh")
fresh_edit <- define_palette(
  swatch = structure(c('#111111', '#65ADC2', '#233B43', '#E84646', 
    '#C29365','#362C21', '#316675','#168E7F', '#109B37')),
  gradient = list(low='#65ADC2', high='#362C21')
)
fresh_edit$text <- list(inner = '#000000', outer = '#000000')
fresh_edit$line <- list(inner = '#000000', outer = '#000000')
# fresh_edit$gridline = '#C3C3C3'
fresh_edit$gridline = '#EEE4DA'
ggthemr(fresh_edit)

# # # loading functions

label_app_stan <- "1000GP_Phase1_stan"
label_app_main <- "1000GP_Phase1_main"
label_app_file <- "1000GPP1"

label_sim_stan  = "2017111501"
label_sim_train = "2017111502"  # first half
label_sim_test  = "2017111502"  # second half
label_sim_link  = "2017111503"

# list of simulated files
list_sim_files <- function(label, start, end, pop) {
  my_files_sim <- NULL
  my_files_sim[["metadata"]] <- 
    mapply(function(i) paste("../../data/simulations/", label, "/", label, "_", i, "_", pop, "_metadata.txt", sep=""), start:end)
  my_files_sim[["features"]] <- 
    mapply(function(i) paste("../../data/simulations/", label, "/", label, "_", i, "_", pop, "_features.txt", sep=""), start:end)
  return(my_files_sim) 
}

# list of application files
list_app_files <- function(label, label_x, chr_list, pop) {
  my_files_app <- NULL
  my_files_app[["metadata"]] <- mapply(function(chr) 
    paste("../../data/applications/", label, "/", label_x, "_chr", chr, "_pop", pop, "_metadata.txt", sep=""), chr_list)
  my_files_app[["features"]] <- mapply(function(chr) 
    paste("../../data/applications/", label, "/", label_x, "_chr", chr, "_pop", pop, "_features.txt", sep=""), chr_list)
  return(my_files_app) 
}

# save train and test files
my_files_sim_stan <- NULL
my_files_sim_train <- NULL
my_files_sim_test <- NULL
my_files_sim_link <- NULL
for (pop in pop_list) {
  my_files_sim_stan[[pop]]  <- list_sim_files(label_sim_stan, 0, 0, pop)  # 0, 3
  my_files_sim_train[[pop]] <- list_sim_files(label_sim_train, 0, 199, pop)  # 0, 199
  my_files_sim_test[[pop]]  <- list_sim_files(label_sim_test, 200, 399, pop)  # 200, 399
  my_files_sim_link[[pop]]  <- list_sim_files(label_sim_link, 0, 39, pop)  # 0, 39
}

# save application files
my_files_app_stan <- NULL
my_files_app_main <- NULL
chr_list = c(1:22)  # 1, 22
for (pop in pop_list) {
  my_files_app_stan[[pop]] <- list_app_files(label_app_stan, label_app_file, chr_list, pop)
  my_files_app_main[[pop]] <- list_app_files(label_app_main, label_app_file, chr_list, pop)
}

# load data from tsv files
load_data <- function(my_files, pop) {
  my_table <- NULL
  for (m in my_files[[pop]][["metadata"]]) {
    if (!exists("my_metadata")) {my_metadata <- read_tsv(m)}
    else {my_metadata <- read_tsv(m) %>% full_join(my_metadata, .)}
  }
  for (f in my_files[[pop]][["features"]]) {
    if (!exists("my_features")) {my_features <- read_tsv(f)}
    else {my_features <- read_tsv(f) %>% full_join(my_features, .)}
  }
  my_table <- inner_join(my_metadata, my_features, by="key")
  return(my_table)
}

# # # preprocessing functions

# list of statistics

stat_misc <- c("segsites_win", "theta_Pi_win", "K_win", "recrates_win")

stat_list <- c("FST_win", "Taj_D_win", "FayWu_H_win", "Zeng_E_win", "H2_win", "H12_win", 
               "H2H1_win", "ZA_win", "Wall_B_win", "iHS_win", "XPEHH_win", "D_theta_Pi_win")

# preprocessing steps
preprocess_data <- function(my_table) {
  # convert cols to factors
  my_table <- my_table %>% 
    mutate_at(vars(key, chr, pop, stype, model, spop, complete, sweep_freq_bin, site_class), funs(as.factor)) 
  # round dist_physpos values
  my_table <- my_table %>% 
    mutate(dist_physpos=ifelse(!is.na(dist_physpos), plyr::round_any(as.numeric(dist_physpos), 100000), dist_physpos)) %>% 
  # filter rows with low recrates and low segsites
    filter(segsites_win >= 2) %>%
    filter(recrates_win >= 0.2) %>%
  # # filter rows with high recrates and high segsites
    filter(recrates_win <= 20) %>%
    filter(segsites_win <= 1000) %>%
  # # removed, naive bayes handles missing values
  # # filter rows with missing values in stat_list
  #   filter_at(vars(stat_list), all_vars(!is.na(.))) %>% 
  # log transform iHS_win and XPEHH_win
    # mutate(iHS_win = log10(iHS_win)) %>% 
    # mutate(XPEHH_win = log10(XPEHH_win))
  # replace zeros with -Inf, but do not log transform
    mutate(iHS_win = ifelse(my_tol(iHS_win, 0), -Inf, iHS_win)) %>% 
    mutate(XPEHH_win = ifelse(my_tol(XPEHH_win, 0), -Inf, XPEHH_win)) %>% 
    mutate(iHS_win = ifelse(my_tol(iHS_win, 1), Inf, iHS_win)) %>% 
    mutate(XPEHH_win = ifelse(my_tol(XPEHH_win, 1), Inf, XPEHH_win)) %>% 
  # bin selection parameters with three bins
    mutate(stime_bin=as.factor(ifelse(!is.na(stime), ifelse(stime>=15000, ifelse(stime>=25000, "3", "2"), "1"), "0"))) %>% 
    mutate(scoeff_bin=as.factor(ifelse(!is.na(scoeff), ifelse(scoeff>=0.005, ifelse(scoeff>=0.05, "3", "2"), "1"), "0"))) %>% 
    mutate(sfreq_bin=as.factor(ifelse(!is.na(sfreq), ifelse(sfreq>=0.0005, ifelse(sfreq>=0.005, "3", "2"), "1"), "0"))) %>% 
    mutate(complete_bin=as.factor(ifelse(!is.na(complete), ifelse(complete==1, "2", "1"), "0"))) %>% 
  # bin pre and post sweep parameters
    mutate(sprefix_freq=as.factor(ifelse(complete==0, ifelse(sweep_freq_pop>=0.33, ifelse(sweep_freq_pop>=0.67, "3", "2"), "1"), NA))) %>% 
    mutate(spostfix_time=as.factor(ifelse(complete==1, ifelse(sweep_fix_time>=10000, ifelse(sweep_fix_time>=20000, "3", "2"), "1"), NA)))
  # remove regions around acen and gvar
  # return(preprocess_cyto(my_table))
  return(my_table)
}

# # # mapping objects

# stable names
stat_code <- c("key"="key", "chr"="chr", "pop"="pop", "stype"="stype", "model"="model", "seed"="seed",
               "mutrate"="Mutation rate", "recrate"="Recombination rate", "spop"="Sweep population", 
               "stime"="Start time of sweep", "scoeff"="Selection coefficient", 
               "sfreq"="Initial frequency", "complete"="Completeness of sweep", 
               "sweep_freq_bin"="Final sample frequency", "sweep_freq_pop"="Final population frequency",
               "sweep_fix_time"="Fixation time of sweep", "sweep_dur_time"="Duration time of sweep",
               "origins"="Number of distinct origins", "physpos"="Physical position", 
               "dist_physpos"="Distance to sweep (Kb)", "site_class"="Sweep class",
               "segsites_win"="Number of segregating sites", "recrates_win"="Average recombination rate", 
               "theta_Pi_win"="Average nucleotide diversity", "K_win"="Number of distinct haplotypes", 
               "FST_win"="FST", "Taj_D_win"="Tajima's D", "FayWu_H_win"="Fay and Wu's H", "Zeng_E_win"="Zeng's E",
               "H2_win"="H2", "H12_win"="H12", "H2H1_win"="H2/H1", "ZA_win"="ZA", "Wall_B_win"="Wall's B", 
               "iHS_win"="|iHS|>2", "XPEHH_win"="|XPEHH|>2.5", "D_theta_Pi_win"="Delta Pi")

class_code <- c("neutral"="Neutral", "missing"="Missing", 
                "sweepweak"="Sweep-ish", "sweep"="Sweep", 
                "softsweep"="Soft sweep (BF>=10)", "hardsweep"="Hard sweep (BF>=10)", 
                "softlinked"="Soft linked", "hardlinked"="Hard linked", 
                "softsweepweak"="Soft sweep (BF<10)", "hardsweepweak"="Hard sweep (BF<10)", 
                "nonneutral_alt"="Non-neutral", "noneutral"="Non-neutral")

# my_palette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
#                 "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B2B2B2", "#000000")
# barplot(rep(1, length(my_palette)), col=my_palette)

color_code <- c("neutral"="#33A02C", "missing"="#B2B2B2", 
                "sweepweak"="#CAB2D6", "sweep"="#6A3D9A", 
                "softsweep"="#1F78B4", "hardsweep"="#FF7F00", 
                "softlinked"="#A6CEE3", "hardlinked"="#FDBF6F", 
                "softsweepweak"="#A6CEE3", "hardsweepweak"="#FDBF6F", 
                "nonneutral_alt"="#FB9A99", "nonneutral"="#E31A1C")

stime_code <- c("0"="neutral", "1"="[5kya, 15kya)", "2"="[15kya, 25kya)", "3"="[25kya, 35kya)")
scoeff_code <- c("0"="neutral", "1"="[0.0005, 0.005)", "2"="[0.005, 0.05)", "3"="[0.05, 0.5)")
sfreq_code <- c("0"="neutral", "1"="[1/(2Ne), 0.0005)", "2"="[0.0005, 0.005)", "3"="[0.005, 0.05)")

# multiplot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout==i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

# extend genomic range
# https://support.bioconductor.org/p/78652/
extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x)=="*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x)=="+" | strand(x)=="*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

my_tol <- function(i, j) {
  return(abs(i-j)<=.Machine$double.eps)
}
