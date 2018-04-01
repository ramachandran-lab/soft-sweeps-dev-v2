source("functions_helper.R")
set.seed(2017070901)

# # # 

# # # standardization

my_table_sim_stan <- NULL
my_table_app_stan <- NULL
for (pop in pop_list){
  # app/sim load
  my_table_app_stan[[pop]] <- load_data(my_files_app_stan, pop)
  my_table_sim_stan[[pop]] <- load_data(my_files_sim_stan, pop)
  my_table_app_stan[[pop]] <- my_table_app_stan[[pop]] %>% 
    mutate(XPEHH_12_snp=as.numeric(XPEHH_12_snp)) %>% 
    mutate(XPEHH_23_snp=as.numeric(XPEHH_23_snp)) %>% 
    mutate(XPEHH_31_snp=as.numeric(XPEHH_31_snp))
  my_table_sim_stan[[pop]] <- my_table_sim_stan[[pop]] %>% 
    mutate(XPEHH_12_snp=as.numeric(XPEHH_12_snp)) %>% 
    mutate(XPEHH_23_snp=as.numeric(XPEHH_23_snp)) %>% 
    mutate(XPEHH_31_snp=as.numeric(XPEHH_31_snp))
}
my_table_app_stan <- bind_rows(my_table_app_stan[pop_list])
my_table_sim_stan <- bind_rows(my_table_sim_stan[pop_list])

cone = color_code["neutral"]
# SFS histogram
ggplot() +
  geom_density(data=my_table_sim_stan, mapping=aes(x=DAF_snp, y=(..count..)/sum(..count..)), color=cone, fill=cone, alpha=0.2) +
  geom_density(data=my_table_app_stan, mapping=aes(x=DAF_snp, y=(..count..)/sum(..count..)), color="black", linetype="dashed") + 
  guides(color=F) + xlab("DAF") + ylab("Proportion") + facet_wrap(~pop) + theme(aspect.ratio=1) # + ggtitle(pop) +
ggsave(file=paste("../../results/summary_plots_stan/", "SFS_histogram", plout, sep=""))

# iHS density
ggplot() +
  geom_density(data=my_table_sim_stan, mapping=aes(x=iHS_snp), color=cone, fill=cone, alpha=0.2) +
  geom_density(data=my_table_app_stan, mapping=aes(x=iHS_snp), color="black", linetype="dashed") +
  guides(color=F) + xlab("iHS") + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) # + ggtitle(pop) +
ggsave(file=paste("../../results/summary_plots_stan/", "iHS_density", plout, sep=""))

# XPEHH density
ggplot() +
  geom_density(data=my_table_sim_stan, mapping=aes(x=XPEHH_12_snp), color=cone, fill=cone, alpha=0.2) +
  geom_density(data=my_table_app_stan, mapping=aes(x=XPEHH_12_snp), color="black", linetype="dashed") +
  guides(color=F) + xlab("XPEHH (AF-EU)") + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) # + ggtitle(pop) +
ggsave(file=paste("../../results/summary_plots_stan/", "XPEHH_12_density", plout, sep=""))

ggplot() +
  geom_density(data=my_table_sim_stan, mapping=aes(x=XPEHH_23_snp), color=cone, fill=cone, alpha=0.2) +
  geom_density(data=my_table_app_stan, mapping=aes(x=XPEHH_23_snp), color="black", linetype="dashed") +
  guides(color=F) + xlab("XPEHH (EU-AS)") + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) # + ggtitle(pop) +
ggsave(file=paste("../../results/summary_plots_stan/", "XPEHH_23_density", plout, sep=""))

ggplot() +
  geom_density(data=my_table_sim_stan, mapping=aes(x=XPEHH_31_snp), color=cone, fill=cone, alpha=0.2) +
  geom_density(data=my_table_app_stan, mapping=aes(x=XPEHH_31_snp), color="black", linetype="dashed") +
  guides(color=F) + xlab("XPEHH (AS-AF)") + ylab("Density") + facet_wrap(~pop) + theme(aspect.ratio=1) # + ggtitle(pop) +
ggsave(file=paste("../../results/summary_plots_stan/", "XPEHH_31_density", plout, sep=""))
# }
my_table_sim_stan <- NULL
my_table_app_stan <- NULL
