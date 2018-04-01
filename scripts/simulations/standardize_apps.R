library(tidyverse)
label = "1000GPP1"
pop_list = c("YRI", "CEU", "CHB")

for (pop in pop_list) {
  # load features and metadata
  if (exists("my_metadata_apps_stan")) {rm(my_metadata_apps_stan)}
  if (exists("my_features_apps_stan")) {rm(my_features_apps_stan)}
  if (exists("my_table_apps_stan")) {rm(my_table_apps_stan)}
  for (chr in 1:22) {
    if (!exists("my_metadata_apps_stan")) {
      my_metadata_apps_stan <- read_tsv(paste("../../data/applications/1000GP_Phase1_stan/1000GPP1_chr", chr, "_pop", pop, "_metadata.txt", sep=""))
    }
    else {
      my_metadata_apps_stan <- read_tsv(paste("../../data/applications/1000GP_Phase1_stan/1000GPP1_chr", chr, "_pop", pop, "_metadata.txt", sep="")) %>% full_join(my_metadata_apps_stan)
    }
    if (!exists("my_features_apps_stan")) {
      my_features_apps_stan <- read_tsv(paste("../../data/applications/1000GP_Phase1_stan/1000GPP1_chr", chr, "_pop", pop, "_features.txt", sep=""))
    }
    else {
      my_features_apps_stan <- read_tsv(paste("../../data/applications/1000GP_Phase1_stan/1000GPP1_chr", chr, "_pop", pop, "_features.txt", sep="")) %>% full_join(my_features_apps_stan)
    }
  }
  
  # join metadata to features
  my_table_apps_stan <- inner_join(my_metadata_apps_stan, my_features_apps_stan, by="key")
  
  # add DAF_bins
  bins = 40
  my_table_apps_stan <- my_table_apps_stan %>% mutate(DAF_bins = cut(DAF_snp, breaks=bins, labels=c(1:bins)))
  table(my_table_apps_stan$DAF_bins); summarise(group_by(my_table_apps_stan, DAF_bins), min(DAF_snp), max(DAF_snp))
  
  # get means and sd
  my_applications_iHH <- 
    my_table_apps_stan %>% group_by(DAF_bins) %>% summarise(
      iHS_snp_mean=round(mean(iHS_snp, na.rm=TRUE), 4), 
      iHS_snp_sd=round(sd(iHS_snp, na.rm=TRUE), 4) # ,
      # D_iHH_snp_mean=round(mean(D_iHH_snp, na.rm=TRUE), 4), 
      # D_iHH_snp_sd=round(sd(D_iHH_snp, na.rm=TRUE), 4)
    )
  
  # save means and sd
  write.table(mutate(my_applications_iHH, DAF_bins=as.numeric(DAF_bins)), paste("../../scripts/simulations/standardize/my_applications_iHH_", pop, ".txt", sep=""), sep='\t', row.names=FALSE)
  
  # get means and sds
  my_applications_XPEHH <- 
    my_table_apps_stan %>% summarise(
      XPEHH_12_snp_mean=round(mean(XPEHH_12_snp, na.rm=TRUE), 4), 
      XPEHH_12_snp_sd=round(sd(XPEHH_12_snp, na.rm=TRUE), 4),
      XPEHH_23_snp_mean=round(mean(XPEHH_23_snp, na.rm=TRUE), 4), 
      XPEHH_23_snp_sd=round(sd(XPEHH_23_snp, na.rm=TRUE), 4), 
      XPEHH_31_snp_mean=round(mean(XPEHH_31_snp, na.rm=TRUE), 4), 
      XPEHH_31_snp_sd=round(sd(XPEHH_31_snp, na.rm=TRUE), 4)
    )
  
  # save means and sds
  write.table(my_applications_XPEHH, paste("../../scripts/simulations/standardize/my_applications_XPEHH_", pop, ".txt", sep=""), sep='\t', row.names=FALSE)
}
