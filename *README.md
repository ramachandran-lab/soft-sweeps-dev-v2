# Soft but not hard sweep classification is sensitive to unknown class imbalance in Bayesian classifier

Stephen Rong, Sohini Ramachandran, Brown University

Questions: Contact stephen[underscore]rong[at]brown[dot]edu or post a git issue

## Contents:

This git repo contains the following directories:

1. /scripts/simulations/ scripts for performing the simulations and calculating summary statistics, contains a copy of ...

2. /scripts/applications/ scripts for preprocessing the 1000 Genomes Phase 1 raw haplotype data

3. /scripts/analyses/ scripts for generating temporary and final figures and results

4. /results/ used to contain the generated temporary and final figures and results 

## Software Requirements:

Python v2.7.13 -- R v3.4.1 -- 

<!-- 5. /data/applications/1000GP_Phase1_raw/ contains the 1000 Genomes Phase 1 raw haplotype data, downloaded from ...

6. /1000GP_Phase1_raw/ contains the preprocessed 1000 Genomes Phase 1 haplotypes ...

7. 1000GP_Phase1_misc/ contains extra files used in the analysis, explained in README.txt ...

8. 1000GP_Phase1_stan/ contains the SNP-level summary statistics used to standardize iHS, DiHH, and XP-EHH for the application data ...

9. 1000GP_Phase1_main/ contains the window-level summary statistics used as the application data ...

10. /hapbin_out_temp/ are used for intermediate hapbin files ...

11. /data/simulations/2017111501/ contains the SNP-level summary statistics used to standardize iHS, DiHH, and XP-EHH for the simulation data ...

12. /2017111502/ contains the window-level summary statistics split into two disjoing sets for the training data and the test data

13. /2017111502/ contains the window-level summary statistics used to assess performance on neutral sits linked to a beneficial mutation ...

14. /hapbin_in_temp/ and /hapbin_out_temp/ are used for intermediate hapbin files ... -->

<!-- ## Software Requirements:

Python v2.7.13

numpy v1.11.3

pandas v0.19.2

pygg 0.1.7

networkx v1.11

matplotlib v2.0.0

scipy v0.18.1

sklearn v0.18.1

R v3.4.0

RColorBrewer v1.1.2

tidyverse v1.1.1

wrapr v0.1.2

ks v1.10.6

gplots v3.0.1

ggpubr v0.1.2

scales v0.4.1

plyr v1.8.4

msms (https://github.com/delt0r/msms)

hapbin (https://github.com/evotools/hapbin) -->

## Generating Training, Test, and Application Data:

Download 1000 Genomes Phase 1 raw haplotypes from [here](https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono.tgz) (accessed on 12/8/2016), and unzip contents to /data/applictions/1000GP_Phase1_raw/. To generate contents of /data/applications/1000GP_Phase1_post/, go to /scripts/applications/ and run preprocess_data.py. Make sure modules in /scripts/simulations/ are in path and hapbin is installed. 

To generate contents of /data/applications/1000GP_Phase1_stan/, go to /scripts/simulations/slurm_scripts/ and run scripts_apps_stan.py, and then run 1000GPP1_stan.sh from /scripts/simulations/. To generate contents of /data/simulations/2017111501/, go to /scripts/simulations/slurm_scripts/ and run scripts_sims_stan.py, and then run 2017111501.sh from /scripts/simulations/. 

Then go to /scripts/analyses/ and run standardize_apps.R and standardize_sims.R. 

To generate contents of /data/applications/1000GP_Phase1_main/, go to /scripts/simulations/slurm_scripts/ and run scripts_apps_main.py, and then run 1000GPP1_main.sh from /scripts/simulations/. To generate contents of /data/simulations/2017111502 and /data/simulations/2017111502, go to /scripts/simulations/slurm_scripts/ and run scripts_sims_train.py and scripts_sims_test.py, and then run 2017111502.sh and 2017111503.sh from /scripts/simulations/.

To generate result figures and tables, go to the /scripts/analyses/ directory and run master_script.R from within R, with the first two lines uncommented.

## Generating Result Figures from Saved Above Data:

<!-- We have included the training and test data required for generating result figures and tables as temp.zip.  -->

To generate result figures and tables, unzip these folders and replace the existing placeholder folders. Then go to the /scripts/analyses/ directory and run master_script.R from within R.
