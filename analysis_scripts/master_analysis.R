###########################################################################################
### SAMSARA - Result output analysis                                                    ###
###########################################################################################
library(dplyr)
library(vegan)
library(reshape2)

##############################################################################################################################
##### 1. Get passed arguments, load analysis functions and set working directory
##############################################################################################################################
setwd('/run/user/1000/gvfs/smb-share:domain=share,server=share.univie.ac.at,share=ter,user=seand93/PROJECTS/SomSOM/Step_III/SAMSARA')

###########################################################################################
### This script too was designed to run in a .bash script and takes the same first three inputs as master_slurm_simulation.py
#args <- commandArgs(TRUE)
# args <- commandArgs(trailingOnly=TRUE)
# experiment <- args[1]
# treatment <- args[2]
# simulation <- args[3]

###########################################################################################
### Inputs for manual runs 
experiment <- 'fig2_test'
treatment <- 's4'
simulation <- 's4'

###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_analysis_functions'
for (i in list.files(path=('./analysis_scripts/analysis_functions'),full.names=TRUE)) {
  source(i)
}
###########################################################################################
### Set working directory to simulation
setwd(paste0('./simulation_data/',experiment,'/',treatment,'/',simulation))

##############################################################################################################################
##### 2. Load and manipulate data
##############################################################################################################################
npop_cubd_log <- read.csv("./npop_cubd_log.csv")[,-1]
kabs_cubd_log <- read.csv("./kabs_cubd_log.csv")[,-1]
sp_correl_log <- read.csv("./sp_correl_log.csv")[,-1]
sign_ajd_res_log <- read.csv("./sign_ajd_res_log.csv")[,-1]
sprsc_cor_res_log <- read.csv("./sprsc_cor_res_log.csv")[,-1]
###########################################################################################
### Load simulated data
# See 'master_slurm_simulation.py' for a description of the content of files below
settings <- read.csv("./settings_log.csv")[,-1]
npop <- read.csv("./n_pop_log_log.csv")[,-1]
nxyz <- read.csv("./node_xyz_log.csv")[,-1]
kabs <- read.csv("./n_kabs_log.csv")[,-1]
eucdm <- read.csv("./euc_rsc_prf_sim_mat_log.csv")[,-1]
intm <- read.csv("./sp_int_mat_log.csv")[,-1]
###########################################################################################
### Extract information / Set parameters
cube_dim <- as.numeric(settings[settings$setting == 'shape',2])     # Dimensions of cube from 'shape' input
n_num <- as.numeric(settings[settings$setting == 'n_num',2])        # Number of habitats
sp_num <- as.numeric(settings[settings$setting == 'sp_num',2])      # Number of species
rsc_num <- as.numeric(settings[settings$setting == 'rsc_num',2])    # Number of resources
hab_nms <- unique(kabs$n_nms)           # Names of habitats
sp_nms <- colnames(npop)[1:sp_num]      # Names of habitats
rsc_nms <- colnames(kabs)[1:rsc_num]    # Names of habitats
#settings_updated <- read.csv("./settings_updated.csv")[,-1]
#ncubes_vec <- c(settings_updated[settings_updated$setting == 'sdim_vec',2])        # Number of habitats

###########################################################################################
### Delete erroneous runs -> See why in description of gLV integration in './Load_simulation_functions/set_dynamics.py'
err_runs <- npop[npop$err > 0,]
err_runs_idx <- unique(err_runs$repl)
#print(paste("The number of error containing runs was:",length(err_runs_idx)))
npop <- npop[!npop$repl %in% err_runs_idx,] 
nxyz <- nxyz[!nxyz$repl %in% err_runs_idx,]
kabs <- kabs[!kabs$repl %in% err_runs_idx,]
eucdm <- eucdm[!eucdm$repl %in% err_runs_idx,]
intm <- intm[!intm$repl %in% err_runs_idx,]
### Get all run indices without errors in runs
r_num <- unique(npop$repl) 
#r_num <- seq(0,1,1)

##############################################################################################################################
##### 3. Set additional parameters for analysis
##############################################################################################################################
### Figuring out the upper and lower quantiles of the average distribution of environmental preferences similarity E_{ij} 
# The 'get_Qs' function can be found in './Load_analysis_functions/miscelaneous_functions.R'
res <- get_Qs(eucdm, r_num)
euc_thr_pos <- res[1]                 # All pairs with E_{ij} above this value can be considered very similar
euc_thr_neg <- res[2]                 # All pairs with E_{ij} above this value can be considered very dissimilar
# These values give us the thresholds for determining whether species are very similar or dissimilar in their E_{ij}

##############################################################################################################################
##### 3. Set additional parameters for analysis
##############################################################################################################################

###########################################################################################
### Apply all analyses -> get_comp_res found in './Load_analysis_functions/main_functions.R' is the main function in which analysis scripts are called
comp_res <- analysis_habitat(sp_num, rsc_num, r_num, eucdm, ncubes_vec, euc_thr_neg, euc_thr_pos, npop_cubd_log, kabs_cubd_log, sp_correl_log, sign_ajd_res_log, sprsc_cor_res_log)
div_res_comp <- comp_res$div_res_log
div_rsc_res_comp <- comp_res$div_rsc_res_log
infmat_res_comp <- comp_res$infmat_res_log
full_res_subdir_comp <- comp_res$full_res_subdir_comp
infm_cooc_comp <- comp_res$infm_cooc_comp
###########################################################################################
### Save outputs (here still results per run)
write.csv(div_res_comp, "./div_res_comp.csv")
write.csv(div_rsc_res_comp, "./div_rsc_res_comp.csv")
write.csv(infmat_res_comp, "./infmat_res_comp.csv")
write.csv(full_res_subdir_comp, './full_res.csv')
write.csv(infm_cooc_comp, './infm_res.csv')
##############################################################################################################################
### The End
print("We done did it again!!!")
