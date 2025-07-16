###########################################################################################
### SAMSARA - Result output analysis                                                    ###
###########################################################################################
library(dplyr)
library(reshape2)
library(vegan)

##############################################################################################################################
##### 1. Get passed arguments, load analysis functions and set working directory
##############################################################################################################################
#setwd('/run/user/1000/gvfs/smb-share:domain=share,server=share.univie.ac.at,share=ter,user=seand93/PROJECTS/SomSOM/Step_III/SAMSARA')
#setwd('/mnt/share/PROJECTS/SomSOM/Step_III/SAMSARA')

###########################################################################################
### This script too was designed to run in a .bash script and takes the same first three inputs as master_slurm_simulation.py
#args <- commandArgs(TRUE)
args <- commandArgs(trailingOnly=TRUE)
experiment <- args[1]
treatment <- args[2]
simulation <- args[3]
# Extract the remaining arguments as a single string
paramstr <- ifelse(length(args) > 3, args[4], "")
#paramstr <- "--noise=1"      # test

### Now we need to first check whether noise is in the additional arguments and if so we need to extract the value 
#paramstr <- '--n_num:25 --nrep:2 --noise:0.1'
if (grepl('noise', paramstr)) {
  noise_std <- as.numeric(strsplit(strsplit(paramstr,' --')[[1]][grep('noise=',(strsplit(paramstr,' --')[[1]]))],'noise=')[[1]][2])
} else {
  noise_std <- 0
}
###########################################################################################
### Inputs for manual runs 
#experiment <- 'xfig2'
#treatment <- 's1'
#simulation <- 's1'
#noise_std <- 0.25

###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_analysis_functions'
for (i in list.files(path=('./analysis_scripts/processing_functions'),full.names=TRUE)) {
  source(i)
}

###########################################################################################
### Set working directory to simulation
setwd(paste0('./simulation_data/',experiment,'/',treatment,'/',simulation))


##############################################################################################################################
##### 2. Load and manipulate data
##############################################################################################################################

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
#noise_std <- 0                          # as.numeric(settings[settings$setting == 'rsc_num',2])    # Number of resources
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
#r_num <- seq(0,10,1)

###########################################################################################
# In some very rare cases species abundances are fixed at 0.001 for the entire simulation. This is the case when a species is in an
# unfavourable habitat (meaning K_{i,h} < 0.1) and has no interaction partners. As all species are initialised with an abundance 
# of 0.001 (ini_abd) and these species have a growth rate of 0 they remain at this abundance throughout the simulation
ini_abd <- 0.001
npop_copy <- npop[,1:sp_num]
npop_copy[npop_copy <= ini_abd] <- 0
npop <- cbind(npop_copy, npop[,c(sp_num+1:4)])
###########################################################################################
### Adding noise to population data
# There are two options, but only the abundance adjusted noise (noise_mode <- 'habitat_adj') is presented in the paper!
noise_mode <- 'habitat_adj'           # 'fixed_val_norm', 'habitat_adj'
# # The value for noise_std is given above in 1. -> The function that adds noise can be found in './Load_analysis_functions/miscelaneous_functions.py'
# if (noise_std > 0) {
#   npop <- add_noise_per_habitat(npop, noise_std, noise_mode)  
# }


##############################################################################################################################
##### 3. Set additional parameters for analysis
##############################################################################################################################

###########################################################################################
### Settings for analysis functions
sdim_vec <- c(0.333)                  # This value is not needed in this script! It becomes relevant in the 'master_slurm_analysis_cubes.R' script
ncubes_vec <- c(25)                   # Number of habitats that are randomly sampled -> Needs to be equal to or lower than n_num. 
# Here either a single number or a vector can be provided. Results will always have a column 'Hab_subsmplX' with X being the number of habitats sampled 
###########################################################################################
### Figuring out the upper and lower quantiles of the average distribution of environmental preferences similarity E_{ij} 
# The 'get_Qs' function can be found in './Load_analysis_functions/miscelaneous_functions.R'
res <- get_Qs(eucdm, r_num)
euc_thr_pos <- res[1]                 # All pairs with E_{ij} above this value can be considered very similar
euc_thr_neg <- res[2]                 # All pairs with E_{ij} above this value can be considered very dissimilar
# These values give us the thresholds for determining whether species are very similar or dissimilar in their E_{ij}

###########################################################################################
### Set correlation coefficient threshold and parameters for null model testing
corcoeff_thr_neg <- -0.7              # Spearmans correlation coefficient threshold for negative co-occurrence
corcoeff_thr_pos <- 0.7               # Spearmans correlation coefficient threshold for positive co-occurrence
nperm2 <- 1000                        # Number of Nullmodels generated for calculating sign. co-occurrences
p_val_thr2 <- 0.05                    # Significance level for which sign. co-occurrences are tests

###########################################################################################
### Update settings with those from analysis
settings2 <- as.data.frame(cbind(c('sdim_vec','ncubes_vec','nperm2','p_val_thr2','noise_mode','noise_std','euc_thr_pos','euc_thr_neg','corcoeff_thr_neg','corcoeff_thr_pos'),
                                 c(paste(round(sdim_vec,3), collapse = ' '),paste(ncubes_vec, collapse = ' '),nperm2,p_val_thr2,noise_mode,noise_std,euc_thr_pos,euc_thr_neg,corcoeff_thr_neg,corcoeff_thr_pos)))
colnames(settings2) <- colnames(settings)
settings <- rbind(settings,settings2)
write.csv(settings,'settings_updated.csv')


##############################################################################################################################
##### 3. Set additional parameters for analysis
##############################################################################################################################

###########################################################################################
### Apply all analyses -> get_comp_res found in './Load_analysis_functions/main_functions.R' is the main function in which analysis scripts are called
comp_res <- process_habitat(sp_num, rsc_num, r_num, npop, kabs, nxyz, nperm2, p_val_thr2, eucdm, ncubes_vec, euc_thr_neg, euc_thr_pos, corcoeff_thr_neg, corcoeff_thr_pos)
npop_cubd_log <- comp_res$npop_cubd_log
kabs_cubd_log <- comp_res$kabs_cubd_log
sp_correl_log <- comp_res$sp_correl_log
sign_ajd_res_log <- comp_res$sign_ajd_res_log
sprsc_cor_res_log <- comp_res$sprsc_cor_res_log
###########################################################################################
### Save outputs (here still results per run)
write.csv(npop_cubd_log, "./npop_cubd_log.csv")
write.csv(kabs_cubd_log, "./kabs_cubd_log.csv")
write.csv(sp_correl_log, "./sp_correl_log.csv")
write.csv(sign_ajd_res_log, "./sign_ajd_res_log.csv")
write.csv(sprsc_cor_res_log, "./sprsc_cor_res_log.csv")
##############################################################################################################################
### The End
print("We done did it!!!")
