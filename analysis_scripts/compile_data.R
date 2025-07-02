###########################################################################################
### SAMSARA - Collect all data from treatments in an experiment                         ###
###########################################################################################
### Read libraries
rm(workdir)  # Remove any previous definition of 'workdir'

args <- commandArgs(trailingOnly=TRUE)
workdir <- args[1]
parent_set_of_analyses <- args[2]

###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_collectfunctions'
source(paste0(workdir,c('/simulation_code/Load_collect_functions/collectfuns.R')))
### Finally, set the new working directory to the specific run you want to analyse
setwd(paste0(workdir,"/Result_master_dir/",parent_set_of_analyses))

###########################################################################################
### Load data from individual run folders
pd <- get_and_mod_data(parent_set_of_analyses)
write.csv(pd, 'exp_results.csv')

