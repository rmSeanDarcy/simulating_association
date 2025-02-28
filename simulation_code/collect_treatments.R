###########################################################################################
### SAMSARA - Collect all data from treatments in an experiment                         ###
###########################################################################################
### Read libraries
args <- commandArgs(trailingOnly=TRUE)
workdir <- args[1]
parent_set_of_analyses <- args[2]
### Read libraries

###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_collectfunctions'
setwd(workdir)
for (i in list.files(path= paste(workdir,c('/simulation_code/Load_collect_functions'),sep=''), full.names=TRUE)) {
  source(i)
}
### Finally, set the new working directory to the specific run you want to analyse
setwd(paste0(workdir,"/","Result_master_dir","/",parent_set_of_analyses))

###########################################################################################
### Load data from individual run folders
pd <- get_and_mod_data(parent_set_of_analyses)
write.csv(pd, 'exp_results.csv')

