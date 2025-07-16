###########################################################################################
### SAMSARA - Collect all data from treatments and simulations in an experiment         ###
###########################################################################################

##############################################################################################################################
##### 1. Get passed arguments, load analysis functions and set working directory
##############################################################################################################################
#setwd('/run/user/1000/gvfs/smb-share:domain=share,server=share.univie.ac.at,share=ter,user=seand93/PROJECTS/SomSOM/Step_III/SAMSARA')

###########################################################################################
args <- commandArgs(trailingOnly=TRUE)
experiment <- args[1]
#experiment <- 'fig4'
dir.create(file.path(paste0('./analysis_data/',experiment)), showWarnings = FALSE)


###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_analysis_functions'
for (i in list.files(path=('./analysis_scripts/compile_functions'),full.names=TRUE)) {
  source(i)
}

###########################################################################################
### Collect data for all treatments and simulations within an experiment
infm_res_compiled <- compile_data(experiment_name = experiment, "infm_res.csv")
full_res_compiled <- compile_data(experiment_name = experiment, "full_res.csv")
###########################################################################################
### Load data from individual run folders
write.csv(infm_res_compiled, paste0('./analysis_data/', experiment, '/infm_res_compiled.csv'))
write.csv(full_res_compiled, paste0('./analysis_data/', experiment, '/full_res_compiled.csv'))
