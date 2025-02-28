###########################################################################################
### SAMSARA - Loading Set of analyses data and preparing for plotting                   ###
###########################################################################################
#workdir = '/home/swani/Documents/computational_research_tools/homework4/SAMSARA'
#setwd(workdir)
#parent_set_of_analyses = 'Fig2'

###########################################################################################
# Load and prep data
get_and_mod_data <- function(parent_set_of_analyses) {
  runs_incl <- list.dirs(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses), full.names = FALSE, recursive = FALSE)
  # For every treatment rbind results to all results
  pd <- c()
  #i <- 1
  for (i in 1:length(runs_incl)) {
    pred_res <- load_data_fun(parent_set_of_analyses, set_of_analyses = runs_incl[i])
    pd <- rbind(pd,pred_res)
    print(runs_incl[i])
  }
  ###########################################################################################
  ### Add corrected match values
  # All comparisons considered
  pd$eucp_cor <- pd$eucp - (pd$eucshrp - pd$eucshrintp) - (pd$eucintp - pd$eucshrintp) - pd$eucshrintp
  pd$shrp_cor <- pd$shrp - (pd$eucshrp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$intp_cor <- pd$intp - (pd$eucintp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$eucshrp_cor <- pd$eucshrp - pd$eucshrintp
  pd$eucintp_cor <- pd$eucintp - pd$eucshrintp
  pd$shrintp_cor <- pd$shrintp - pd$eucshrintp
  #
  pd$eucn_cor <- pd$eucn - (pd$eucshrn - pd$eucshrintn) - (pd$eucintn - pd$eucshrintn) - pd$eucshrintn
  pd$shrn_cor <- pd$shrn - (pd$eucshrn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$intn_cor <- pd$intn - (pd$eucintn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$eucshrn_cor <- pd$eucshrn - pd$eucshrintn
  pd$eucintn_cor <- pd$eucintn - pd$eucshrintn
  pd$shrintn_cor <- pd$shrintn - pd$eucshrintn
  #####
  # A Resource preference and interactions based matching considered
  pd$eucp_corA <- pd$eucp - pd$eucintp
  pd$intp_corA <- pd$intp - pd$eucintp
  pd$eucintp_corA <- pd$eucintp
  pd$unxpA <- pd$unxp + (pd$shrp - (pd$eucshrp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corA <- pd$eucn - pd$eucintn
  pd$intn_corA <- pd$intn - pd$eucintn
  pd$eucintn_corA <- pd$eucintn
  pd$unxnA <- pd$unxn + (pd$shrn - (pd$eucshrn) - (pd$shrintn - pd$eucshrintn))
  #####
  # B Resource preference and shared resource association based matching considered
  pd$eucp_corB <- pd$eucp - pd$eucshrp
  pd$shrp_corB <- pd$shrp - pd$eucshrp
  pd$eucshrp_corB <- pd$eucshrp
  pd$unxpB <- pd$unxp + (pd$intp - (pd$eucintp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corB <- pd$eucn - pd$eucshrn
  pd$shrn_corB <- pd$shrn - pd$eucshrn
  pd$eucshrn_corB <- pd$eucshrn
  pd$unxnB <- pd$unxn + (pd$intn - (pd$eucintn) - (pd$shrintn - pd$eucshrintn))
  return(pd)
}

###########################################################################################
load_data_fun <- function(parent_set_of_analyses, set_of_analyses) {
  parentdirs <- list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,"/",set_of_analyses,sep=''), full.names=TRUE) #, pattern = "Z"
  pred_res <- read.csv(paste0(parentdirs,'/infm_res.csv'))[-1]
  pred_res$run_nm <- set_of_analyses
  return(pred_res)
}
