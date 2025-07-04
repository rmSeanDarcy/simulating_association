###########################################################################################
### SAMSARA - Result analysis functions - Main analysis function                        ###
###########################################################################################

###########################################################################################
### Data analysis meta-function                                                         
analysis_habitat <- function(sp_num, rsc_num, r_num, eucdm, ncubes_vec, euc_thr_neg, euc_thr_pos, npop_cubd_log, kabs_cubd_log, sp_correl_log, sign_ajd_res_log, sprsc_cor_res_log) {
  ###########################################################################################
  ### The following data frames will be filled with results from respective analyses  
  ## Basic data
  div_res_log <- c()
  div_rsc_res_log <- c()
  infmat_res_log <- c()
  ###########################################################################################
  ### Subsampling habitats and analysing data
  #hab_smpl_num <- ncubes_vec 
  ### Cycle through the specified number of habitats of interest (in the paper only 25, apart from the testing where we analyse 10-300 in Supporting information)
  #d <- unique(npop_cubd_log$d)[1]
  for (d in unique(npop_cubd_log$d)) {
    npopd <- npop_cubd_log[npop_cubd_log$d == d,]
    kabsd <- kabs_cubd_log[kabs_cubd_log$d == d,]
    sp_correl <- sp_correl_log[sp_correl_log$d == d,]
    sign_ajd_res <- sign_ajd_res_log[sign_ajd_res_log$d == d,]
    sprsc_cor_res <- sprsc_cor_res_log[sprsc_cor_res_log$d == d,] 
    ###########################################################################################
    ### Main analysis functions
    # Diversity
    div_res <- div_mtrcs(npopd, sp_num, hab_num = 25, r_num, eucdm)
    div_rsc_res <- div_mtrcs_rsc(kabsd, rsc_num, hab_num = 25, r_num)
    ### Prediction metrics
    infmat_res <- get_pred_stats(sp_correl, sign_ajd_res, sp_num, r_num, eucdm, intm, sprsc_cor_res, rsc_num, sp_nms, pred_shuf, euc_thr_neg, euc_thr_pos) 
    ### Save result data
    ## Saving results into result logs (d is the identifier = 'HabSubsmpl_' + h the number of habitats subsampled)
    d_nm <- d
    div_res_log <- rbind(div_res_log, cbind(div_res, d=d_nm))
    div_rsc_res_log <- rbind(div_rsc_res_log, cbind(div_rsc_res, d=d_nm))
    infmat_res_log <- rbind(infmat_res_log, cbind(infmat_res, d=d_nm))
    print(paste0("Analysis complete for subsampling number: ",d_nm))
  }
  ###########################################################################################
  ## Next we calculate the means of basic diversity metrics
  # Functions can be found in './Load_analysis_functions/means_for_all_runs_functions.R'
  full_res_subdir_comp <- get_means_of_all_results(div_res_log,div_rsc_res_log)
  ## First we quantify how well the co-occurrences match drivers accross all runs 
  # Functions can be found in './Load_analysis_functions/matching_coocurrence_functions.R'
  infm <- infmat_res_log
  infm_cooc_comp <- get_cooc_imfmpreds(infmat_res_log, euc_thr_pos, euc_thr_neg)
  ###########################################################################################
  ### Export results 
  comp_res <- list(div_res_log = div_res_log,
                   div_rsc_res_log = div_rsc_res_log,
                   infmat_res_log = infmat_res_log,
                   full_res_subdir_comp = full_res_subdir_comp,
                   infm_cooc_comp = infm_cooc_comp)
  return(comp_res)
}




