###########################################################################################
### SAMSARA - Result analysis functions - Main analysis function                        ###
###########################################################################################

###########################################################################################
### Data analysis meta-function                                                         
get_comp_res <- function(sp_num, rsc_num, r_num, npop, kabs, nxyz, nperm2, p_val_thr2, eucdm, ncubes_vec, euc_thr_neg, euc_thr_pos, corcoeff_thr_neg, corcoeff_thr_pos) {
  
  ### The following data frames will be filled with results from respective analyses  
  ## Basic data
  npop_cubd_log <- c()            # All species abundance tables analysed
  kabs_cubd_log <- c()            # All resource abundance tables analysed
  
  div_res_log <- c()              # Metrics on community data (richness, beta diversity etc.)
  div_rsc_res_log <- c()          # Metrics on resource data (richness, beta diversity etc.)
  
  sign_ajd_res_log <- c()         # Log file for all significant and thresholded species x species correlations -> Co-occurrence adjacency matrices
  sp_correl_log <- c()            # -> All species x species correlations
  sprsc_cor_res_log <- c()        # Log file for all significant and thresholded species x resource correlations
  sprsc_correl_log <- c()         # -> All species x resource correlations
  
  ## Co-occurrence data  
  infmat_res_log <- c()           # Big log file with all species pairwise combinations and relevant metrics 
  # -> co-occurrence, interaction coefficients, environmental preference similarity, environmental triplets, ...
  infmat_num_res_log <- c()       # Co-occurrences quantities -> Number of positive/negative co-occurrences and positive/negative triplets    
  
  ###########################################################################################
  ### Subsampling habitats and analysing data
  hab_smpl_num <- ncubes_vec 
  ### Cycle through the specified number of habitats of interest (in the paper only 25, apart from the testing where we analyse 10-300 in Supporting information)
  #h <- hab_smpl_num[1]
  for (h in hab_smpl_num) {
    npopd <- c()
    kabsd <- c()
    nxyzd <- c()
    for (r in r_num) {
      ### Here we select a random number of habitats, and subsample every dataset by their indices
      smpl_idx <- sample(nrow(npop[npop$repl == r,]))[1:h]
      npopr <- npop[npop$repl == r,]
      npopd <- as.data.frame(rbind(npopd, npopr[smpl_idx,]))
      kabsr <- kabs[kabs$repl == r,]
      kabsd <- as.data.frame(rbind(kabsd, kabsr[smpl_idx,]))
      nxyzr <- nxyz[nxyz$repl == r,]
      nxyzd <- as.data.frame(rbind(nxyzd, nxyzr[smpl_idx,]))
    }
    
    ###########################################################################################
    ### Main analysis functions
    ###########################################################################################
    ### The R scripts containing each function is added after the #
    # All functions can be found within './Load_analysis_functions'
    
    ## Community and resource diversity metrics
    div_res <- div_mtrcs(npopd, sp_num, hab_num = h, r_num, eucdm)          # diversity_functions.R
    div_rsc_res <- div_mtrcs_rsc(kabsd, rsc_num, hab_num = h, r_num)        # diversity_functions.R

    ## Simple adjacency matrices (Spearmans correlations)
    sp_correl <- get_signif_adj_all(npop, sp_num, r_num)                                            # calculate_coocurrence_functions.R
    ## Calculation of significant associations
    sign_ajd_res <- get_signif_adj(npopd, p_val_thr2, nperm2, sp_num, r_num)                        # calculate_coocurrence_functions.R
    sprsc_cor_res <- get_signif_sprsc(npopd, kabsd, p_val_thr2, nperm2, sp_num, rsc_num, r_num)     # calculate_coocurrence_functions.R
    
    ## Delete associations below defined thresholds (0.7 and -0.7 in paper)
    sign_ajd_res[,1:(ncol(sign_ajd_res)-1)][sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] < corcoeff_thr_pos & sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] > corcoeff_thr_neg] <- 0
    sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)][sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] < corcoeff_thr_pos & sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] > corcoeff_thr_neg] <- 0
    
    ## Prediction metrics
    infmat_res <- get_pred_stats(sp_correl, sign_ajd_res, sp_num, r_num, eucdm, intm, sprsc_cor_res, rsc_num, sp_nms, pred_shuf, euc_thr_neg, euc_thr_pos)  # collecting_coocurrence_functions.R

    ###########################################################################################
    ### Save result data
    ## Saving results into result logs (d is the identifier = 'HabSubsmpl_' + h the number of habitats subsampled)
    d_nm <- paste0('HabSubsmpl_',h)
    # Find definitions of datasets above
    npop_cubd_log <- rbind(npop_cubd_log, cbind(npopd, d=d_nm))
    kabs_cubd_log <- rbind(kabs_cubd_log, cbind(kabsd, d=d_nm))
    
    div_res_log <- rbind(div_res_log, cbind(div_res, d=d_nm))
    div_rsc_res_log <- rbind(div_rsc_res_log, cbind(div_rsc_res, d=d_nm))
    
    sp_correl_log <- rbind(sp_correl_log, cbind(sp_correl, d=d_nm))
    sign_ajd_res_log <- rbind(sign_ajd_res_log, cbind(sign_ajd_res, d=d_nm))
    sprsc_cor_res_log <- rbind(sprsc_cor_res_log, cbind(sprsc_cor_res, d=d_nm))
    
    infmat_res_log <- rbind(infmat_res_log, cbind(infmat_res, d=d_nm))
    print(paste0("Analysis complete for subsampling number: ",h))
  }
  
  ###########################################################################################
  ## Next we calculate the means of basic diversity metrics
  # Functions can be found in './Load_analysis_functions/means_for_all_runs_functions.R'
  full_res_subdir_comp <- get_means_of_all_results(div_res_log,div_rsc_res_log)
  ## First we quantify how well the co-occurrences match drivers accross all runs 
  # Functions can be found in './Load_analysis_functions/matching_coocurrence_functions.R'
  infm <- infmat_res_log
  infm_cooc_comp <- get_cooc_imfmpreds(infm, euc_thr_pos, euc_thr_neg)
  
  ###########################################################################################
  ### Export results 
  comp_res <- list(npop_cubd_log = npop_cubd_log,
                   kabs_cubd_log = kabs_cubd_log,
                   div_res_log = div_res_log,
                   div_rsc_res_log = div_rsc_res_log,
                   sp_correl_log = sp_correl_log,
                   sign_ajd_res_log = sign_ajd_res_log,
                   sprsc_cor_res_log = sprsc_cor_res_log,
                   infmat_res_log = infmat_res_log,
                   full_res_subdir_comp = full_res_subdir_comp,
                   infm_cooc_comp = infm_cooc_comp)
  return(comp_res)
}




