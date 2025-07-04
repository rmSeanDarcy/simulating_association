###########################################################################################
### SAMSARA - Result analysis functions - Main analysis function                        ###
###########################################################################################

###########################################################################################
### Data analysis meta-function                                                         
process_composite <- function(npop_comp, kabs_comp, coords_comp, sdim_vec, sp_num, rsc_num, r_num, rand.time, npop, kabs, nxyz, nperm2, p_val_thr2, eucdm, ncubes_vec, addit_hab_smpl_num, pred_shuf, euc_thr_neg, euc_thr_pos, corcoeff_thr_neg, corcoeff_thr_pos) {
  ###########################################################################################
  ### The following data frames will be filled with results from respective analyses  
  ## Basic data
  npop_cubd_log <- c()            # All species abundance tables analysed
  kabs_cubd_log <- c()            # All resource abundance tables analysed
  ## Adjacency and co-occurrence data
  sign_ajd_res_log <- c()         # Log file for all significant and thresholded species x species correlations -> Co-occurrence adjacency matrices
  sp_correl_log <- c()            # -> All species x species correlations
  sprsc_cor_res_log <- c()        # Log file for all significant and thresholded species x resource correlations
  sprsc_correl_log <- c()         # -> All species x resource correlations
  ###########################################################################################
  ### Load sampling volume scales
  sdim_vec <- unique(npop_comp$d)
  ### Subsampling habitats
  hab_smpl_num <- unique(ncubes_vec)#[1:2]
  #h <- hab_smpl_num[1]
  #r <- 0
  for (h in hab_smpl_num) {
    npopd <- c()
    kabsd <- c()
    nxyzd <- c()
    for (r in r_num) {
      smpl_idx <- sample(nrow(npop[npop$repl == r,]))[1:h]
      npopr <- npop[npop$repl == r,]
      npopd <- as.data.frame(rbind(npopd, npopr[smpl_idx,]))
      kabsr <- kabs[kabs$repl == r,]
      kabsd <- as.data.frame(rbind(kabsd, kabsr[smpl_idx,]))
      nxyzr <- nxyz[nxyz$repl == r,]
      nxyzd <- as.data.frame(rbind(nxyzd, nxyzr[smpl_idx,]))
    }
    ### Add sampling noise here
    if (noise_std > 0) {
      npopd <- add_noise_per_habitat(npopd, noise_std, noise_mode)  
    }
    ### Calculation of significant associations
    sign_ajd_res <- get_signif_adj(npopd, p_val_thr2, nperm2, sp_num, r_num)
    sprsc_cor_res <- get_signif_sprsc(npopd, kabsd, p_val_thr2, nperm2, sp_num, rsc_num, r_num)
    sp_correl <- get_signif_adj_all(npopd, sp_num, r_num)
    ### Delete associations below a certain threshold (for now 0.7)
    sign_ajd_res[,1:(ncol(sign_ajd_res)-1)][sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] < corcoeff_thr_pos & sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] > corcoeff_thr_neg] <- 0
    sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)][sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] < corcoeff_thr_pos & sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] > corcoeff_thr_neg] <- 0
    ###########################################################################################
    ### Saving results into result logs (d is the identifier = 'HabSubsmpl_' + h the number of habitats subsampled)
    d_nm <- paste0('Subsmpl_',h,'__Dim_',0)
    npopd <- npopd[,colnames(npopd) %in% c(sp_nms, 'repl'),]
    npop_cubd_log <- rbind(npop_cubd_log, cbind(npopd, d=d_nm))
    kabsd <- kabsd[,colnames(kabsd) %in% c(rsc_nms, 'repl'),]
    kabs_cubd_log <- rbind(kabs_cubd_log, cbind(kabsd, d=d_nm))
    sp_correl_log <- rbind(sp_correl_log, cbind(sp_correl, d=d_nm))
    sign_ajd_res_log <- rbind(sign_ajd_res_log, cbind(sign_ajd_res, d=d_nm))
    sprsc_cor_res_log <- rbind(sprsc_cor_res_log, cbind(sprsc_cor_res, d=d_nm))
    print(paste0("Analysis complete for subsampling number: ",h))
  }
  ###########################################################################################
  ### Subsampling composite cubes
  # Here we subsample those cubes with composite communities. With the function 'get_scale_coms'
  # we previously aggregated all habitats within certain cube sizes (defined by side lengths in
  # sdim_vec) and saved resulting species abundances in npop_comp, resource abundances in kabs_comp
  # and new sampling cube coordinates in coords_comp (same corner of every cube). We can also subsample
  # these cubes in different numbers which we define in 'ncubes_vec'. Keep in mind that we are not only
  # limited by the total amount of cubes that can be sampled (27 cubes for cube side length = 0.33), 
  # but we also have to exclude those composites where no habitats are present. In most cases this 
  # is not that, but this can lead to issues with reduced test runs! 
  ### Gotta get the colnames into shape so I can merge the popdata with the following
  npop_cubd_log <- npop_cubd_log[,colnames(npop_cubd_log) %in% c(sp_nms, 'repl', 'd'),] 
  #hc <- ncubes_vec[1] 
  #d <- sdim_vec[1]
  #r <- 0
  for (hc in ncubes_vec) {
    for (d in sdim_vec) {
      npop_compd <- npop_comp[npop_comp$d == d,] 
      kabs_compd <- kabs_comp[kabs_comp$d == d,]
      coords_compd <- coords_comp[coords_comp$d == d,]
      ### Remove those runs where a required number of samples (cubes) are not present:
      npopd <- c()
      kabsd <- c()
      nxyzd <- c()
      for (r in r_num) {
        ### If there are fewer samples than cubes:
        cubidx <- sort(sample(npop_compd[npop_compd$repl == r,]$smpl_cub_idx, hc))
        npopd <- as.data.frame(rbind(npopd, npop_compd[npop_compd$repl == r & npop_compd$smpl_cub_idx %in% cubidx,]))
        kabsd <- as.data.frame(rbind(kabsd, kabs_compd[kabs_compd$repl == r & kabs_compd$smpl_cub_idx %in% cubidx,]))
        nxyzd <- as.data.frame(rbind(nxyzd, coords_compd[coords_compd$repl == r & coords_compd$smpl_cubs_idx %in% cubidx,]))
      }
      ###########################################################################################
      if (noise_std > 0) {
        colnames(npopd)[colnames(npopd) == 'smpl_cub_idx'] <- 'Ridx'
        npopd <- add_noise_per_habitat(npopd, noise_std, noise_mode)
        colnames(npopd)[colnames(npopd) == 'Ridx'] <- 'smpl_cub_idx'
      }
      ### Calculation of significant associations
      sign_ajd_res <- get_signif_adj(npopd, p_val_thr2, nperm2, sp_num, r_num)
      sprsc_cor_res <- get_signif_sprsc(npopd, kabsd, p_val_thr2, nperm2, sp_num, rsc_num, r_num)
      sp_correl <- get_signif_adj_all(npopd, sp_num, r_num)
      ### Delete associations below a certain threshold (for now 0.7)
      sign_ajd_res[,1:(ncol(sign_ajd_res)-1)][sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] < corcoeff_thr_pos & sign_ajd_res[,1:(ncol(sign_ajd_res)-1)] > corcoeff_thr_neg] <- 0
      sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)][sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] < corcoeff_thr_pos & sprsc_cor_res[,1:(ncol(sprsc_cor_res)-1)] > corcoeff_thr_neg] <- 0
      ### Saving all results, individually and combined
      dc_nm <- paste0('Subsmpl_',hc,'__Dim_', round(d,3))
      npopd <- npopd[,colnames(npopd) %in% c(sp_nms, 'repl'),]
      npop_cubd_log <- rbind(npop_cubd_log, cbind(npopd, d=dc_nm))
      kabsd <- kabsd[,colnames(kabsd) %in% c(rsc_nms, 'repl'),]
      kabs_cubd_log <- rbind(kabs_cubd_log, cbind(kabsd, d=dc_nm))
      sp_correl_log <- rbind(sp_correl_log, cbind(sp_correl, d=dc_nm))
      sign_ajd_res_log <- rbind(sign_ajd_res_log, cbind(sign_ajd_res,d=dc_nm))
      sprsc_cor_res_log <- rbind(sprsc_cor_res_log, cbind(sprsc_cor_res,d=dc_nm))
      print(paste0("Analysis complete for subsampling number: ",hc," at sample cube side length: ",d))
    }
  }
  ###########################################################################################
  ### Export results 
  comp_res <- list(npop_cubd_log = npop_cubd_log,
                   kabs_cubd_log = kabs_cubd_log,
                   sp_correl_log = sp_correl_log,
                   sign_ajd_res_log = sign_ajd_res_log,
                   sprsc_cor_res_log = sprsc_cor_res_log)
  return(comp_res)
}
