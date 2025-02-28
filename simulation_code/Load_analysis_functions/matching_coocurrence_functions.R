###########################################################################################
### SAMSARA - Functions for quantifying matches between co-occurrences and drivers      ###
###########################################################################################

### Get number of co-occurrences predicted by different drivers
get_cooc_imfmpreds <- function(infm, euc_thr_pos, euc_thr_neg) {
  
  resl <- c()
  ## We cycle through each number of habitats sampled (d)
  #j <- unique(infm$d)[1]
  for (j in unique(infm$d)) {
    infmj <- infm[infm$d == j,]
    res <- c()
    ## Now we assess matches in each individual run 
    #r <- 0 
    for (r in unique(infm$repl)) {
      infmr <- infmj[infmj$repl == r,]
      
      ###########################################################################################
      ### Get positive matching data
      ###########################################################################################
      infmrp <- infmr[infmr$coocpos == 1,]
      totp <- nrow(infmrp)
      eucp <- (sum(infmrp$euc < euc_thr_pos))/totp
      intp <- (sum(infmrp$intpos == 1))/totp
      shrp <- (sum(infmrp$shrd_rsc_pos == 1))/totp
      eucshrp <- (sum(infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1))/totp
      eucintp <- (sum(infmrp$euc < euc_thr_pos & infmrp$intpos == 1))/totp
      shrintp <- (sum(infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1))/totp
      eucshrintp <- (sum(infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1))/totp
      unxp <- sum(infmrp$euc > euc_thr_pos & infmrp$intpos == 0 & infmrp$shrd_rsc_pos == 0)/totp
      eucp_cor <- eucp - (eucshrp - eucshrintp) - (eucintp - eucshrintp) - eucshrintp
      shrp_cor <- shrp - (eucshrp - eucshrintp) - (shrintp - eucshrintp) - eucshrintp
      intp_cor <- intp - (eucintp - eucshrintp) - (shrintp - eucshrintp) - eucshrintp
      eucshrp_cor <- eucshrp - eucshrintp
      eucintp_cor <- eucintp - eucshrintp
      shrintp_cor <- shrintp - eucshrintp
      print(paste('In', j, 'rep', r, 'matching positvives in nrep',r,'sum to:', c(eucp_cor + shrp_cor + intp_cor + eucshrp_cor + eucintp_cor + shrintp_cor + eucshrintp + unxp))) 
      ###########################################################################################
      ## Get corresponding data on mean interaction coefficients
      # All below is not used in the paper
      mnint_intp <- mean(infmrp[infmrp$intpos == 1,]$int)
      mnint_eucp <- mean(infmrp[infmrp$euc < euc_thr_pos,]$int)
      mnint_shrp <- mean(infmrp[infmrp$shrd_rsc_pos == 1,]$int)
      mnint_eucshrp <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1,]$int)
      mnint_eucintp <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$intpos == 1,]$int)
      mnint_shrintp <- mean(infmrp[infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1,]$int)
      ## Get corresponding data on mean environmental preference similairty  
      mneuc_intp <- mean(abs(infmrp[infmrp$intpos == 1,]$euc))
      mneuc_eucp <- mean(abs(infmrp[infmrp$euc < euc_thr_pos,]$euc))
      mneuc_shrp <- mean(abs(infmrp[infmrp$shrd_rsc_pos == 1,]$euc))
      mneuc_eucshrp <- mean(abs(infmrp[infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1,]$euc))
      mneuc_eucintp <- mean(abs(infmrp[infmrp$euc < euc_thr_pos & infmrp$intpos == 1,]$euc))
      mneuc_shrintp <- mean(abs(infmrp[infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1,]$euc))
      ## Get more info
      mneuc_p <- mean(abs(infmrp$euc))
      mnint_p <- mean(infmrp$int)
      mneuc_nonintp <- mean(abs(infmrp[infmrp$intpos == 0,]$euc))
      mneuc_noneucp <- mean(abs(infmrp[infmrp$euc < euc_thr_pos,]$euc))
      mneuc_nonshrp <- mean(abs(infmrp[infmrp$shrd_rsc_pos == 0,]$euc))
      mneuc_noneucshrp <- mean(abs(infmrp[infmrp$euc > euc_thr_pos & infmrp$shrd_rsc_pos == 0,]$euc))
      mnint_nonintp <- mean(infmrp[infmrp$intpos == 0,]$int)
      mnint_noneucp <- mean(infmrp[infmrp$euc > euc_thr_pos,]$int)
      mnint_nonshrp <- mean(infmrp[infmrp$shrd_rsc_pos == 0,]$int)
      mnint_noneucintp <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$intpos == 0,]$int)
      mneuc_noneucintp <- mean(abs(infmrp[infmrp$euc > euc_thr_pos & infmrp$intpos == 0,]$euc)) 
      ## Get even more info
      mnint_onlyintp_spsort <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$intpos == 1,]$int)
      mnint_onlyeucp_spsort <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$intpos == 0,]$int)
      mnint_onlyshrp_rscasm <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$shrd_rsc_pos == 1,]$int)
      mnint_onlyeucp_rscasm <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 0,]$int)
      mneuc_onlyintp_spsort <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$intpos == 1,]$euc)
      mneuc_onlyeucp_spsort <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$intpos == 0,]$euc) 
      mneuc_onlyshrp_rscasm <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$shrd_rsc_pos == 1,]$euc)
      mneuc_onlyeucp_rscasm <- mean(infmrp[infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 0,]$euc)
      mnint_noneucshrp <- mean(infmrp[infmrp$euc > euc_thr_pos & infmrp$shrd_rsc_pos == 0,]$int)
      
      ###########################################################################################
      ### Get negative matching data
      ###########################################################################################
      infmrn <- infmr[infmr$coocneg == 1,]
      totn <- nrow(infmrn)
      eucn <- (sum(infmrn$euc > euc_thr_neg))/totn
      intn <- (sum(infmrn$intneg == 1))/totn
      shrn <- (sum(infmrn$shrd_rsc_neg == 1))/totn
      eucshrn <- (sum(infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1))/totn
      eucintn <- (sum(infmrn$euc > euc_thr_neg & infmrn$intneg == 1))/totn
      shrintn <- (sum(infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1))/totn
      eucshrintn <- (sum(infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1))/totn
      unxn <- sum(infmrn$euc < euc_thr_neg & infmrn$intneg == 0 & infmrn$shrd_rsc_neg == 0)/totn
      eucn_cor <- eucn - (eucshrn - eucshrintn) - (eucintn - eucshrintn) - eucshrintn
      shrn_cor <- shrn - (eucshrn - eucshrintn) - (shrintn - eucshrintn) - eucshrintn
      intn_cor <- intn - (eucintn - eucshrintn) - (shrintn - eucshrintn) - eucshrintn
      eucshrn_cor <- eucshrn - eucshrintn
      eucintn_cor <- eucintn - eucshrintn
      shrintn_cor <- shrintn - eucshrintn
      print(paste('In', j, 'rep', r, 'matching negatives in nrep',r,'sum to:', c(eucn_cor + shrn_cor + intn_cor + eucshrn_cor + eucintn_cor + shrintn_cor +  eucshrintn + unxn)))
      ## Get corresponding data on mean interaction coefficients  
      mnint_intn <- mean(infmrn[infmrn$intneg == 1,]$int)
      mnint_eucn <- mean(infmrn[infmrn$euc > euc_thr_neg,]$int)
      mnint_shrn <- mean(infmrn[infmrn$shrd_rsc_neg == 1,]$int)
      mnint_eucshrn <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1,]$int)
      mnint_eucintn <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$intneg == 1,]$int)
      mnint_shrintn <- mean(infmrn[infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1,]$int)
      ## Get corresponding data on mean environmental preference similairty  
      mneuc_intn <- mean(abs(infmrn[infmrn$intneg == 1,]$euc))
      mneuc_eucn <- mean(abs(infmrn[infmrn$euc > euc_thr_neg,]$euc))
      mneuc_shrn <- mean(abs(infmrn[infmrn$shrd_rsc_neg == 1,]$euc))
      mneuc_eucshrn <- mean(abs(infmrn[infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1,]$euc))
      mneuc_eucintn <- mean(abs(infmrn[infmrn$euc > euc_thr_neg & infmrn$intneg == 1,]$euc))
      mneuc_shrintn <- mean(abs(infmrn[infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1,]$euc))
      ## Get more info
      mneuc_n <- mean(abs(infmrn$euc))
      mnint_n <- mean(infmrn$int)
      mneuc_nonintn <- mean(abs(infmrn[infmrn$intneg == 0,]$euc))
      mneuc_noneucn <- mean(abs(infmrn[infmrn$euc < euc_thr_neg,]$euc))
      mneuc_nonshrn <- mean(abs(infmrn[infmrn$shrd_rsc_neg == 0,]$euc))
      mneuc_noneucshrn <- mean(abs(infmrn[infmrn$euc < euc_thr_neg & infmrn$shrd_rsc_neg == 0,]$euc))
      mnint_nonintn <- mean(infmrn[infmrn$intneg == 0,]$int)
      mnint_noneucn <- mean(infmrn[infmrn$euc < euc_thr_pos,]$int)
      mnint_nonshrn <- mean(infmrn[infmrn$shrd_rsc_neg == 0,]$int)
      ## Get even more info
      mnint_noneucintn <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$intneg == 0,]$int)
      mneuc_noneucintn <- mean(abs(infmrn[infmrn$euc < euc_thr_neg & infmrn$intneg == 0,]$euc))
      ### Get corresponding info
      mnint_onlyintn_spsort <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$intneg == 1,]$int)
      mnint_onlyeucn_spsort <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$intneg == 0,]$int)
      mnint_onlyshrn_rscasm <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$shrd_rsc_neg == 1,]$int)
      mnint_onlyeucn_rscasm <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 0,]$int)
      mneuc_onlyintn_spsort <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$intneg == 1,]$euc)
      mneuc_onlyeucn_spsort <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$intneg == 0,]$euc)
      mneuc_onlyshrn_rscasm <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$shrd_rsc_neg == 1,]$euc)
      mneuc_onlyeucn_rscasm <- mean(infmrn[infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 0,]$euc)
      mnint_noneucshrn <- mean(infmrn[infmrn$euc < euc_thr_neg & infmrn$shrd_rsc_neg == 0,]$int)
      
      ### Get corresponding interaction data 
      res <- rbind(res,c(totp,eucp,intp,shrp,eucshrp,eucintp,shrintp,eucshrintp,unxp,eucp_cor,shrp_cor,intp_cor,eucshrp_cor,eucintp_cor,shrintp_cor,
                         mneuc_p,mnint_p,mneuc_nonintp,mneuc_noneucp,mneuc_nonshrp,mneuc_noneucshrp,mnint_nonintp,mnint_noneucp,mnint_nonshrp,
                         mnint_intp,mnint_eucp,mnint_shrp,mnint_eucshrp,mnint_eucintp,mnint_shrintp,mneuc_intp,mneuc_eucp,mneuc_shrp,mneuc_eucshrp,mneuc_eucintp,mneuc_shrintp,mnint_noneucintp,mneuc_noneucintp,
                         mnint_onlyintp_spsort,mnint_onlyeucp_spsort,mnint_onlyshrp_rscasm,mnint_onlyeucp_rscasm,mneuc_onlyintp_spsort,mneuc_onlyeucp_spsort,mneuc_onlyshrp_rscasm,mneuc_onlyeucp_rscasm,mnint_noneucshrp,
                         totn,eucn,intn,shrn,eucshrn,eucintn,shrintn,eucshrintn,unxn,eucn_cor,shrn_cor,intn_cor,eucshrn_cor,eucintn_cor,shrintn_cor,
                         mneuc_n,mnint_n,mneuc_nonintn,mneuc_noneucn,mneuc_nonshrn,mneuc_noneucshrn,mnint_nonintn,mnint_noneucn,mnint_nonshrn,                         
                         mnint_intn,mnint_eucn,mnint_shrn,mnint_eucshrn,mnint_eucintn,mnint_shrintn,mneuc_intn,mneuc_eucn,mneuc_shrn,mneuc_eucshrn,mneuc_eucintn,mneuc_shrintn,mnint_noneucintn,mneuc_noneucintn,
                         mnint_onlyintn_spsort,mnint_onlyeucn_spsort,mnint_onlyshrn_rscasm,mnint_onlyeucn_rscasm,mneuc_onlyintn_spsort,mneuc_onlyeucn_spsort,mneuc_onlyshrn_rscasm,mneuc_onlyeucn_rscasm,mnint_noneucshrn))
    }
    resx <- c()
    #k <- 1
    # Get means for all runs -> Ignore Runs with NaN values (f.ex. where there were no negative co-occurrences)
    for (k in 1:ncol(res)) {
      resx <- c(resx, sum(res[!is.na(res[,k]),k])/length(res[!is.na(res[,k]),k]))
    }
    resl <- rbind(resl, c(resx,j))
  }
  resl <- as.data.frame(resl)
  colnames(resl) <- c('totp','eucp','intp','shrp','eucshrp','eucintp','shrintp','eucshrintp','unxp','eucp_cor','shrp_cor','intp_cor','eucshrp_cor','eucintp_cor','shrintp_cor',
                      'mneuc_p','mnint_p','mneuc_nonintp','mneuc_noneucp','mneuc_nonshrp','mneuc_noneucshrp','mnint_nonintp','mnint_noneucp','mnint_nonshrp',
                      'mnint_intp','mnint_eucp','mnint_shrp','mnint_eucshrp','mnint_eucintp','mnint_shrintp','mneuc_intp','mneuc_eucp','mneuc_shrp','mneuc_eucshrp','mneuc_eucintp','mneuc_shrintp','mnint_noneucintp','mneuc_noneucintp',
                      'mnint_onlyintp_spsort','mnint_onlyeucp_spsort','mnint_onlyshrp_rscasm','mnint_onlyeucp_rscasm','mneuc_onlyintp_spsort','mneuc_onlyeucp_spsort','mneuc_onlyshrp_rscasm','mneuc_onlyeucp_rscasm','mnint_noneucshrp',
                      'totn','eucn','intn','shrn','eucshrn','eucintn','shrintn','eucshrintn','unxn','eucn_cor','shrn_cor','intn_cor','eucshrn_cor','eucintn_cor','shrintn_cor',
                      'mneuc_n','mnint_n','mneuc_nonintn','mneuc_noneucn','mneuc_nonshrn','mneuc_noneucshrn','mnint_nonintn','mnint_noneucn','mnint_nonshrn',                         
                      'mnint_intn','mnint_eucn','mnint_shrn','mnint_eucshrn','mnint_eucintn','mnint_shrintn','mneuc_intn','mneuc_eucn','mneuc_shrn','mneuc_eucshrn','mneuc_eucintn','mneuc_shrintn','mnint_noneucintn','mneuc_noneucintn',
                      'mnint_onlyintn_spsort','mnint_onlyeucn_spsort','mnint_onlyshrn_rscasm','mnint_onlyeucn_rscasm','mneuc_onlyintn_spsort','mneuc_onlyeucn_spsort','mneuc_onlyshrn_rscasm','mneuc_onlyeucn_rscasm','mnint_noneucshrn')
  resl <- as.data.frame(resl)
  resl[,1:20] <- sapply(resl[,1:20], as.numeric)
  infm_cooc_comp <- resl
  return(resl)
}

