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
    ## Now we assess matches in each individual run 
    #r <- 0 
    res <- c()
    for (r in unique(infm$repl)) {
      infmr <- infmj[infmj$repl == r,]
      
      ###########################################################################################
      ### Get positive matching data
      ###########################################################################################
      infmrp <- infmr[infmr$coocpos == 1,]
      totp <- nrow(infmrp)
      pd <- list(totp = totp,
                 eucp = (sum(infmrp$euc < euc_thr_pos))/totp,
                 intp = (sum(infmrp$intpos == 1))/totp,
                 shrp = (sum(infmrp$shrd_rsc_pos == 1))/totp,
                 eucshrp = (sum(infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1))/totp,
                 eucintp = (sum(infmrp$euc < euc_thr_pos & infmrp$intpos == 1))/totp,
                 shrintp = (sum(infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1))/totp,
                 eucshrintp = (sum(infmrp$euc < euc_thr_pos & infmrp$shrd_rsc_pos == 1 & infmrp$intpos == 1))/totp,
                 unxp = sum(infmrp$euc > euc_thr_pos & infmrp$intpos == 0 & infmrp$shrd_rsc_pos == 0)/totp)
      # All comparisons considered
      pd$eucp_cor <- pd$eucp - (pd$eucshrp - pd$eucshrintp) - (pd$eucintp - pd$eucshrintp) - pd$eucshrintp
      pd$shrp_cor <- pd$shrp - (pd$eucshrp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
      pd$intp_cor <- pd$intp - (pd$eucintp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
      pd$eucshrp_cor <- pd$eucshrp - pd$eucshrintp
      pd$eucintp_cor <- pd$eucintp - pd$eucshrintp
      pd$shrintp_cor <- pd$shrintp - pd$eucshrintp
      # A Resource preference and interactions based matching considered
      pd$eucp_corA <- pd$eucp - pd$eucintp
      pd$intp_corA <- pd$intp - pd$eucintp
      pd$eucintp_corA <- pd$eucintp
      pd$unxpA <- pd$unxp + (pd$shrp - (pd$eucshrp) - (pd$shrintp - pd$eucshrintp))
      # B Resource preference and shared resource association based matching considered
      pd$eucp_corB <- pd$eucp - pd$eucshrp
      pd$shrp_corB <- pd$shrp - pd$eucshrp
      pd$eucshrp_corB <- pd$eucshrp
      pd$unxpB <- pd$unxp + (pd$intp - (pd$eucintp) - (pd$shrintp - pd$eucshrintp))
      resp <- unlist(pd)
      
      ###########################################################################################
      ### Get negative matching data
      ###########################################################################################
      infmrn <- infmr[infmr$coocneg == 1,]
      totn <- nrow(infmrn)
      pd <- list(totn = totn,
                 eucn = (sum(infmrn$euc > euc_thr_neg))/totn,
                 intn = (sum(infmrn$intneg == 1))/totn,
                 shrn = (sum(infmrn$shrd_rsc_neg == 1))/totn,
                 eucshrn = (sum(infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1))/totn,
                 eucintn = (sum(infmrn$euc > euc_thr_neg & infmrn$intneg == 1))/totn,
                 shrintn = (sum(infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1))/totn,
                 eucshrintn = (sum(infmrn$euc > euc_thr_neg & infmrn$shrd_rsc_neg == 1 & infmrn$intneg == 1))/totn,
                 unxn = sum(infmrn$euc < euc_thr_neg & infmrn$intneg == 0 & infmrn$shrd_rsc_neg == 0)/totn)
      # All comparisons considered
      pd$eucn_cor <- pd$eucn - (pd$eucshrn - pd$eucshrintn) - (pd$eucintn - pd$eucshrintn) - pd$eucshrintn
      pd$shrn_cor <- pd$shrn - (pd$eucshrn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
      pd$intn_cor <- pd$intn - (pd$eucintn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
      pd$eucshrn_cor <- pd$eucshrn - pd$eucshrintn
      pd$eucintn_cor <- pd$eucintn - pd$eucshrintn
      pd$shrintn_cor <- pd$shrintn - pd$eucshrintn
      # B Resource preference and shared resource association based matching considered
      pd$eucn_corA <- pd$eucn - pd$eucintn
      pd$intn_corA <- pd$intn - pd$eucintn
      pd$eucintn_corA <- pd$eucintn
      pd$unxnA <- pd$unxn + (pd$shrn - (pd$eucshrn) - (pd$shrintn - pd$eucshrintn))
      #
      pd$eucn_corB <- pd$eucn - pd$eucshrn
      pd$shrn_corB <- pd$shrn - pd$eucshrn
      pd$eucshrn_corB <- pd$eucshrn
      pd$unxnB <- pd$unxn + (pd$intn - (pd$eucintn) - (pd$shrintn - pd$eucshrintn))
      resn <- unlist(pd)
      res <- rbind(res, c(resp,resn))
    }
    resx <- c()
    resy <- c()
    #k <- 1
    # Get means for all runs -> Ignore Runs with NaN values (f.ex. where there were no negative co-occurrences)
    for (k in 1:ncol(res)) {
      # Extract the column
      col_vals <- res[, k]
      # Remove NAs
      valid_vals <- col_vals[!is.na(col_vals)]
      # Compute mean
      mean_val <- mean(valid_vals)
      resx <- c(resx, mean_val)
      # Compute standard deviation
      sd_val <- sd(valid_vals)
      resy <- c(resy, sd_val)
    }    
    names(resx) <- colnames(res)
    names(resy) <- paste0(colnames(res),'_sd')
    resx <- c(resx,resy)
    ###########################################################################################
    ### Correct for number of runs with zero co-occurrences -> Required for calculations later!
    # Number of runs with zeros
    resx['nzerop'] <- sum(res[,'totp'] == 0)
    resx['nzeron'] <- sum(res[,'totn'] == 0)
    # Average number of co-occurrences including zeros
    resx['cor_totp'] <- sum(res[res[,'totp'] > 0,][,'totp'])/sum(res[,'totp'] > 0)
    resx['cor_totn'] <- sum(res[res[,'totn'] > 0,][,'totn'])/sum(res[,'totn'] > 0)
    resl <- rbind(resl, c(resx, d = j))
  }
  resl <- as.data.frame(resl)
  return(resl)
}

