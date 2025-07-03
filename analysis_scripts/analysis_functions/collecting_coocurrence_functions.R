#########################################################################################
### SAMSARA - Assessing how well co-occurrences match either drivers                  ###
#########################################################################################

#########################################################################################
### Main analysis function
# Here we assess how well co-occurrences match either interactions or environmental prefrence similarity

###########################################################################################
### Applying all predictions in one function
get_pred_stats <- function(sp_correl, sign_ajd_res, sp_num, r_num, eucdm, intm, sprsc_cor_res, rsc_num, sp_nms, pred_shuf, euc_thr_neg, euc_thr_pos) { 
  infmat_res <- c()
  #r <- 0
  for (r in r_num){
    
    ###########################################################################################
    ### Get all matrices into the same shape (three column file) for further analysis
    # Spearmans correlations adjacency matrix
    spadjmr <- sp_correl[sp_correl$repl == r,1:sp_num]
    spadjmr[is.na(spadjmr)] <- 0
    diag(spadjmr) <- NA
    spadjmr[lower.tri(spadjmr)] <- NA
    rownames(spadjmr) <- colnames(spadjmr) 
    spadjmr <- reshape2::melt(as.matrix(spadjmr))
    # Significant and trimmed co-occurrences matrix
    adjmr <- sign_ajd_res[sign_ajd_res$repl == r,1:sp_num]
    diag(adjmr) <- NA
    adjmr[lower.tri(adjmr)] <- NA
    rownames(adjmr) <- colnames(adjmr) 
    adjmr <- reshape2::melt(as.matrix(adjmr))
    # Interaction matrix
    intmr <- intm[intm$repl == r,1:sp_num]
    diag(intmr) <- NA
    intmr[lower.tri(intmr)] <- NA
    rownames(intmr) <- colnames(intmr) 
    intmr <- reshape2::melt(as.matrix(intmr))
    # Environmental preference euclidean distance matrix 
    eucdmr <- eucdm[eucdm$repl == r,1:sp_num]
    diag(eucdmr) <- NA
    eucdmr[lower.tri(eucdmr)] <- NA
    rownames(eucdmr) <- colnames(eucdmr) 
    eucdmr <- reshape2::melt(as.matrix(eucdmr))
    # Environmental preference euclidean distance matrix -> Trimmed by thresholds calculated in master_slurm_analysis.R  
    eucdmtrimr <- eucdm[eucdm$repl == r,1:sp_num]
    eucdmtrimr[eucdmtrimr > euc_thr_pos & eucdmtrimr < euc_thr_neg] <- 0
    diag(eucdmtrimr) <- NA
    eucdmtrimr[lower.tri(eucdmtrimr)] <- NA
    rownames(eucdmtrimr) <- colnames(eucdmtrimr) 
    eucdmtrimr <- reshape2::melt(as.matrix(eucdmtrimr))
    
    ###########################################################################################
    ### Combine into single dataframe
    infm <- cbind(spadjmr,adjmr$value,intmr$value,eucdmr$value,eucdmtrimr$value)
    colnames(infm) <- c('N1','N2','correl','cooc','int','euc','euctrim')
    infm <- infm[!is.na(infm$cooc),]
    
    ###########################################################################################
    ### Add booleans for matching 
    ###########################################################################################
    infm$coocpos <- ifelse(infm$cooc > 0, 1, 0)
    infm$coocneg <- ifelse(infm$cooc < 0, 1, 0)
    infm$intpos <- ifelse(infm$int < 0, 1, 0)     # Positive coefficients are competition!
    infm$intneg <- ifelse(infm$int > 0, 1, 0)     # Negative coefficients are mutualism!
    infm$euctrimpos <- ifelse(infm$euctrim <= euc_thr_pos & infm$euctrim != 0, 1, 0)
    infm$euctrimneg <- ifelse(infm$euctrim >= euc_thr_neg, 1, 0)

    ###########################################################################################
    ### Assessing environmental triplets
    sprsc_cor_resr <- sprsc_cor_res[sprsc_cor_res$repl == r,1:rsc_num]
    rownames(sprsc_cor_resr) <- sp_nms
    # Function can be found below
    infm <- get_shared_resource_assoc(infm, sprsc_cor_resr)
    
    infmat_res <- rbind(infmat_res, cbind(infm, repl=r))
    print(paste0("Matching has been made: ",r))
  }
  return(infmat_res)
}


#########################################################################################
### Determining environmental triplets (mostly referred to as 'shared resource associations' in code)  
# Get binary on whether two species share an association or inverse association with a resource
get_shared_resource_assoc <- function(infm, sprsc_cor_resr) {
  
  shrd_rsc_assoc_pos <- c()
  shrd_rsc_assoc_neg <- c()
  
  ### Cycle through every species pair
  #k <- 1
  for (k in 1:nrow(infm)) {
    n1 <- as.character(infm[k,]$N1)
    n2 <- as.character(infm[k,]$N2)
    # Assess for either species (1 and 2), whether they have positive or negative associations with resources 
    pos1 <- ifelse(sprsc_cor_resr[rownames(sprsc_cor_resr) == n1,] > 0, 1, 0)
    pos2 <- ifelse(sprsc_cor_resr[rownames(sprsc_cor_resr) == n2,] > 0, 1, 0)
    neg1 <- ifelse(sprsc_cor_resr[rownames(sprsc_cor_resr) == n1,] < 0, 1, 0)
    neg2 <- ifelse(sprsc_cor_resr[rownames(sprsc_cor_resr) == n2,] < 0, 1, 0)
    ### Now we need to apply conditions to determine whether two species share associations with the same resource 
    ## 1st step: Two species have the same association with a resource -> Positive environmental triplet
    # If two species share the same association with a resource - be it positive or negative - they  
    # might be constrained by similar environmental pressures, thereby showing positive co-abundance
    if (max(pos1 + pos2) == 2 | max(neg1 + neg2) == 2) {  # In this condition we check whether there are cases where both species have an association with the same resource(s)
      assoc_shrd_pos <- 1 
    } else {
      assoc_shrd_pos <- 0
    }
    ## 2nd step: Two species have the opposite association with a resource -> Negative environmental triplet
    # When one species associates positively and the other negatively with a resources abundance
    # their might show opposing environmental constraints and be negatively linked in their abundances  
    if (max(pos1 + neg2) == 2 | max(pos2 + neg1) == 2) {  # Now we check whether one species might have a positive association with the one resource while the other has a negative association with the same resource
      assoc_shrd_neg <- 1 
    } else {
      assoc_shrd_neg <- 0
    }
    ## 3rd step: In the very rare cases where two species both have a positive and a negative environmental triplet we determine them not associated (no triplets) 
    # This has never occurred in testing, but is a precaution we have implemented nonetheless
    if (assoc_shrd_neg == 1 & assoc_shrd_pos == 1) {
      assoc_shrd_pos <- 0
      assoc_shrd_neg <- 0
    }
    shrd_rsc_assoc_pos <- rbind(shrd_rsc_assoc_pos, c(n1,n2,assoc_shrd_pos))
    shrd_rsc_assoc_neg <- rbind(shrd_rsc_assoc_neg, c(n1,n2,assoc_shrd_neg))      
  }
  infm$shrd_rsc_pos <-  as.numeric(shrd_rsc_assoc_pos[,3])
  infm$shrd_rsc_neg <-  as.numeric(shrd_rsc_assoc_neg[,3])
  return(infm)
}

