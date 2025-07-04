###########################################################################################
### SAMSARA - Calculate significant associations (sp x sp and sp x rsc)                 ###
###########################################################################################

###########################################################################################
### Calculates all species x species Spearmans correlations (without significance testing) 
get_signif_adj_all <- function(npop, sp_num, r_num) {
  sign_ajd_res <- c()
  #r <- 0
  for (r in r_num){
    npopr <- npop[npop$repl == r,] 
    npop_cor <- cor(npopr[,1:sp_num], method = 'spearman') #, use = "complete.obs"
    sign_ajd_res <- rbind(sign_ajd_res, cbind(npop_cor,repl = r))
    print(paste("Exporting species correlations: ", r))
  }
  sign_ajd_res <- as.data.frame(sign_ajd_res)
  return(sign_ajd_res)
}

###########################################################################################
### Calculates significant species x species associations with significance testing
get_signif_adj <- function(npop, p_val_thr2, nperm2, sp_num, r_num) {
  sign_ajd_res <- c()
  
  # Cycles through all runs without errors (r_num)
  #r <- 0
  for (r in r_num){
    npopr <- npop[npop$repl == r,] 
    # Calculates Spearmans correlation coefficients for all pairs (ignore 'the standard deviation is zero' error message)
    npop_cor <- cor(npopr[,1:sp_num], method = 'spearman') 
    npop_cor[is.na(npop_cor)] <- 0
    diag(npop_cor) <- 0
    # Divide into positive and negative correlations
    npop_cor_pos <- ifelse(npop_cor > 0, npop_cor,0)
    npop_cor_neg <- ifelse(npop_cor < 0, npop_cor,0)
    # Select nullmodel (from vegan package) -> 'r0_samp' preserves row sums and shuffles values within row
    nm <- nullmodel(npopr[,1:sp_num], "r0_samp")  
    # Create empty matrices to fill with results from each permutation
    signifm <- matrix(0,ncol(npop_cor),ncol(npop_cor))
    signifm_pos <- matrix(0,ncol(npop_cor),ncol(npop_cor))
    signifm_neg <- matrix(0,ncol(npop_cor),ncol(npop_cor))
    for (i in 1:nperm2) { 
      sign <- matrix(0,ncol(npop_cor),ncol(npop_cor))
      sign_pos <- matrix(0,ncol(npop_cor),ncol(npop_cor))
      sign_neg <- matrix(0,ncol(npop_cor),ncol(npop_cor))
      # Generate a randomly shuffled abundance table using the nm set above
      rdm_com <- as.data.frame(simulate(nm)) 
      # Calculate Spearmans correlation coefficients from this new random table
      rdm_cor <- cor(rdm_com, method = 'spearman')
      # Divide into adjacency matrices with positive and negative values 
      rdm_cor_pos <- ifelse(rdm_cor > 0, rdm_cor,0)
      rdm_cor_neg <- ifelse(rdm_cor < 0, rdm_cor,0)
      # Now we check whether the actual coefficients were larger than the randomly shuffled ones (smaller for negative)  
      sign_pos[npop_cor_pos > rdm_cor_pos] <- 1 
      sign_neg[npop_cor_neg < rdm_cor_neg] <- 1
      # Then we count whether this was the case for every permutation
      signifm_pos <- signifm_pos + sign_pos
      signifm_neg <- signifm_neg + sign_neg
      # A robust correlation should be significantly better than correlations that emerge from chance
    }
    # Combine both positive and negative permuted tests into a single matrix and divide by the number of permutations
    # -> Gives the percentage of tests where that were not 'robust' (1-)
    sign_ajd <- 1-((signifm_pos + signifm_neg)/nperm2)
    # We adjust our p-Value (initially set to 0.05) for multiple testing via Bonferroni correction
    # The number of associations we test for 10 species is 45 (half of adjacency matrix, without diagonal) 
    n_tests <-  ((ncol(npop_cor)**2)-ncol(npop_cor))/2
    # Here we trim the actual adjacency matrix to only those values that were tested significant  
    sign_ajd <- ifelse(sign_ajd < p_val_thr2/n_tests, npop_cor, 0)  # p_val_thr2/n_tests -> Corrected p-Value
    sign_ajd <- cbind(sign_ajd, r)
    colnames(sign_ajd) <- c(colnames(npopr[,1:sp_num]), 'repl')
    sign_ajd_res <- rbind(sign_ajd_res, sign_ajd)
    print(paste("Generating null models for significant species co-occurrences : ", r))
  }
  sign_ajd_res <- as.data.frame(sign_ajd_res)
  return(sign_ajd_res)
}


###########################################################################################
### Calculates significant species x resource associations across habitats
# Follows the same principle as above (species x species)
get_signif_sprsc <- function(npop, kabs, p_val_thr2, nperm2, sp_num, rsc_num, r_num) {
  sprsc_cor_res <- c()
  
  for (r in r_num){
    npopr <- npop[npop$repl == r,] 
    kabsr <- kabs[kabs$repl == r,]
    sprsc_cor <- cor(npopr[,1:sp_num], kabsr[,1:rsc_num], method = 'spearman') #, use = "complete.obs"
    sprsc_cor[is.na(sprsc_cor)] <- 0
    
    sp_nm <- nullmodel(npopr[,1:sp_num], "r0_samp")
    rsc_nm <- nullmodel(kabsr[,1:rsc_num], "r0_samp")
    sprsc_signifm_pos <- matrix(0,sp_num,rsc_num)
    sprsc_signifm_neg <- matrix(0,sp_num,rsc_num)
    
    for (i in 1:nperm2) {
      sign_pos <- matrix(0,sp_num,rsc_num)
      sign_neg <- matrix(0,sp_num,rsc_num)
      rdm_sp <- as.data.frame(simulate(sp_nm)) 
      rdm_rsc <- as.data.frame(simulate(rsc_nm)) 
      rdm_sprsc_cor <- cor(rdm_sp, rdm_rsc, method = 'spearman')
      rdm_sprsc_cor_pos <- ifelse(rdm_sprsc_cor > 0, rdm_sprsc_cor, 0)
      rdm_sprsc_cor_neg <- ifelse(rdm_sprsc_cor < 0, rdm_sprsc_cor, 0)
      sign_pos[sprsc_cor > rdm_sprsc_cor_pos] <- 1
      sign_neg[sprsc_cor < rdm_sprsc_cor_neg] <- 1
      sprsc_signifm_pos <- sprsc_signifm_pos + sign_pos 
      sprsc_signifm_neg <- sprsc_signifm_neg + sign_neg 
    }
    sign_sprsc <- 1-((sprsc_signifm_pos+sprsc_signifm_neg)/nperm2)    
    n_tests <-  ncol(sign_sprsc)*nrow(sign_sprsc)
    sign_sprsc <- ifelse(sign_sprsc < p_val_thr2/n_tests, sprsc_cor, 0) 
    sign_sprsc <- cbind(sign_sprsc, r)
    colnames(sign_sprsc) <- c(colnames(kabsr[,1:rsc_num]), 'repl')
    sprsc_cor_res <- rbind(sprsc_cor_res, sign_sprsc)
    print(paste("Generating null models for significant species resource associations : ", r))
  }
  sprsc_cor_res <- as.data.frame(sprsc_cor_res)
  return(sprsc_cor_res)
}

