#########################################################################################
### SAMSARA - Miscelaneous analysis functions                                         ###
#########################################################################################

#########################################################################################
### Get quartiles threshold values for environmental preference similarity
get_Qs <- function(eucdm, r_num) {
  mean_Q1 <- c()
  mean_Q3 <- c()
  #r <- 0
  for (r in r_num) {
    eucdmr <- eucdm[eucdm$repl == r,] 
    mean_Q1 <- c(mean_Q1, quantile(as.numeric(eucdmr[lower.tri(eucdmr)]), c(.25)))
    mean_Q3 <- c(mean_Q3, quantile(as.numeric(eucdmr[lower.tri(eucdmr)]), c(.75)))
  }
  bc_thr_neg <- mean(mean_Q1)           # Lower threshold for testing predictions
  bc_thr_pos <- mean(mean_Q3)           # Upper threshold for testing predictions
  res <- c(bc_thr_neg,bc_thr_pos)
  return(res)
}
