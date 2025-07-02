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


###########################################################################################
### We add noise to our simulated abundances to remove sligh resource preference based correlations 
add_noise_per_habitat <- function(npop, noise_std, noise_mode) {
  if (noise_mode == 'habitat_adj') {
    for (i in r_num) {
      for (j in npop[npop$repl == i,]$Ridx) {
        valvec <- as.numeric(npop[npop$repl == i & npop$Ridx == j,c(1:sp_num)])
        mn <- mean(valvec)
        stdev <- sd(valvec)
        noise_vals <- rnorm(sp_num, mean = 0, sd = stdev*noise_std)
        noise_vals[valvec == 0] <- 0
        npop[npop$repl == i & npop$Ridx == j,c(1:sp_num)] <- valvec + noise_vals 
      }
    }
  } else if (noise_mode == 'species_adj') {
    for (i in r_num) {
      for (j in npop[npop$repl == i,]$Ridx) {
        valvec <- as.numeric(npop[npop$repl == i & npop$Ridx == j,c(1:sp_num)])
        valvec_std <- valvec*noise_std
        #k <- (1:length(valvec))[2] 
        noise_vals <- c()
        for (k in 1:length(valvec)) {
          if (valvec[k] == 0) {
            noise_vals[k] <- 0
          } else {
            noise_vals[k] <- rnorm(1, mean = 0, sd = valvec_std[k]) #valvec[2]
          }
        }
        npop[npop$repl == i & npop$Ridx == j,c(1:sp_num)] <- valvec + noise_vals        
      }
    }
  } else if (noise_mode == 'fixed_val_norm') {
    noise_mat <- matrix(rnorm(nrow(npop)*sp_num, mean = 0, sd = noise_std), ncol = sp_num)
    noise_mat[npop[,c(1:sp_num)] == 0] <- 0
    npop[,c(1:sp_num)] <- npop[,c(1:sp_num)] + noise_mat 
  } else if (noise_mode == 'fixed_val_unif') {
    noise_mat <- matrix(runif(nrow(npop)*sp_num, min = -noise_std, max = noise_std), ncol = sp_num)
    noise_mat[npop[,c(1:sp_num)] == 0] <- 0
    npop[,c(1:sp_num)] <- npop[,c(1:sp_num)] + noise_mat 
  } else {print('Noise input was incorrect')}
  npop[,c(1:sp_num)][npop[,c(1:sp_num)] < 0] <- 0
  return(npop)
}




