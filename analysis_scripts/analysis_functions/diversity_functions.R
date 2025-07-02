###########################################################################################
### SAMSARA - Result analysis functions                                                 ###
###########################################################################################

###########################################################################################
### 1: Diversity metrics for species                                                    ###
div_mtrcs <- function(npop, sp_num, hab_num, r_num, eucdm) {
  com_metrics <- c()
  
  # Cycles through every run (without errors -> r_num)
  for (r in r_num) {
    npopr <- npop[npop$repl == r,]
    npopr <- npopr[,1:sp_num]
    adiv <- c()
    # Calculates richness, Simpsons and Shannon diversity for each habitat individually
    for (i in 1:hab_num){
      adivm <- c(sum(npopr[i,] > 0), vegan::diversity(npopr[i,], index = "simpson"), vegan::diversity(npopr[i,], index = "shannon"))
      adiv <- rbind(adiv,adivm)
    }
    # Calculates Bray Curtis dissimilarity with all habitats
    bdiv <- as.vector(vegan::vegdist(npopr, method = 'bray'))  
    bdiv <- bdiv[bdiv<1 & !is.na(bdiv)]
    # Mean environmental preference similarity 
    mn_rscprf <- mean(as.matrix(eucdm[,1:sp_num]))
    com_metrics <- rbind(com_metrics, c(r, mean(adiv[,1]), sd(adiv[,1]), mean(adiv[,2]), sd(adiv[,2]),
                                        mean(adiv[,3]), sd(adiv[,3]), mean(bdiv), sd(bdiv), mn_rscprf))
    print(paste("Diversity metrics: r = ",c(r),sep = ''))
  }
  com_metrics <- as.data.frame(com_metrics)
  colnames(com_metrics) <- c("repl","mn_rch","sd_rch","mn_simp","sd_simp","mn_shan","sd_shan","mn_bc","sd_bc","mn_rscprf")
  return(com_metrics)
}

###########################################################################################
### 2: Diversity metrics for resources -> Same exact calculations as done for species above
div_mtrcs_rsc <- function(kabs, rsc_num, hab_num, r_num) {
  com_metrics <- c()
  for (r in r_num) {
    kabsr <- kabs[kabs$repl == r,]
    kabsr <- kabsr[,1:rsc_num]
    adiv <- c()
    for (i in 1:hab_num){
      adivm <- c(sum(kabsr[i,] > 0), vegan::diversity(kabsr[i,], index = "simpson"), vegan::diversity(kabsr[i,], index = "shannon"))
      adiv <- rbind(adiv,adivm)
    }
    bdiv <- as.vector(vegan::vegdist(kabsr, method = 'bray'))  
    bdiv <- bdiv[bdiv<1 & !is.na(bdiv)]
    com_metrics <- rbind(com_metrics, c(r, mean(adiv[,1]), sd(adiv[,1]), mean(adiv[,2]), sd(adiv[,2]),
                                        mean(adiv[,3]), sd(adiv[,3]), mean(bdiv), sd(bdiv)))
    print(paste("Diversity resource metrics: r = ",c(r),sep = ''))
  }
  com_metrics <- as.data.frame(com_metrics)
  colnames(com_metrics) <- c("repl","mn_rch_rsc","sd_rch_rsc","mn_simp_rsc","sd_simp_rsc","mn_shan_rsc","sd_shan_rsc","mn_bc_rsc","sd_bc_rsc")
  return(com_metrics)
}



