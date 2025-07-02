###########################################################################################
### SAMSARA - Calculating means of all diversity metrics                                ###
###########################################################################################

### Get means of all results
get_means_of_all_results <- function(div_res_comp,div_rsc_res_comp) { #nst_res_comp
  
  full_res_subdir_comp <- c()
  #j <- unique(div_res_comp$d)[1]
  for (j in unique(div_res_comp$d)) {
    div_res_comp_mn <- apply(div_res_comp[div_res_comp$d == j,colnames(div_res_comp)!='repl'&colnames(div_res_comp)!='d'],2,mean)
    div_rsc_res_comp_mn <- apply(div_rsc_res_comp[div_rsc_res_comp$d == j,colnames(div_rsc_res_comp)!='repl'&colnames(div_rsc_res_comp)!='d'],2,mean)
    full_res_subdir_comp <- rbind(full_res_subdir_comp, c(div_res_comp_mn, div_rsc_res_comp_mn, d=j))
  }
  full_res_subdir_comp <- as.data.frame(full_res_subdir_comp)
  return(full_res_subdir_comp)  
}
