###########################################################################################
### SAMSARA - Data preparation for subsampling habitats within volumes                  ###
###########################################################################################

#nxyz <- nxyzr
#npop <- npopr
#kabs <- kabsr
#i <- 1
#j <- 1
#l <- 1
#sdim <- sdim_vec[1]
#sdim <- 0.4 ### test

###########################################################################################
### Subsample function:
subsampl <- function(sdim, nxyz, npop, kabs, sp_nms, rsc_nms) {
  sdimlen <- sdim
  # Define possible squares to sample: 
  cubs <- c()
  for (i in 1:length(seq(0,1-sdim,sdim))) {
    for (j in 1:length(seq(0,1-sdim,sdim))) {
      for (l in 1:length(seq(0,1-sdim,sdim))) {
        cubs <- rbind(cubs, cbind(seq(0,1-sdim,sdim)[i], seq(0,1-sdim,sdim)[j], seq(0,1-sdim,sdim)[l]))
      }  
    }
  }
  cubs <- as.data.frame(cubs)
  smpl_cubs <- as.data.frame(cubs)     # Cubes extents (:+sdimlen) to be samples
  smpl_cubs$smpl_cubs_idx <- rownames(smpl_cubs)
  #smpl_idx <- sample(nrow(cubs))[1:nsmpl]
  #smpl_cubs <- as.data.frame(cubs[smpl_idx,])     # Cubes extents (:+sdimlen) to be samples 
  #smpl_cubs$smpl_cubs_idx <- seq(1,nsmpl,1)
  # Collect all samples within sample cubes 
  npop_cubs <- c()
  kabs_cubs <- c()
  #i <- 27
  #j <- 142
  for (i in 1:nrow(smpl_cubs)) {
    for (j in 1:nrow(nxyz)) {
      if (nxyz[j,1] > smpl_cubs[i,1] && nxyz[j,1] < smpl_cubs[i,1]+sdimlen && 
          nxyz[j,2] > smpl_cubs[i,2] && nxyz[j,2] < smpl_cubs[i,2]+sdimlen &&
          nxyz[j,3] > smpl_cubs[i,3] && nxyz[j,3] < smpl_cubs[i,3]+sdimlen) {
        npop_cubs <- rbind(npop_cubs, cbind(npop[j,],i, j)) #npop[j,]
        kabs_cubs <- rbind(kabs_cubs, cbind(kabs[j,],i, j)) #npop[j,]
      }
    }
  }
  if (is.null(npop_cubs)) {
    out <- "no_cubes_with_habs"
  } else if (nrow(npop_cubs) <= 1) {
    out <- "only_one_cube_with_one_hab"
  } else {
    colnames(npop_cubs) <- c(colnames(npop), 'smpl_cubs_idx','hab_idx')
    colnames(kabs_cubs) <- c(colnames(kabs), 'smpl_cubs_idx','hab_idx')
    # Aggregate samples to composite measurements  
    npop_comp <- c()
    kabs_comp <- c()
    for (i in unique(npop_cubs$smpl_cubs_idx)) {
      nhab <- sum(npop_cubs$smpl_cubs_idx == i)
      pdat <- colSums(npop_cubs[npop_cubs$smpl_cubs_idx == i, colnames(npop_cubs) %in% sp_nms,])
      npop_comp <- rbind(npop_comp, c(pdat, i, nhab, unique(npop_cubs$repl)))
      kdat <- colSums(kabs_cubs[kabs_cubs$smpl_cubs_idx == i, colnames(kabs_cubs) %in% rsc_nms,])
      kabs_comp <- rbind(kabs_comp, c(kdat, i, nhab, unique(kabs_cubs$repl)))
    }
    if (nrow(kabs_comp) == 1) {
      out <- "only_one_cubes_with_multiple_habs"
    } else {  
      colnames(npop_comp)[length(sp_nms)+1:3] <- c('smpl_cub_idx','comp_hab_num','repl')
      colnames(kabs_comp)[length(rsc_nms)+1:3] <- c('smpl_cub_idx','comp_hab_num','repl')
      npop_compra <- npop_comp
      npop_compra[,1:length(sp_nms)] <- npop_comp[,1:length(sp_nms)]/rowSums(npop_comp[,1:length(sp_nms)]) 
      kabs_compra <- kabs_comp
      kabs_compra[,1:length(rsc_nms)] <- kabs_comp[,1:length(rsc_nms)]/rowSums(kabs_comp[,1:length(rsc_nms)])
      # Delete those coordinates which are not captured in the data
      npop_comp <- as.data.frame(npop_comp)
      #smpl_cubs <- smpl_cubs[smpl_cubs$smpl_cubs_idx %in% unique(npop_comp$smpl_cub_idx),]
      out <- list(smpl_cubs, npop_cubs, npop_comp, kabs_comp, npop_compra, kabs_compra)
    }
  }
  #out
  return(out)
}

###########################################################################################
### Applies subsample function at multiple scales of the cube:
# sdim_vec = A vector of all different scales (side length) samples are groupes (comp = composite) 
#r <- 2
#d <- sdim_vec[1]
get_scale_coms <- function(sdim_vec, nxyz, npop, sp_nms, rsc_nms, rsc_num, r_num) { 
  coords_comp <- c()
  npop_comp_log <- c()
  kabs_comp_log <- c()
  npop_compra_log <- c()
  kabs_compra_log <- c()
  for (d in sdim_vec) {
    for (r in r_num) {
      nxyzr <- nxyz[nxyz$repl == r,]
      npopr <- npop[npop$repl == r,] 
      kabsr <- kabs[kabs$repl == r,]
      
      out <- subsampl(d, nxyzr, npopr, kabsr, sp_nms, rsc_nms)
      if (length(out) == 1) {
        print(paste("Cube sampling",d,r,out))
        next
      } else {
        coords_comp <- rbind(coords_comp, cbind(as.data.frame(out[1]), r, d))
        npop_comp_log <- rbind(npop_comp_log, cbind(as.data.frame(out[3]), d))
        kabs_comp_log <- rbind(kabs_comp_log, cbind(as.data.frame(out[4]), d))
        npop_compra_log <- rbind(npop_compra_log, cbind(as.data.frame(out[5]), d))
        kabs_compra_log <- rbind(kabs_compra_log, cbind(as.data.frame(out[6]), d))
        print(paste("Cube sampling",d,r))
      }
    }
  }
  colnames(coords_comp)[colnames(coords_comp) == 'r'] <- c('repl') 
  scale_com <- list(npop_comp_log,kabs_comp_log,coords_comp,npop_compra_log,kabs_compra_log)
  return(scale_com)
}
