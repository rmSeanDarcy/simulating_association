###########################################################################################
### SAMSARA - Data preparation for subsampling habitats within volumes                  ###
###########################################################################################

#nxyz <- nxyzr
#npop <- npopr
#kabs <- kabsr
#i <- 1
#j <- 1
#l <- 1
#sdim <- sdim_vec[8]

###########################################################################################
### Subsample function:
subsampl_controlled <- function(sdim, nxyz, npop, kabs, sp_nms, rsc_nms) {
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
  colnames(npop_cubs) <- c(colnames(npop), 'smpl_cubs_idx','hab_idx')
  colnames(kabs_cubs) <- c(colnames(kabs), 'smpl_cubs_idx','hab_idx')
  npop_cubs$dim <- sdim
  kabs_cubs$dim <- sdim
  out <- list(smpl_cubs,npop_cubs,kabs_cubs)
  
  return(out)
}

###########################################################################################
### Applies subsample function at multiple scales of the cube:
# sdim_vec = A vector of all different scales (side length) samples are groupes (comp = composite) 
#r <- 2
#d <- sdim_vec[8]
get_scale_coms_controlled <- function(sdim_vec, nxyz, npop, sp_nms, rsc_nms, rsc_num, r_num) { 
  coords_comp <- c()
  npop_comp_log <- c()
  kabs_comp_log <- c()
  #d <- sdim_vec[1]
  for (d in sdim_vec) {
    #r <- r_num
    for (r in r_num) {
      nxyzr <- nxyz[nxyz$repl == r,]
      npopr <- npop[npop$repl == r,] 
      kabsr <- kabs[kabs$repl == r,]
      
      out <- subsampl_controlled(d, nxyzr, npopr, kabsr, sp_nms, rsc_nms)
      
      coords_comp <- rbind(coords_comp, cbind(as.data.frame(out[1]), r, d))
      npop_comp_log <- rbind(npop_comp_log, cbind(as.data.frame(out[2])))
      kabs_comp_log <- rbind(kabs_comp_log, cbind(as.data.frame(out[3])))
      print(paste("Cube sampling",d,r))
    }
  }
  scale_com <- list(coords_comp,npop_comp_log,kabs_comp_log)
  return(scale_com)
}



###########################################################################################
get_equal_cubsums <- function(npop_cub,kabs_cub,sp_num,rsc_num,coord_cub) {
  npop_cubd_sums <- c()
  kabs_cubd_sums <- c()
  coord_cubd_sums <- c()
  
  #d <- sdim_vec[1]
  for (d in sdim_vec) {
    npop_cubd <- npop_cub[npop_cub$dim == d,]
    kabs_cubd <- kabs_cub[kabs_cub$dim == d,]
    coord_cubd <- coord_cub[coord_cub$d == d,]
    
    npop_cubr_sums <- c()
    kabs_cubr_sums <- c()
    coord_cubr_sums <- c()
    #r <- 0
    for (r in unique(npop_cubd$repl)) {
      npop_cubr <- npop_cubd[npop_cubd$repl == r,]
      kabs_cubr <- kabs_cubd[kabs_cubd$repl == r,]
      coord_cubr <- coord_cubd[coord_cubd$r == r,]
      ###########################################################################################
      ##### I will use some logic here to be able to sample lower volumes too
      ### In stead of randomly sampling 25 cubes and then checking for the lowest habitat number as a maximum test
      ### I take the 25 cubes with the highest amount of habitats. This does not introduce a bias as there is no
      ### determinism to where habitats are located. In effect it should be just as random as were I to select 25 at random.
      #sel_cub <- sample(unique(npop_cubr$smpl_cubs_idx),25)
      # Check minimum value for max 25 samples
      tab_idx <- table(npop_cubr$smpl_cubs_idx)
      nms_idx <- names(tab_idx)[rev(order(tab_idx))][1:25]
      npop_cubr <- npop_cubr[npop_cubr$smpl_cubs_idx %in% nms_idx,]
      
      hts <- min(table(npop_cubr$smpl_cubs_idx))
      
      npop_cubi_sums <- c()
      kabs_cubi_sums <- c()
      coord_cubi_sums <- c()
      #i <- unique(npop_cubr$smpl_cubs_idx)[1]
      for (i in unique(npop_cubr$smpl_cubs_idx)) {
        npop_cubi <- npop_cubr[npop_cubr$smpl_cubs_idx == i,]
        
        smplidx <- sample(nrow(npop_cubi),hts)
        
        npop_cubi <- npop_cubi[smplidx,]
        npop_cubi_sums <- rbind(npop_cubi_sums, c(colSums(npop_cubi[,1:sp_num]), d, r, hts, i))
        
        kabs_cubi <- kabs_cubr[kabs_cubr$smpl_cubs_idx == i,]
        kabs_cubi <- kabs_cubi[smplidx,]
        kabs_cubi_sums <- rbind(kabs_cubi_sums, c(colSums(kabs_cubi[,1:rsc_num]), d, r, hts, i))
        
        coord_cubi <- coord_cubr[coord_cubr$smpl_cubs_idx == i,]
        coord_cubi_sums <- rbind(coord_cubi_sums, coord_cubi)
        
        print(paste('Done for dim:', d, 'for rep:', r, 'for cubeindex:', i, 'with n_hab:', hts))
      }
      npop_cubr_sums <- rbind(npop_cubr_sums, npop_cubi_sums)
      kabs_cubr_sums <- rbind(kabs_cubr_sums, kabs_cubi_sums)
      coord_cubr_sums <- rbind(coord_cubr_sums, coord_cubi_sums)
    }
    npop_cubd_sums <- rbind(npop_cubr_sums, npop_cubd_sums)
    kabs_cubd_sums <- rbind(kabs_cubr_sums, kabs_cubd_sums)
    coord_cubd_sums <- rbind(coord_cubr_sums, coord_cubd_sums)
  }
  npop_cubd_sums <- as.data.frame(npop_cubd_sums)
  kabs_cubd_sums <- as.data.frame(kabs_cubd_sums)
  coord_cubd_sums <- as.data.frame(coord_cubd_sums)
  colnames(npop_cubd_sums) <- c(colnames(npop_cub[,1:sp_num]), 'd', 'repl', 'hab', 'smpl_cub_idx')
  colnames(kabs_cubd_sums) <- c(colnames(kabs_cub[,1:rsc_num]), 'd', 'repl', 'hab', 'smpl_cub_idx')
  colnames(coord_cubd_sums) <- c('V1','V2','V3','smpl_cubs_idx','repl','d')
  return(list(npop_cubd_sums,kabs_cubd_sums,coord_cubd_sums))    
}

###########################################################################################
get_equal_cubsums_fullyequal <- function(npop_cub,kabs_cub,sp_num,rsc_num,coord_cub) {
  npop_cubd_sums <- c()
  kabs_cubd_sums <- c()
  coord_cubd_sums <- c()
  
  #d <- sdim_vec[1]
  for (d in sdim_vec) {
    npop_cubd <- npop_cub[npop_cub$dim == d,]
    kabs_cubd <- kabs_cub[kabs_cub$dim == d,]
    coord_cubd <- coord_cub[coord_cub$d == d,]
    
    npop_cubr_sums <- c()
    kabs_cubr_sums <- c()
    coord_cubr_sums <- c()
    
    hts_all <- c()
    #r <- 0
    for (r in unique(npop_cubd$repl)) {
      npop_cubr <- npop_cubd[npop_cubd$repl == r,]
      ###########################################################################################
      ### Same as above: Now I check the minimum value for max 25 samples from all runs and use this value for all runs
      tab_idx <- table(npop_cubr$smpl_cubs_idx)
      nms_idx <- names(tab_idx)[rev(order(tab_idx))][1:25]
      npop_cubr <- npop_cubr[npop_cubr$smpl_cubs_idx %in% nms_idx,]
      
      hts_all <- c(hts_all,min(table(npop_cubr$smpl_cubs_idx)))
    }
    hts <- min(hts_all)
    #r <- 0
    for (r in unique(npop_cubd$repl)) {
      npop_cubr <- npop_cubd[npop_cubd$repl == r,]
      tab_idx <- table(npop_cubr$smpl_cubs_idx)
      nms_idx <- names(tab_idx)[rev(order(tab_idx))][1:25]
      npop_cubr <- npop_cubr[npop_cubr$smpl_cubs_idx %in% nms_idx,]
      
      kabs_cubr <- kabs_cubd[kabs_cubd$repl == r,]
      coord_cubr <- coord_cubd[coord_cubd$r == r,]
      
      npop_cubi_sums <- c()
      kabs_cubi_sums <- c()
      coord_cubi_sums <- c()
      #i <- unique(npop_cubr$smpl_cubs_idx)[1]
      for (i in unique(npop_cubr$smpl_cubs_idx)) {
        # Take one cube and list all habitats within
        npop_cubi <- npop_cubr[npop_cubr$smpl_cubs_idx == i,]
        # Select (hts) number of habitats (indices) within
        smplidx <- sample(nrow(npop_cubi),hts)
        npop_cubi <- npop_cubi[smplidx,]
        # Aggregate pop sums
        npop_cubi_sums <- rbind(npop_cubi_sums, c(colSums(npop_cubi[,1:sp_num]), d, r, hts, i))
        # Same for resources
        kabs_cubi <- kabs_cubr[kabs_cubr$smpl_cubs_idx == i,]
        kabs_cubi <- kabs_cubi[smplidx,]
        kabs_cubi_sums <- rbind(kabs_cubi_sums, c(colSums(kabs_cubi[,1:rsc_num]), d, r, hts, i))
        # Collect the cubes coordinates        
        coord_cubi <- coord_cubr[coord_cubr$smpl_cubs_idx == i,]
        coord_cubi_sums <- rbind(coord_cubi_sums, coord_cubi)
        
        print(paste('Done for dim:', d, 'for rep:', r, 'for cubeindex:', i, 'with n_hab:', hts))
      }
      npop_cubr_sums <- rbind(npop_cubr_sums, npop_cubi_sums)
      kabs_cubr_sums <- rbind(kabs_cubr_sums, kabs_cubi_sums)
      coord_cubr_sums <- rbind(coord_cubr_sums, coord_cubi_sums)
    }
    npop_cubd_sums <- rbind(npop_cubr_sums, npop_cubd_sums)
    kabs_cubd_sums <- rbind(kabs_cubr_sums, kabs_cubd_sums)
    coord_cubd_sums <- rbind(coord_cubr_sums, coord_cubd_sums)
  }
  npop_cubd_sums <- as.data.frame(npop_cubd_sums)
  kabs_cubd_sums <- as.data.frame(kabs_cubd_sums)
  coord_cubd_sums <- as.data.frame(coord_cubd_sums)
  colnames(npop_cubd_sums) <- c(colnames(npop_cub[,1:sp_num]), 'd', 'repl', 'hab', 'smpl_cub_idx')
  colnames(kabs_cubd_sums) <- c(colnames(kabs_cub[,1:rsc_num]), 'd', 'repl', 'hab', 'smpl_cub_idx')
  colnames(coord_cubd_sums) <- c('V1','V2','V3','smpl_cubs_idx','repl','d')
  return(list(npop_cubd_sums,kabs_cubd_sums,coord_cubd_sums))    
}


