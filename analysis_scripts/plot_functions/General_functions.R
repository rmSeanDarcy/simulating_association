###########################################################################################
### SAMSARA - Loading Set of analyses data and preparing for plotting                   ###
###########################################################################################

###########################################################################################
### ggplot2 plotting modifications
cxy <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
             axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())
cxyl <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")
cx <- theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())
cxl <- theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")
cy <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())
cyl <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")
cl <- theme(legend.position = "none")

###########################################################################################
### Get info about quartiles in BC resource preference data
get_Qs <- function(bcdm, r_num) {
  mean_Q1 <- c()
  mean_Q3 <- c()
  for (r in r_num) {
    bcdmr <- bcdm[bcdm$repl == r,] 
    mean_Q1 <- c(mean_Q1, quantile(as.numeric(bcdmr[lower.tri(bcdmr)]), c(.25)))
    mean_Q3 <- c(mean_Q3, quantile(as.numeric(bcdmr[lower.tri(bcdmr)]), c(.75)))
  }
  bc_thr_neg <- mean(mean_Q1)           # Lower threshold for testing predictions
  bc_thr_pos <- mean(mean_Q3)           # Upper threshold for testing predictions
  res <- c(bc_thr_neg,bc_thr_pos)
  return(res)
}

###########################################################################################
##### Loading data for plotting
#i <- 7
get_runs <- function(parent_set_of_analyses, runs_incl) {
  resl <- load_data_fun(parent_set_of_analyses, set_of_analyses = runs_incl[1])
  pd <- resl[[1]]
  fd <- resl[[2]]
  sd <- resl[[3]]
  if (length(runs_incl) >= 2) {
    for (i in 2:length(runs_incl)) {
      
      resl <- load_data_fun(parent_set_of_analyses, set_of_analyses = runs_incl[i])
      pd <- rbind(pd,resl[[1]])
      fd <- rbind(fd,resl[[2]])
      sd <- rbind(sd,resl[[3]])
      print(runs_incl[i])
    }
  }
  out <- list(pd,fd,sd)
  return(out)
}
###########################################################################################
load_data_fun <- function(parent_set_of_analyses, set_of_analyses) {
  parentdirs <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,"/",set_of_analyses,sep=''), full.names=TRUE) #, pattern = "Z"
  parentdirnms <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,"/",set_of_analyses,sep='')) #, pattern = "Z"
  pd <- c()
  fd <- c()
  sd <- c()
  #i <- 4
  for (i in 1:length(parentdirnms)) {
    if (file.exists(paste0(parentdirs[i],'/infm_res.csv')) & file.exists(paste0(parentdirs[i],'/full_res.csv')) & file.exists(paste0(parentdirs[i],'/pred_res.csv'))) {
      print(parentdirnms[i])
      pred_res <- read.csv(paste0(parentdirs[i],'/infm_res.csv'))[-1]
      pred_res$run_nm <- parentdirnms[i]
      full_res <- read.csv(paste0(parentdirs[i],'/full_res.csv'))[-1]
      full_res$run_nm <- parentdirnms[i]
      scor_res <- read.csv(paste0(parentdirs[i],'/pred_res.csv'))[-1]
      scor_res$run_nm <- parentdirnms[i]
      pd <- rbind(pd,pred_res)
      fd <- rbind(fd,full_res)
      sd <- rbind(sd,scor_res)
    }
    else {
      print(paste0('Required files do not exist in: ',parentdirnms[i]))
    }
  }
  resl <- list(pd,fd,sd)
  return(resl)
}
###########################################################################################
##### Prep functions for plotting
prep_plot_data <- function(fd,pd,sd) {
  pd$run_nm <- as.factor(pd$run_nm)
  pd$d <- as.factor(pd$d)
  hab_levels <- levels(pd$d)[grep('HabSubsmpl', levels(pd$d))]
  cub_levels <- levels(pd$d)[grep('HabSubsmpl', levels(pd$d))]
  hab_levels_info <- as.data.frame(cbind(hab_levels, unlist(strsplit(hab_levels, '_'))[2*(1:length(hab_levels))]))
  colnames(hab_levels_info) <- c('d','hab_num')
  cub_levels_info <- as.data.frame(cbind(cub_levels, unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))],
                                         unlist(strsplit(unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))], '_'))[2*(1:length(cub_levels))],
                                         unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))-1],
                                         unlist(strsplit(unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))-1], '_'))[2*(1:length(cub_levels))]))
  colnames(cub_levels_info) <- c('d','dim_lab','dim','hab_num_lab','hab_num')
  ### Create data subset for community metrics
  fd_hab <- fd[fd$d %in% hab_levels,]
  fd_hab <- left_join(hab_levels_info,fd_hab, by = 'd', multiple = 'all')
  fd_hab$hab_num <- factor(as.numeric(fd_hab$hab_num), levels = sort(as.numeric(unique(fd_hab$hab_num))))
  fd_cub <- fd[fd$d %in% cub_levels,]
  fd_cub <- left_join(cub_levels_info, fd_cub, by = 'd', multiple = 'all')
  fd_cub$hab_num <- factor(as.numeric(fd_cub$hab_num), levels = sort(as.numeric(unique(fd_cub$hab_num))))
  fd_cub$dim <- factor(as.numeric(fd_cub$dim), levels = sort(as.numeric(unique(fd_cub$dim))))
  ### Create data subset for co-occurence matching data
  pd_hab <- pd[pd$d %in% hab_levels,]
  pd_hab <- left_join(hab_levels_info,pd_hab, by = 'd', multiple = 'all')
  pd_hab$hab_num <- factor(as.numeric(pd_hab$hab_num), levels = sort(as.numeric(unique(pd_hab$hab_num))))
  pd_cub <- pd[pd$d %in% cub_levels,]
  pd_cub <- left_join(cub_levels_info,pd_cub, by = 'd', multiple = 'all')
  pd_cub$hab_num <- factor(as.numeric(pd_cub$hab_num), levels = sort(as.numeric(unique(pd_cub$hab_num))))
  pd_cub$dim <- factor(as.numeric(pd_cub$dim), levels = sort(as.numeric(unique(pd_cub$dim))))
  ### Create data subset for prediction score metrics
  sd_hab <- sd[sd$d %in% hab_levels,]
  sd_hab <- left_join(hab_levels_info,sd_hab, by = 'd', multiple = 'all')
  sd_hab$hab_num <- factor(as.numeric(sd_hab$hab_num), levels = sort(as.numeric(unique(sd_hab$hab_num))))
  sd_cub <- sd[sd$d %in% cub_levels,]
  sd_cub <- left_join(cub_levels_info,sd_cub, by = 'd', multiple = 'all')
  sd_cub$hab_num <- factor(as.numeric(sd_cub$hab_num), levels = sort(as.numeric(unique(sd_cub$hab_num))))
  sd_cub$dim <- factor(as.numeric(sd_cub$dim), levels = sort(as.numeric(unique(sd_cub$dim))))
  resl <- list(fd_hab,pd_hab,sd_hab,fd_cub,pd_cub,sd_cub)
  return(resl)
}

prep_plot_data_cubes <- function(fd,pd,sd) {
  pd$run_nm <- as.factor(pd$run_nm)
  pd$d <- as.factor(pd$d)
  hab_levels <- levels(pd$d)[grep('HabSubsmpl', levels(pd$d))]
  cub_levels <- levels(pd$d)[grep('Dim', levels(pd$d))]
  hab_levels_info <- as.data.frame(cbind(hab_levels, unlist(strsplit(hab_levels, '_'))[2*(1:length(hab_levels))]))
  colnames(hab_levels_info) <- c('d','hab_num')
  cub_levels_info <- as.data.frame(cbind(cub_levels, unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))],
                                         unlist(strsplit(unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))], '_'))[2*(1:length(cub_levels))],
                                         unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))-1],
                                         unlist(strsplit(unlist(strsplit(cub_levels, '__'))[2*(1:length(cub_levels))-1], '_'))[2*(1:length(cub_levels))]))
  colnames(cub_levels_info) <- c('d','dim_lab','dim','hab_num_lab','hab_num')
  ### Create data subset for community metrics
  fd_hab <- fd[fd$d %in% hab_levels,]
  fd_hab <- left_join(hab_levels_info,fd_hab, by = 'd', multiple = 'all')
  fd_hab$hab_num <- factor(as.numeric(fd_hab$hab_num), levels = sort(as.numeric(unique(fd_hab$hab_num))))
  fd_cub <- fd[fd$d %in% cub_levels,]
  fd_cub <- left_join(cub_levels_info, fd_cub, by = 'd', multiple = 'all')
  fd_cub$hab_num <- factor(as.numeric(fd_cub$hab_num), levels = sort(as.numeric(unique(fd_cub$hab_num))))
  fd_cub$dim <- factor(as.numeric(fd_cub$dim), levels = sort(as.numeric(unique(fd_cub$dim))))
  ### Create data subset for co-occurence matching data
  pd_hab <- pd[pd$d %in% hab_levels,]
  pd_hab <- left_join(hab_levels_info,pd_hab, by = 'd', multiple = 'all')
  pd_hab$hab_num <- factor(as.numeric(pd_hab$hab_num), levels = sort(as.numeric(unique(pd_hab$hab_num))))
  pd_cub <- pd[pd$d %in% cub_levels,]
  pd_cub <- left_join(cub_levels_info,pd_cub, by = 'd', multiple = 'all')
  pd_cub$hab_num <- factor(as.numeric(pd_cub$hab_num), levels = sort(as.numeric(unique(pd_cub$hab_num))))
  pd_cub$dim <- factor(as.numeric(pd_cub$dim), levels = sort(as.numeric(unique(pd_cub$dim))))
  ### Create data subset for prediction score metrics
  sd_hab <- sd[sd$d %in% hab_levels,]
  sd_hab <- left_join(hab_levels_info,sd_hab, by = 'd', multiple = 'all')
  sd_hab$hab_num <- factor(as.numeric(sd_hab$hab_num), levels = sort(as.numeric(unique(sd_hab$hab_num))))
  sd_cub <- sd[sd$d %in% cub_levels,]
  sd_cub <- left_join(cub_levels_info,sd_cub, by = 'd', multiple = 'all')
  sd_cub$hab_num <- factor(as.numeric(sd_cub$hab_num), levels = sort(as.numeric(unique(sd_cub$hab_num))))
  sd_cub$dim <- factor(as.numeric(sd_cub$dim), levels = sort(as.numeric(unique(sd_cub$dim))))
  resl <- list(fd_hab,pd_hab,sd_hab,fd_cub,pd_cub,sd_cub)
  return(resl)
}

###########################################################################################
# Load and prep data
get_and_mod_data <- function(parent_set_of_analyses) {
  runs_incl <- list.dirs(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses), full.names = FALSE, recursive = FALSE)
  out <- get_runs(parent_set_of_analyses, runs_incl)
  pd <- out[[1]]
  fd <- out[[2]]
  sd <- out[[3]]
  ###########################################################################################
  ### Add corrected match values
  colnames(pd)
  # All comparisons considered
  pd$eucp_cor <- pd$eucp - (pd$eucshrp - pd$eucshrintp) - (pd$eucintp - pd$eucshrintp) - pd$eucshrintp
  pd$shrp_cor <- pd$shrp - (pd$eucshrp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$intp_cor <- pd$intp - (pd$eucintp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$eucshrp_cor <- pd$eucshrp - pd$eucshrintp
  pd$eucintp_cor <- pd$eucintp - pd$eucshrintp
  pd$shrintp_cor <- pd$shrintp - pd$eucshrintp
  #
  pd$eucn_cor <- pd$eucn - (pd$eucshrn - pd$eucshrintn) - (pd$eucintn - pd$eucshrintn) - pd$eucshrintn
  pd$shrn_cor <- pd$shrn - (pd$eucshrn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$intn_cor <- pd$intn - (pd$eucintn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$eucshrn_cor <- pd$eucshrn - pd$eucshrintn
  pd$eucintn_cor <- pd$eucintn - pd$eucshrintn
  pd$shrintn_cor <- pd$shrintn - pd$eucshrintn
  #####
  # A Resource preference and interactions based matching considered
  pd$eucp_corA <- pd$eucp - pd$eucintp
  pd$intp_corA <- pd$intp - pd$eucintp
  pd$eucintp_corA <- pd$eucintp
  pd$unxpA <- pd$unxp + (pd$shrp - (pd$eucshrp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corA <- pd$eucn - pd$eucintn
  pd$intn_corA <- pd$intn - pd$eucintn
  pd$eucintn_corA <- pd$eucintn
  pd$unxnA <- pd$unxn + (pd$shrn - (pd$eucshrn) - (pd$shrintn - pd$eucshrintn))
  #####
  # B Resource preference and shared resource association based matching considered
  pd$eucp_corB <- pd$eucp - pd$eucshrp
  pd$shrp_corB <- pd$shrp - pd$eucshrp
  pd$eucshrp_corB <- pd$eucshrp
  pd$unxpB <- pd$unxp + (pd$intp - (pd$eucintp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corB <- pd$eucn - pd$eucshrn
  pd$shrn_corB <- pd$shrn - pd$eucshrn
  pd$eucshrn_corB <- pd$eucshrn
  pd$unxnB <- pd$unxn + (pd$intn - (pd$eucintn) - (pd$shrintn - pd$eucshrintn))
  ### Get descriptors for different runs
  run_nm_strs <- c('same0','same5','same10','diffL','knorm','inv')
  pd_nw <- c()
  fd_nw <- c()
  sd_nw <- c()
  #i <- run_nm_strs[1]
  for (i in run_nm_strs) {
    pd_nw <- rbind(pd_nw, pd[grepl(i, pd$run_nm),])
    fd_nw <- rbind(fd_nw, fd[grepl(i, fd$run_nm),])
    sd_nw <- rbind(sd_nw, sd[grepl(i, sd$run_nm),])
  }
  
  pd <- pd_nw
  fd <- fd_nw
  sd <- sd_nw
  # Get better axis names
  #param_ref <- get_param_ref()
  ### Prep data for subdividing and better plotting
  resl <- prep_plot_data(fd,pd,sd)
  fd_hab <- resl[[1]]
  pd_hab <- resl[[2]] 
  sd_hab <- resl[[3]] 
  fd_cub <- resl[[4]]
  pd_cub <- resl[[5]]
  sd_cub <- resl[[6]]
  
  fd_hab$rscdist <- NA 
  pd_hab$rscdist <- NA
  sd_hab$rscdist <- NA
  fd_cub$rscdist <- NA
  pd_cub$rscdist <- NA
  sd_cub$rscdist <- NA

  #j <- run_nm_strs[4]
  for (j in run_nm_strs) {
    if (sum(grepl(j , fd_hab$run_nm)) > 0) {
      fd_hab[grepl(j , fd_hab$run_nm),]$rscdist <- j
      pd_hab[grepl(j , pd_hab$run_nm),]$rscdist <- j
      sd_hab[grepl(j , sd_hab$run_nm),]$rscdist <- j
      #fd_cub[grepl(j , fd_cub$run_nm),]$rscdist <- j
      #pd_cub[grepl(j , pd_cub$run_nm),]$rscdist <- j
      #sd_cub[grepl(j , sd_cub$run_nm),]$rscdist <- j
    }
  }
  
  fd_hab$hab_num <- as.numeric(as.character(fd_hab$hab_num))
  pd_hab$hab_num <- as.numeric(as.character(pd_hab$hab_num))
  sd_hab$hab_num <- as.numeric(as.character(sd_hab$hab_num))
  out <- list(fd_hab, pd_hab, sd_hab, fd_cub, pd_cub, sd_cub)
  return(out)
}


###########################################################################################
# Load and prep data
get_and_mod_data_cubes <- function(parent_set_of_analyses) {
  runs_incl <- list.dirs(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses), full.names = FALSE, recursive = FALSE)
  out <- get_runs(parent_set_of_analyses, runs_incl)
  pd <- out[[1]]
  fd <- out[[2]]
  sd <- out[[3]]
  #unique(sd$run_nm)
  ###########################################################################################
  ### Add corrected match values
  colnames(pd)
  # All comparisons considered
  pd$eucp_cor <- pd$eucp - (pd$eucshrp - pd$eucshrintp) - (pd$eucintp - pd$eucshrintp) - pd$eucshrintp
  pd$shrp_cor <- pd$shrp - (pd$eucshrp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$intp_cor <- pd$intp - (pd$eucintp - pd$eucshrintp) - (pd$shrintp - pd$eucshrintp) - pd$eucshrintp
  pd$eucshrp_cor <- pd$eucshrp - pd$eucshrintp
  pd$eucintp_cor <- pd$eucintp - pd$eucshrintp
  pd$shrintp_cor <- pd$shrintp - pd$eucshrintp
  #
  pd$eucn_cor <- pd$eucn - (pd$eucshrn - pd$eucshrintn) - (pd$eucintn - pd$eucshrintn) - pd$eucshrintn
  pd$shrn_cor <- pd$shrn - (pd$eucshrn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$intn_cor <- pd$intn - (pd$eucintn - pd$eucshrintn) - (pd$shrintn - pd$eucshrintn) - pd$eucshrintn
  pd$eucshrn_cor <- pd$eucshrn - pd$eucshrintn
  pd$eucintn_cor <- pd$eucintn - pd$eucshrintn
  pd$shrintn_cor <- pd$shrintn - pd$eucshrintn
  #####
  # A Resource preference and interactions based matching considered
  pd$eucp_corA <- pd$eucp - pd$eucintp
  pd$intp_corA <- pd$intp - pd$eucintp
  pd$eucintp_corA <- pd$eucintp
  pd$unxpA <- pd$unxp + (pd$shrp - (pd$eucshrp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corA <- pd$eucn - pd$eucintn
  pd$intn_corA <- pd$intn - pd$eucintn
  pd$eucintn_corA <- pd$eucintn
  pd$unxnA <- pd$unxn + (pd$shrn - (pd$eucshrn) - (pd$shrintn - pd$eucshrintn))
  #####
  # B Resource preference and shared resource association based matching considered
  pd$eucp_corB <- pd$eucp - pd$eucshrp
  pd$shrp_corB <- pd$shrp - pd$eucshrp
  pd$eucshrp_corB <- pd$eucshrp
  pd$unxpB <- pd$unxp + (pd$intp - (pd$eucintp) - (pd$shrintp - pd$eucshrintp))
  #
  pd$eucn_corB <- pd$eucn - pd$eucshrn
  pd$shrn_corB <- pd$shrn - pd$eucshrn
  pd$eucshrn_corB <- pd$eucshrn
  pd$unxnB <- pd$unxn + (pd$intn - (pd$eucintn) - (pd$shrintn - pd$eucshrintn))
  ### Get descriptors for different runs
  run_nm_strs <- c('same0','same5','same10','diffL','knorm','inv') ############# Changed this!!!!
  pd_nw <- c()
  fd_nw <- c()
  sd_nw <- c()
  #i <- run_nm_strs[1]
  for (i in run_nm_strs) {
    pd_nw <- rbind(pd_nw, pd[grepl(i, pd$run_nm),])
    fd_nw <- rbind(fd_nw, fd[grepl(i, fd$run_nm),])
    sd_nw <- rbind(sd_nw, sd[grepl(i, sd$run_nm),])
  }
  
  pd <- pd_nw
  fd <- fd_nw
  sd <- sd_nw
  # Get better axis names
  #param_ref <- get_param_ref()
  ### Prep data for subdividing and better plotting
  resl <- prep_plot_data_cubes(fd,pd,sd)
  #fd_hab <- resl[[1]]
  #pd_hab <- resl[[2]] 
  #sd_hab <- resl[[3]] 
  fd_cub <- resl[[4]]
  pd_cub <- resl[[5]]
  sd_cub <- resl[[6]]
  
  #fd_hab$rscdist <- NA 
  #pd_hab$rscdist <- NA
  #sd_hab$rscdist <- NA
  fd_cub$rscdist <- NA
  pd_cub$rscdist <- NA
  sd_cub$rscdist <- NA
  
  #j <- run_nm_strs[2]
  for (j in run_nm_strs) {
    if (sum(grepl(j , fd_cub$run_nm)) > 0) {
      #fd_hab[grepl(j , fd_hab$run_nm),]$rscdist <- j
      #pd_hab[grepl(j , pd_hab$run_nm),]$rscdist <- j
      #sd_hab[grepl(j , sd_hab$run_nm),]$rscdist <- j
      fd_cub[grepl(j , fd_cub$run_nm),]$rscdist <- j
      pd_cub[grepl(j , pd_cub$run_nm),]$rscdist <- j
      sd_cub[grepl(j , sd_cub$run_nm),]$rscdist <- j
    }
  }
  
  #fd_hab$hab_num <- as.numeric(as.character(fd_hab$hab_num))
  #pd_hab$hab_num <- as.numeric(as.character(pd_hab$hab_num))
  #sd_hab$hab_num <- as.numeric(as.character(sd_hab$hab_num))
  out <- list(fd_cub, pd_cub, sd_cub) #fd_hab, pd_hab, sd_hab, 
  return(out)
}


###########################################################################################
# Load and abd data
get_abddata <- function(parent_set_of_analyses) {

  runs_incl <- list.dirs(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses), full.names = FALSE, recursive = FALSE)
  #j <- runs_incl[1]
  #i <- 1
  abd_dat <- c()
  for (j in runs_incl) {
    parentdirs <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,"/",j,sep=''), full.names=TRUE) #, pattern = "Z"
    parentdirnms <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,"/",j,sep='')) #, pattern = "Z"
    for (i in 1:length(parentdirnms)) {
      if (file.exists(paste0(parentdirs[i],'/infm_res.csv')) & file.exists(paste0(parentdirs[i],'/full_res.csv')) & file.exists(paste0(parentdirs[i],'/pred_res.csv'))) {
        print(parentdirnms[i])
        abd_dat_tmp <- read.csv(paste0(parentdirs[i],'/abd_res_comp.csv'))[-1]
        abd_dat <- rbind(abd_dat,cbind(abd_dat_tmp, run_nm = parentdirnms[i]))      
      }
      else {
        print(paste0('Required files do not exist in: ',parentdirnms[i]))
      }
    }
  }
  return(abd_dat)
}


