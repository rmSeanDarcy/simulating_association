###########################################################################################
### SAMSARA - Figure 2 and plots for related Supporting information                     ###
###########################################################################################

###########################################################################################
### Fig 2A
get_correlation_driver_map <- function(infm, plotdims, titltxt) {
  
  dx <- infm[infm$d == 'HabSubsmpl_25',]
  dx$corbin <- 'Not associated'
  dx$corbin[dx$cooc > 0] <- 'Positive'
  dx$corbin[dx$cooc < 0] <- 'Negative'
  dx$corbin <- factor(dx$corbin, levels = c('Positive','Negative','Not associated'))
  dx$shrbin <- 'No relationship'
  dx$shrbin[dx$shrd_rsc_pos == 1] <- 'Shared'
  dx$shrbin[dx$shrd_rsc_neg == 1] <- 'Inverse'
  dx$shrbin <- factor(dx$shrbin, levels = c('No relationship','Shared','Inverse'))
  dx$int <- dx$int*-1
  dx$euc <- dx$euc*-1
  
  ggplot(dx, aes(x = int, y = euc, col = corbin)) + ylim(-1.5,0) +
    geom_point(size = 2,alpha = 0.3) + xlab(TeX('Interaction coefficient $\\alpha_{ij}$')) + ylab('Environmental preference similarity $E_{ij}$') + theme_bw() +
    scale_color_manual(values = c('red','blue','grey')) + labs(color = 'Species pair') + #, shape = 'Resource\nrelationship'
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}

###########################################################################################
### Fig. 2 B
plot_match_allruns_eucintp <- function(pd_hab,label_mode,titltxt,hnum,plotdims) {
  
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA')]
  pd_hab$run_nm <- pd_hab$simulation
  rownames(pd_haby) <- pd_hab$run_nm
  # intp_corA = Precision for co-occurrences matching mutualism (without intersection with preference similarity)
  # eucp_corA = Precision for co-occurrences matching preference similarity (without intersection with mutualsim)
  # eucintp_corA = Precision for co-occurrences matching intersection preference similarity and mutualsim
  # unxpA = Remaining 'precision' for co-occurrences that don't match either
  ### For individual driver runs the unexplained needs to be adjusted!
  pd_haby[rownames(pd_haby) == 'c2',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'c2',]$eucp
  pd_haby[rownames(pd_haby) == 'c2',]$eucp_corA <- pd_hab[rownames(pd_haby) == 'c2',]$eucp
  pd_haby[rownames(pd_haby) == 'c1',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'c1',]$intp
  pd_haby[rownames(pd_haby) == 'c1',]$intp_corA <- pd_hab[rownames(pd_haby) == 'c1',]$intp
  pd_haby[rownames(pd_haby) == 'c1',]$eucintp_corA <- 0
  pd_haby[rownames(pd_haby) == 'c1',]$eucp_corA <- 0
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  colnames(pd) <- c('matches','rscdist','perc')
  comp <- as.data.frame(cbind(c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$rscdist,pd_hab$run_nm),]$totp
  pd$perc_n <- pd$perc*pd_hab[match(pd$rscdist,pd_hab$run_nm),]$totp
  
  sumsx <- c()
  for (i in unique(pd$rscdist)) {
    sumsx <- c(sumsx,sum(pd[pd$rscdist == i,]$perc_n))
  }
  
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  pd$rscdist <- factor(pd$rscdist, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5')) #, 'cf'
  levels(pd$rscdist)[levels(pd$rscdist)=='c1'] <- 'Interactions only'
  levels(pd$rscdist)[levels(pd$rscdist)=='c2'] <- 'Environment only'
  levels(pd$rscdist)[levels(pd$rscdist)=='s1'] <- 'Environment & Interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='s4'] <- 'Sampling noise'
  levels(pd$rscdist)[levels(pd$rscdist)=='s2'] <- 'Equal interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='s5'] <- 'Random extinctions (20%)'
  levels(pd$rscdist)[levels(pd$rscdist)=='s3'] <- 'Resources anti-correlated'
  
  ggplot(pd, aes(y=perc_n, x=rscdist, fill = lab)) + #, col = lab
    geom_bar(position="stack", stat="identity", lwd = 0) +
    theme_bw() + xlab('Scenarios') + ylab('Positive\nassociations') +
    scale_fill_manual(values = c(rest_col,euc_col,eucint_col,int_col), name = 'Association\nmatches:') +
    #scale_color_manual(values = c(darken(rest_col, darkfac),darken(euc_col, darkfac),darken(eucint_col, darkfac),darken(int_col, darkfac)), name = 'Association\nmatches:') +
    
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_col_pattern(aes(pattern = lab),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    #theme(axis.text.x=element_text(angle=45)) +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    coord_flip() + scale_x_discrete(limits = rev(levels(pd$rscdist))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Fig. 2C
plot_match_allruns_eucintn <- function(pd_hab,label_mode,titltxt,hnum,plotdims) {
  
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA')]
  pd_hab$run_nm <- pd_hab$simulation
  rownames(pd_haby) <- pd_hab$run_nm
  ### For individual driver runs the unexplained needs to be adjusted!
  pd_haby[rownames(pd_haby) == 'c2',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'c2',]$eucn
  pd_haby[rownames(pd_haby) == 'c2',]$eucn_corA <- pd_hab[rownames(pd_haby) == 'c2',]$eucn
  pd_haby[rownames(pd_haby) == 'c1',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'c1',]$intn
  pd_haby[rownames(pd_haby) == 'c1',]$intn_corA <- pd_hab[rownames(pd_haby) == 'c1',]$intn
  pd_haby[rownames(pd_haby) == 'c1',]$eucintn_corA <- 0
  pd_haby[rownames(pd_haby) == 'c1',]$eucn_corA <- 0
  
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  
  colnames(pd) <- c('matches','rscdist','perc')
  comp <- as.data.frame(cbind(c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$rscdist,pd_hab$run_nm),]$totn
  pd$perc_n <- pd$perc*pd_hab[match(pd$rscdist,pd_hab$run_nm),]$totn
  
  sumsx <- c()
  for (i in unique(pd$rscdist)) {
    sumsx <- c(sumsx,sum(pd[pd$rscdist == i,]$perc_n))
  }
  
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  pd$rscdist <- factor(pd$rscdist, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5')) #, 'cf'
  levels(pd$rscdist)[levels(pd$rscdist)=='c1'] <- 'Interactions only'
  levels(pd$rscdist)[levels(pd$rscdist)=='c2'] <- 'Environment only'
  levels(pd$rscdist)[levels(pd$rscdist)=='s1'] <- 'Environment & Interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='s4'] <- 'Sampling noise'
  levels(pd$rscdist)[levels(pd$rscdist)=='s2'] <- 'Equal interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='s5'] <- 'Random extinctions (20%)'
  levels(pd$rscdist)[levels(pd$rscdist)=='s3'] <- 'Resources anti-correlated'
  
  ggplot(pd, aes(y=perc_n, x=rscdist, fill = lab)) + #, col = lab
    geom_bar(position="stack", stat="identity", lwd = 0) +
    theme_bw() + xlab('Scenarios') + ylab('Negative\nassociations') +
    scale_fill_manual(values = c(rest_col,euc_col,eucint_col,int_col), name = 'Association\nmatches:') +
    #scale_color_manual(values = c(darken(rest_col, darkfac),darken(euc_col, darkfac),darken(eucint_col, darkfac),darken(int_col, darkfac)), name = 'Association\nmatches:') +
    
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_col_pattern(aes(pattern = lab),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    #theme(axis.text.x=element_text(angle=45)) +
    coord_flip() + scale_x_discrete(limits = rev(levels(pd$rscdist))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}

###########################################################################################
##### Recovery of intial settings (i.e. recall or sensitivity)
# params = c('intp','totp')
# total = c(4,0,4,11,4,4,4)
plot_recovery_int <- function(pd_hab, total, params) {
  pd <- pd_hab[,colnames(pd_hab) %in% c(params)]
  pd <- as.data.frame(sapply(pd, as.numeric))
  colnames(pd) <- c('match','total_cooc')
  pd$scenario <- pd_hab$simulation
  pd$total_match <- pd$match*pd$total_cooc
  pd$total_settings <- total
  pd$total_nonmatch <- pd$total_settings - pd$total_match
  pdx <- pd[,colnames(pd) %in% c('total_match', 'total_nonmatch')]
  rownames(pdx) <- pd$scenario 
  pd <- as.data.frame(melt(t(pdx)))
  colnames(pd) <- c('matches', 'scenario', 'abs')
  pd$scenario <- factor(pd$scenario, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5'))
  pd$matches <- factor(pd$matches, levels = c('total_nonmatch','total_match'))
  ggplot(pd, aes(x = scenario, y = abs, fill = matches)) +
    geom_bar(position="stack", stat="identity", lwd = 0) + theme_bw() +
    scale_fill_manual(values = c(rest_col,int_col), name = 'Association\nmatches') +
    xlab('Scenarios') + ylab('Potential\nmatches') +
    coord_flip() + scale_x_discrete(limits = rev(levels(pd$scenario))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}  
### Euclidean distance
# total = c(0,11,11,11,11,11,11)
# params = c('eucp','totp')
# params = c('eucn','totn')
plot_recovery_euc <- function(pd_hab, total, params) {
  pd <- pd_hab[,colnames(pd_hab) %in% c(params)]
  pd <- as.data.frame(sapply(pd, as.numeric))
  # Have to set the control for environment to zero here -> In the simulations I give species preferences, but overrider their K_{i,h} later on. When analysing co-occurrecnces by chance their shared preference can still match co-occurrecnes though
  if ('eucp' %in% colnames(pd)) {
    pd[pd_hab$simulation == 'c1',c(1,2)] <- 0 
  } else if ('eucn' %in% colnames(pd)) {
    pd[pd_hab$simulation == 'c1',c(1,2)] <- 0
    pd[pd_hab$simulation == 'c2',c(1,2)] <- 0
  }
  colnames(pd) <- c('match','total_cooc')
  pd$scenario <- pd_hab$simulation
  pd$total_match <- pd$match*pd$total_cooc
  pd$total_settings <- total
  pd$total_nonmatch <- pd$total_settings - pd$total_match
  pdx <- pd[,colnames(pd) %in% c('total_match', 'total_nonmatch')]
  rownames(pdx) <- pd$scenario 
  pd <- as.data.frame(melt(t(pdx)))
  colnames(pd) <- c('matches', 'scenario', 'abs')
  pd$scenario <- factor(pd$scenario, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5'))
  pd$matches <- factor(pd$matches, levels = c('total_nonmatch','total_match'))
  ggplot(pd, aes(x = scenario, y = abs, fill = matches)) +
    geom_bar(position="stack", stat="identity", lwd = 0) + theme_bw() +
    scale_fill_manual(values = c(rest_col,euc_col), name = 'Association\nmatches') +
    xlab('Scenarios') + ylab('Potential\nmatches') +
    coord_flip() + scale_x_discrete(limits = rev(levels(pd$scenario))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}  

