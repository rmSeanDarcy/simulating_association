###########################################################################################
### SAMSARA - Figure 2 and plots for related Supporting information                     ###
###########################################################################################
sml_font <- 8
med_font <- 10
big_font <- 11
font_size_control <- theme(text=element_text(size=sml_font), #change font size of all text
                           axis.text=element_text(size=sml_font), #change font size of axis text
                           axis.title=element_text(size=med_font), #change font size of axis titles
                           plot.title=element_text(size=big_font), #change font size of plot title
                           legend.text=element_text(size=med_font), #change font size of legend text
                           legend.title=element_text(size=med_font)) #change font size of legend title   

###########################################################################################
###########################################################################################
##### Main manuscript: Fig 2
###########################################################################################
###########################################################################################

###########################################################################################
### Matching co-occurrences in a simplified system with initial drivers in various scenarios
###########################################################################################
##### Fig. 2 B) 
###########################################################################################
### Plot for negative co-occurrences
plot_match_allruns_eucintp <- function(pd_hab,label_mode,titltxt,hnum,plotdims) {
  
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA')]
  rownames(pd_haby) <- pd_hab$run_nm
  # intp_corA = Precision for co-occurrences matching mutualism (without intersection with preference similarity)
  # eucp_corA = Precision for co-occurrences matching preference similarity (without intersection with mutualsim)
  # eucintp_corA = Precision for co-occurrences matching intersection preference similarity and mutualsim
  # unxpA = Remaining 'precision' for co-occurrences that don't match either
  ### For individual driver runs the unexplained needs to be adjusted!
  #pd_haby[rownames(pd_haby) == 'c2',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'c2',]$eucp
  #pd_haby[rownames(pd_haby) == 'c2',]$eucp_corA <- pd_hab[rownames(pd_haby) == 'c2',]$eucp
  #pd_haby[rownames(pd_haby) == 'c1',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'c1',]$intp
  #pd_haby[rownames(pd_haby) == 'c1',]$intp_corA <- pd_hab[rownames(pd_haby) == 'c1',]$intp
  #pd_haby[rownames(pd_haby) == 'c1',]$eucintp_corA <- 0 
  #pd_haby[rownames(pd_haby) == 'c1',]$eucp_corA <- 0
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

  #pd$rscdist <- factor(pd$rscdist, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5')) #, 'cf'
  #levels(pd$rscdist)[levels(pd$rscdist)=='c1'] <- 'Interactions only'  
  #levels(pd$rscdist)[levels(pd$rscdist)=='c2'] <- 'Environment only'
  #levels(pd$rscdist)[levels(pd$rscdist)=='s1'] <- 'Environment & Interactions'
  #levels(pd$rscdist)[levels(pd$rscdist)=='s4'] <- 'Sampling noise'
  #levels(pd$rscdist)[levels(pd$rscdist)=='s2'] <- 'Equal interactions'
  #levels(pd$rscdist)[levels(pd$rscdist)=='s5'] <- 'Random extinctions (20%)'
  #levels(pd$rscdist)[levels(pd$rscdist)=='s3'] <- 'Resources anti-correlated'  

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
### Plot for negative co-occurrences
plot_match_allruns_eucintn <- function(pd_hab,label_mode,titltxt,hnum,plotdims) {
  
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA')]
  rownames(pd_haby) <- pd_hab$run_nm
  ### For individual driver runs the unexplained needs to be adjusted!
  # pd_haby[rownames(pd_haby) == 'c2',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'c2',]$eucn
  # pd_haby[rownames(pd_haby) == 'c2',]$eucn_corA <- pd_hab[rownames(pd_haby) == 'c2',]$eucn
  # pd_haby[rownames(pd_haby) == 'c1',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'c1',]$intn
  # pd_haby[rownames(pd_haby) == 'c1',]$intn_corA <- pd_hab[rownames(pd_haby) == 'c1',]$intn
  # pd_haby[rownames(pd_haby) == 'c1',]$eucintn_corA <- 0 
  # pd_haby[rownames(pd_haby) == 'c1',]$eucn_corA <- 0
  
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
  
  # pd$rscdist <- factor(pd$rscdist, levels = c('s1', 'c1', 'c2', 's2', 's3', 's4', 's5')) #, 'cf'
  # levels(pd$rscdist)[levels(pd$rscdist)=='c1'] <- 'Interactions only'  
  # levels(pd$rscdist)[levels(pd$rscdist)=='c2'] <- 'Environment only'
  # levels(pd$rscdist)[levels(pd$rscdist)=='s1'] <- 'Environment & Interactions'
  # levels(pd$rscdist)[levels(pd$rscdist)=='s4'] <- 'Sampling noise'
  # levels(pd$rscdist)[levels(pd$rscdist)=='s2'] <- 'Equal interactions'
  # levels(pd$rscdist)[levels(pd$rscdist)=='s5'] <- 'Random extinctions (20%)'
  # levels(pd$rscdist)[levels(pd$rscdist)=='s3'] <- 'Resources anti-correlated'  
  
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




