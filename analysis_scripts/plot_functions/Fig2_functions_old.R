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
##### Interaction x environmental preference plots
###########################################################################################
##### Fig. 2 A) 
###########################################################################################
get_correlation_driver_map <- function(infmat_res_comp, plotdims, titltxt) {
  
  dx <- infmat_res_comp[infmat_res_comp$d == 'HabSubsmpl_25',]
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
### Matching co-occurrences in a simplified system with initial drivers in various scenarios
###########################################################################################
##### Fig. 2 B) 
###########################################################################################
### Plot for negative co-occurrences
plot_match_allruns_eucintp <- function(pd_hab,label_mode,titltxt,hnum,plotdims) {
  
  runs <- c('basic_F2_standard_same0__X__Z','null_int_F2_intonly_same0__X__Z','null_rsc_F2_envonly_same0__X__Z','basic_F2_noise_same0__X__Znoise0.25',
            'basic_F2_randomextinct_same0__X__Z','basic_F2_equalint_same0__X__Z','basic_F2_invdist_invdist0__X__Z')
  pd_hab <- pd_hab[pd_hab$run_nm %in% runs,]
  pd_hab[grepl('noise',pd_hab$run_nm),]$rscdist <- 'noise' 
  pd_hab[grepl('equalint',pd_hab$run_nm),]$rscdist <- 'equalint' 
  pd_hab[grepl('envonly',pd_hab$run_nm),]$rscdist <- 'envonly' 
  pd_hab[grepl('randomextinct',pd_hab$run_nm),]$rscdist <- 'randomextinct' 
  pd_hab[grepl('intonly',pd_hab$run_nm),]$rscdist <- 'intonly' 
  ### Subset for required columns
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA')]
  rownames(pd_haby) <- pd_hab$rscdist
  # intp_corA = Precision for co-occurrences matching mutualism (without intersection with preference similarity)
  # eucp_corA = Precision for co-occurrences matching preference similarity (without intersection with mutualsim)
  # eucintp_corA = Precision for co-occurrences matching intersection preference similarity and mutualsim
  # unxpA = Remaining 'precision' for co-occurrences that don't match either
  ### For individual driver runs the unexplained needs to be adjusted!
  pd_haby[rownames(pd_haby) == 'envonly',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'envonly',]$eucp
  pd_haby[rownames(pd_haby) == 'envonly',]$eucp_corA <- pd_hab[rownames(pd_haby) == 'envonly',]$eucp
  pd_haby[rownames(pd_haby) == 'intonly',]$unxpA <- 1-pd_hab[rownames(pd_haby) == 'intonly',]$intp
  pd_haby[rownames(pd_haby) == 'intonly',]$intp_corA <- pd_hab[rownames(pd_haby) == 'intonly',]$intp
  pd_haby[rownames(pd_haby) == 'intonly',]$eucintp_corA <- 0 
  pd_haby[rownames(pd_haby) == 'intonly',]$eucp_corA <- 0
  
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  colnames(pd) <- c('matches','rscdist','perc')
  comp <- as.data.frame(cbind(c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$rscdist,pd_hab$rscdist),]$totp
  pd$perc_n <- pd$perc*pd_hab[match(pd$rscdist,pd_hab$rscdist),]$totp
  
  sumsx <- c()
  for (i in unique(pd$rscdist)) {
    sumsx <- c(sumsx,sum(pd[pd$rscdist == i,]$perc_n))
  }
  
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"

  pd$rscdist <- factor(pd$rscdist, levels = c('same0', 'intonly', 'envonly', 'equalint', 'inv', 'noise', 'randomextinct')) #, 'cf'
  levels(pd$rscdist)[levels(pd$rscdist)=='intonly'] <- 'Interactions only'  
  levels(pd$rscdist)[levels(pd$rscdist)=='envonly'] <- 'Environment only'
  levels(pd$rscdist)[levels(pd$rscdist)=='same0'] <- 'Environment & Interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='noise'] <- 'Sampling noise'
  levels(pd$rscdist)[levels(pd$rscdist)=='equalint'] <- 'Equal interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='randomextinct'] <- 'Random extinctions (20%)'
  levels(pd$rscdist)[levels(pd$rscdist)=='inv'] <- 'Resources anti-correlated'  

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
  
  runs <- c('basic_F2_standard_same0__X__Z','null_int_F2_intonly_same0__X__Z','null_rsc_F2_envonly_same0__X__Z','basic_F2_noise_same0__X__Znoise0.25',
            'basic_F2_randomextinct_same0__X__Z','basic_F2_equalint_same0__X__Z','basic_F2_invdist_invdist0__X__Z')
  pd_hab <- pd_hab[pd_hab$run_nm %in% runs,]
  pd_hab[grepl('noise',pd_hab$run_nm),]$rscdist <- 'noise'       
  pd_hab[grepl('equalint',pd_hab$run_nm),]$rscdist <- 'equalint'        
  pd_hab[grepl('envonly',pd_hab$run_nm),]$rscdist <- 'envonly'            
  pd_hab[grepl('randomextinct',pd_hab$run_nm),]$rscdist <- 'randomextinct' 
  pd_hab[grepl('intonly',pd_hab$run_nm),]$rscdist <- 'intonly' 
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA')]
  rownames(pd_haby) <- pd_hab$rscdist
  ### For individual driver runs the unexplained needs to be adjusted!
  pd_haby[rownames(pd_haby) == 'envonly',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'envonly',]$eucn
  pd_haby[rownames(pd_haby) == 'envonly',]$eucn_corA <- pd_hab[rownames(pd_haby) == 'envonly',]$eucn
  pd_haby[rownames(pd_haby) == 'intonly',]$unxnA <- 1-pd_hab[rownames(pd_haby) == 'intonly',]$intn
  pd_haby[rownames(pd_haby) == 'intonly',]$intn_corA <- pd_hab[rownames(pd_haby) == 'intonly',]$intn
  pd_haby[rownames(pd_haby) == 'intonly',]$eucintn_corA <- 0 
  pd_haby[rownames(pd_haby) == 'intonly',]$eucn_corA <- 0
  
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  
  colnames(pd) <- c('matches','rscdist','perc')
  comp <- as.data.frame(cbind(c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$rscdist,pd_hab$rscdist),]$totn
  pd$perc_n <- pd$perc*pd_hab[match(pd$rscdist,pd_hab$rscdist),]$totn
  
  sumsx <- c()
  for (i in unique(pd$rscdist)) {
    sumsx <- c(sumsx,sum(pd[pd$rscdist == i,]$perc_n))
  }
  
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  pd$rscdist <- factor(pd$rscdist, levels = c('same0', 'intonly', 'envonly', 'equalint', 'inv', 'noise', 'randomextinct')) #, 'cf'
  levels(pd$rscdist)[levels(pd$rscdist)=='intonly'] <- 'Interactions only'  
  levels(pd$rscdist)[levels(pd$rscdist)=='envonly'] <- 'Environment only'
  levels(pd$rscdist)[levels(pd$rscdist)=='same0'] <- 'Environment & Interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='noise'] <- 'Sampling noise'
  levels(pd$rscdist)[levels(pd$rscdist)=='equalint'] <- 'Equal interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='randomextinct'] <- 'Random extinctions (20%)'
  levels(pd$rscdist)[levels(pd$rscdist)=='inv'] <- 'Resources anti-correlated'  
  
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
###########################################################################################
### General function takes different inputs
#params <- c('intp','totp')
plot_recovery <- function(pd_hab, total, params, plt_col, plt_col_lab, titltxt, plotdims) {
  
  runs <- c('basic_F2_standard_same0__X__Z','null_int_F2_intonly_same0__X__Z','null_rsc_F2_envonly_same0__X__Z','basic_F2_noise_same0__X__Znoise0.25',
            'basic_F2_randomextinct_same0__X__Z','basic_F2_equalint_same0__X__Z','basic_F2_invdist_invdist0__X__Z')
  pd_hab <- pd_hab1[pd_hab1$run_nm %in% runs,]
  pd_hab[grepl('noise',pd_hab$run_nm),]$rscdist <- 'noise' 
  pd_hab[grepl('equalint',pd_hab$run_nm),]$rscdist <- 'equalint' 
  pd_hab[grepl('envonly',pd_hab$run_nm),]$rscdist <- 'envonly' 
  pd_hab[grepl('randomextinct',pd_hab$run_nm),]$rscdist <- 'randomextinct' 
  pd_hab[grepl('intonly',pd_hab$run_nm),]$rscdist <- 'intonly'
  pd <- as.data.frame(cbind(pd_hab[,colnames(pd_hab) == params[1]],
                            pd_hab[,colnames(pd_hab) == params[2]],
                            pd_hab[,colnames(pd_hab) == 'rscdist']))
  pd[,1:2] <- sapply(pd[,1:2], as.numeric)
  colnames(pd) <- c('match','total_cooc','rscdist')
  pd$total_match <- pd$match*pd$total_cooc
  pd$total_settings <- total
  pd$total_nonmatch <- pd$total_settings - pd$total_match
  pdx <- pd[,colnames(pd) %in% c('total_match', 'total_nonmatch')]
  rownames(pdx) <- pd$rscdist 
  pd <- as.data.frame(melt(t(pdx)))
  colnames(pd) <- c('matches', 'rscdist', 'abs')
  pd$rscdist <- factor(pd$rscdist, levels = c('same0', 'intonly', 'envonly', 'equalint', 'inv', 'noise', 'randomextinct')) #, 'cf'
  levels(pd$rscdist)[levels(pd$rscdist)=='intonly'] <- 'Interactions only'  
  levels(pd$rscdist)[levels(pd$rscdist)=='envonly'] <- 'Environment only'
  levels(pd$rscdist)[levels(pd$rscdist)=='same0'] <- 'Environment & Interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='noise'] <- 'Sampling noise'
  levels(pd$rscdist)[levels(pd$rscdist)=='equalint'] <- 'Equal interactions'
  levels(pd$rscdist)[levels(pd$rscdist)=='randomextinct'] <- 'Random extinctions (20%)'
  levels(pd$rscdist)[levels(pd$rscdist)=='inv'] <- 'Resources anti-correlated'  
  
  if (params[1] == 'eucn') {
    pd[pd$rscdist == 'Environment only' & pd$matches == 'total_match',]$abs <- 0
    pd[pd$rscdist == 'Environment only' & pd$matches == 'total_nonmatch',]$abs <- 11
  }
  
  pd$matches <- factor(pd$matches, levels = c('total_nonmatch','total_match'))
  ggplot(pd, aes(x = rscdist, y = abs, fill = matches)) +
    geom_bar(position="stack", stat="identity", lwd = 0) + theme_bw() +
    scale_fill_manual(values = c(rest_col,plt_col), name = 'Association\nmatches', labels = c('Unmatched',plt_col_lab)) +
    xlab('Scenarios') + ylab('Potential\nmatches') +
    coord_flip() + scale_x_discrete(limits = rev(levels(pd$rscdist))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}  


###########################################################################################
###########################################################################################
##### Fig.2: Supporting information                                                   #####
###########################################################################################
###########################################################################################

###########################################################################################
##### SFig2_matches_vs_chance                                                         #####
###########################################################################################
### Plot intersections of Venn diagrams random data
plt_random_euler <- function(prob_int, prob_euc) {
  intersecting_inside <- euler(c("A" = 1-prob_euc-prob_int+(prob_int*prob_euc) , "B" = 0, "C" = 0,
                                 "A&B" = prob_euc, "A&C" = prob_int, "B&C" = 0,
                                 "A&B&C" = prob_int*prob_euc))
  plot(intersecting_inside,
       fills = list(fill = c(NA,euc_col,int_col)), quantities = list(type = "counts", round=2),
       legend = list(side = "right", labels = c("All associations","Env. preference","Interactions")))
}
### Plot intersections of Venn diagrams actual data positive
plt_euler_pos <- function(pd_hab) {
  pd_hab <- pd_hab[pd_hab$hab_num == 25,]
  pd <- pd_hab[,colnames(pd_hab) %in% c('eucp_corA','intp_corA','eucintp_corA')]
  pd_perc <- round(pd,2)
  pd_perc$totp_subtr <- 1 - sum(pd_perc)
  pd_perc$totp <- 1
  values_perc <- c("A" = pd_perc$totp_subtr, "B" = 0, "C" = 0, 
                   "A&B" = pd_perc$eucp_corA,"A&C" = pd_perc$intp_corA, "B&C" = 0,
                   "A&B&C" = pd_perc$eucintp_corA)
  euler_plot_perc <- euler(values_perc)
  plot(euler_plot_perc, shape = "ellipse", main = list(label = paste0("Total negative associations = ",pd_hab$totn), cex = 0.9),
       fills = c(NA, euc_col, int_col), quantities = list(type = "counts",font = 3, round=2),
       legend = list(side = "right", labels = c("All associations","Env. preference","Interactions")))    
}
### Plot intersections of Venn diagrams actual data negative
plt_euler_neg <- function(pd_hab) {
  pd_hab <- pd_hab[pd_hab$hab_num == 25,]
  pd <- pd_hab[,colnames(pd_hab) %in% c('eucn_corA','intn_corA','eucintn_corA')]
  pd_perc <- round(pd,2)
  pd_perc$totn_subtr <- 1 - sum(pd_perc)
  pd_perc$totn <- 1
  values_perc <- c("A" = pd_perc$totn_subtr, "B" = 0, "C" = 0, 
                   "A&B" = pd_perc$eucn_corA,"A&C" = pd_perc$intn_corA, "B&C" = 0,
                   "A&B&C" = pd_perc$eucintn_corA)
  euler_plot_perc <- euler(values_perc)
  plot(euler_plot_perc, shape = "ellipse", main = list(label = paste0("Total negative associations = ",pd_hab$totn), cex = 0.9),
       fills = c(NA, euc_col, int_col), quantities = list(type = "counts",font = 3, round=2),
       legend = list(side = "right", labels = c("All associations","Env. preference","Interactions")))    
}
###########################################################################################
plt_nonrdm_pos <- function(pd_hab, prob_int, prob_euc, hnum, plotdims, titltxt) {
  pd_hab <- pd_hab[pd_hab$hab_num == 25,]
  pd <- pd_hab[,colnames(pd_hab) %in% c('eucp','intp','eucintp','unxp')]
  pd_perc <- round(pd,2)
  ### Here I need to take away the intersections from the prob_int and prob_euc too
  ei_intersect <- prob_euc*prob_int
  pd_rdm_perc <- c(eucp = prob_euc, intp = prob_int, eucintp = ei_intersect, unxp = 1-(prob_euc+prob_int-ei_intersect)) #, totp = 1
  ### Relative improvement:
  pd_relimpr <- ((pd_perc/pd_rdm_perc)-pd_rdm_perc)*100
  pd_relimpr <- as.data.frame(t(pd_relimpr))
  pd_relimpr$cat <- rownames(pd_relimpr)
  colnames(pd_relimpr) <- c('val', 'cat')
  pd_relimpr <- pd_relimpr[pd_relimpr$cat != 'unxp',]
  pd_relimpr$cat <- factor(pd_relimpr$cat, levels = c('intp','eucp','eucintp')) #'unxp',
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="intp"] <- 'Interactions'
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="eucp"] <- "Resource pref.\nsimilarity"
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="eucintp"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  ggplot(data = pd_relimpr, aes(x = cat, y = val, fill = cat)) +  
    geom_bar(stat = 'identity') + xlab('Association matches') + ylab('Deviates from\nchance [%]') +
    scale_fill_manual(values = c(int_col,euc_col,euc_col), name = 'Association\nmatches:') +
    scale_pattern_manual('stack', values = c('none','none','stripe'), guide = "none") +
    geom_col_pattern(aes(pattern = cat),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    theme_minimal() + theme(axis.text.x=element_text(angle=45)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
plt_nonrdm_neg <- function(pd_hab, prob_int, prob_euc, hnum, plotdims, titltxt) {
  pd_hab <- pd_hab[pd_hab$hab_num == 25,]
  pd <- pd_hab[,colnames(pd_hab) %in% c('eucn','intn','eucintn','unxn')]
  pd_perc <- round(pd,2)
  ### Here I need to take away the intersections from the prob_int and prob_euc too
  ei_intersect <- prob_euc*prob_int
  pd_rdm_perc <- c(eucn = prob_euc, intn = prob_int, eucintn = ei_intersect, unxn = 1-(prob_euc+prob_int-ei_intersect)) #, totp = 1
  ### Relative improvement:
  pd_relimpr <- ((pd_perc/pd_rdm_perc)-pd_rdm_perc)*100
  pd_relimpr <- as.data.frame(t(pd_relimpr))
  pd_relimpr$cat <- rownames(pd_relimpr)
  colnames(pd_relimpr) <- c('val', 'cat')
  pd_relimpr <- pd_relimpr[pd_relimpr$cat != 'unxn',]
  pd_relimpr$cat <- factor(pd_relimpr$cat, levels = c('intn','eucn','eucintn')) #'unxp',
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="intn"] <- 'Interactions'
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="eucn"] <- "Resource pref.\nsimilarity"
  levels(pd_relimpr$cat)[levels(pd_relimpr$cat)=="eucintn"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  ggplot(data = pd_relimpr, aes(x = cat, y = val, fill = cat)) +  
    geom_bar(stat = 'identity') + xlab('Association matches') + ylab('Deviates from\nchance [%]') +
    scale_fill_manual(values = c(int_col,euc_col,euc_col), name = 'Association\nmatches:') +
    scale_pattern_manual('stack', values = c('none','none','stripe'), guide = "none") +
    geom_col_pattern(aes(pattern = cat),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    theme_minimal() + theme(axis.text.x=element_text(angle=45)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}




















###########################################################################################
##### Sfig2_number_of_habitats
###########################################################################################

###########################################################################################
plot_match_eucintp <- function(pd_hab,label_mode,titltxt,plotdims) {
  rownames(pd_hab) <- pd_hab$hab_num
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA')]
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  colnames(pd) <- c('matches','n_hab','perc')
  comp <- as.data.frame(cbind(c('eucp_corA', 'intp_corA', 'eucintp_corA', 'unxpA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$n_hab,pd_hab$hab_num),]$totp
  pd$perc_n <- pd$perc*pd_hab[match(pd$n_hab,pd_hab$hab_num),]$totp
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  pd$n_hab <- factor(pd$n_hab)
  
  #colnames(pd_hab)
  pd$mnint_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_onlyintp_spsort, 2) 
  pd$mnint_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$mneuc_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_onlyintp_spsort, 2) 
  pd$mneuc_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$mnG_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_onlyintp_spsort, 2)
  pd$mnG_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$vrG_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_onlyintp_spsort, 2)
  pd$vrG_onlyintp_spsort[pd$lab != 'Int'] <- NA
  
  pd$mnint_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_onlyeucp_spsort, 2)
  pd$mnint_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$mneuc_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_onlyeucp_spsort, 2)
  pd$mneuc_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$mnG_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_onlyeucp_spsort, 2)
  pd$mnG_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$vrG_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_onlyeucp_spsort, 2)
  pd$vrG_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  
  pd$mnint_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_eucintp, 2)
  pd$mnint_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$mneuc_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_eucintp, 2)
  pd$mneuc_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$mnG_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_eucintp, 2)
  pd$mnG_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$vrG_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_eucintp, 2)
  pd$vrG_eucintp[pd$lab != 'RP & Int'] <- NA
  
  pd$mnint_noneucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_noneucintp, 2)
  pd$mnint_noneucintp[pd$lab != 'Rest'] <- NA
  pd$mneuc_noneucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_noneucintp, 2)
  pd$mneuc_noneucintp[pd$lab != 'Rest'] <- NA
  
  sumsx <- c()
  for (i in unique(pd$n_hab)) {
    sumsx <- c(sumsx,sum(pd[pd$n_hab == i,]$perc_n))
  }
  
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  pd$mneuc_onlyintp_spsort <- pd$mneuc_onlyintp_spsort*-1
  pd$mneuc_onlyeucp_spsort <- pd$mneuc_onlyeucp_spsort*-1
  pd$mneuc_eucintp <- pd$mneuc_eucintp*-1
  pd$mneuc_noneucintp <- pd$mneuc_noneucintp*-1
  
  p <- ggplot(pd, aes(y=perc_n, x=n_hab, fill = lab)) + #, col = lab
    geom_bar(position="stack", stat="identity", lwd = 0) +
    theme_test() + xlab('Number of habitats') + ylab('Positive associations') +
    scale_fill_manual(values = c(rest_col,euc_col,eucint_col,int_col), name = 'Association\nmatches:') +
    #scale_color_manual(values = c(darken(rest_col, darkfac),darken(euc_col, darkfac),darken(eucint_col, darkfac),darken(int_col, darkfac)), name = 'Association\nmatches:') +
    
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_col_pattern(aes(pattern = lab),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  
  if (label_mode == 'int'){
    p + geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_onlyintp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_onlyeucp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_eucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_noneucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col)    
  } else if (label_mode == 'euc'){
    p + geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_onlyintp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_onlyeucp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), colour = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_eucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), colour = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_noneucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = euccoef_col)
  } else {p}
}
###########################################################################################
plot_match_eucintn <- function(pd_hab,label_mode,titltxt,plotdims) {
  
  rownames(pd_hab) <- pd_hab$hab_num
  pd_haby <- pd_hab[colnames(pd_hab) %in% c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA')]
  pdy <- as.data.frame(melt(t(pd_haby)))
  pd <- c()
  for (i in unique(pdy$Var2)) {
    pdyx <- pdy[pdy$Var2 == i,] 
    pd <- rbind(pd, pdyx[c(2,3,1,4),])    
  }
  colnames(pd) <- c('matches','n_hab','perc')
  comp <- as.data.frame(cbind(c('eucn_corA', 'intn_corA', 'eucintn_corA', 'unxnA'), c('RP', 'Int', 'RP & Int', 'Rest')))
  pd$lab <- comp[match(pd$matches,comp$V1),]$V2
  pd$n_assoc <- pd_hab[match(pd$n_hab,pd_hab$hab_num),]$totn
  pd$perc_n <- pd$perc*pd_hab[match(pd$n_hab,pd_hab$hab_num),]$totn
  pd$lab <- factor(pd$lab, levels = c('Rest', 'RP', 'RP & Int', 'Int'))
  pd$n_hab <- factor(pd$n_hab)
  
  #colnames(pd_hab)
  pd$mnint_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_onlyintn_spsort, 2) 
  pd$mnint_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$mneuc_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_onlyintn_spsort, 2) 
  pd$mneuc_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$mnG_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_onlyintn_spsort, 2)
  pd$mnG_onlyintp_spsort[pd$lab != 'Int'] <- NA
  pd$vrG_onlyintp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_onlyintn_spsort, 2)
  pd$vrG_onlyintp_spsort[pd$lab != 'Int'] <- NA
  
  pd$mnint_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_onlyeucn_spsort, 2)
  pd$mnint_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$mneuc_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_onlyeucn_spsort, 2)
  pd$mneuc_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$mnG_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_onlyeucn_spsort, 2)
  pd$mnG_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  pd$vrG_onlyeucp_spsort <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_onlyeucn_spsort, 2)
  pd$vrG_onlyeucp_spsort[pd$lab != 'RP'] <- NA
  
  pd$mnint_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_eucintn, 2)
  pd$mnint_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$mneuc_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_eucintn, 2)
  pd$mneuc_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$mnG_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnG_eucintn, 2)
  pd$mnG_eucintp[pd$lab != 'RP & Int'] <- NA
  pd$vrG_eucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$vrG_eucintn, 2)
  pd$vrG_eucintp[pd$lab != 'RP & Int'] <- NA
  
  pd$mnint_noneucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mnint_noneucintn, 2)
  pd$mnint_noneucintp[pd$lab != 'Rest'] <- NA
  pd$mneuc_noneucintp <- round(pd_hab[match(pd$n_hab,pd_hab$hab_num),]$mneuc_noneucintn, 2)
  pd$mneuc_noneucintp[pd$lab != 'Rest'] <- NA
  
  sumsx <- c()
  for (i in unique(pd$n_hab)) {
    sumsx <- c(sumsx,sum(pd[pd$n_hab == i,]$perc_n))
  }
  levels(pd$lab)[levels(pd$lab)=="Rest"] <- "Rest"
  levels(pd$lab)[levels(pd$lab)=="RP"] <- "Resource pref.\nsimilarity"
  levels(pd$lab)[levels(pd$lab)=="Int"] <- "Interactions"
  levels(pd$lab)[levels(pd$lab)=="RP & Int"] <- "Resource pref.\nsimilarity and\ninteractions"
  
  pd$mneuc_onlyintp_spsort <- pd$mneuc_onlyintp_spsort*-1
  pd$mneuc_onlyeucp_spsort <- pd$mneuc_onlyeucp_spsort*-1
  pd$mneuc_eucintp <- pd$mneuc_eucintp*-1
  pd$mneuc_noneucintp <- pd$mneuc_noneucintp*-1
  
  ### Messed up in analysis function -> These are the true values 
  pd$mneuc_noneucintp <- NA#pd$mneuc_noneucintp*-1
  pd[pd$n_hab == 10 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.66
  pd[pd$n_hab == 15 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.6
  pd[pd$n_hab == 25 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.62
  pd[pd$n_hab == 50 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.61
  pd[pd$n_hab == 100 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.59
  pd[pd$n_hab == 300 & pd$matches == 'unxnA',]$mneuc_noneucintp <- -0.58
  
  pd$mnint_noneucintp <- NA#pd$mneuc_noneucintp*-1
  pd[pd$n_hab == 10 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.02
  pd[pd$n_hab == 15 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.01
  pd[pd$n_hab == 25 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.01
  pd[pd$n_hab == 50 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.01
  pd[pd$n_hab == 100 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.02
  pd[pd$n_hab == 300 & pd$matches == 'unxnA',]$mnint_noneucintp <- -0.01
  
  
  p <- ggplot(pd, aes(y=perc_n, x=n_hab, fill = lab)) + 
    geom_bar(position="stack", stat="identity", lwd = 0) +
    theme_test() + xlab('Number of habitats') + ylab('Negative associations') +
    scale_fill_manual(values = c(rest_col,euc_col,eucint_col,int_col), name = 'Association\nmatches:') +
    
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_col_pattern(aes(pattern = lab),pattern_color = NA,
                     pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    
    #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  
  if (label_mode == 'int'){
    p + geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_onlyintp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) + #'#E7D4E8FF'
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_onlyeucp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_eucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", -1*(mnint_noneucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = intcoef_col)
  } else if (label_mode == 'euc'){
    p + geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_onlyintp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), colour = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_onlyeucp_spsort))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), colour = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_eucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), colour = euccoef_col) +
      geom_label(aes(label = gsub("0\\.", "\\.", (mneuc_noneucintp))), fontface = "bold", label.size = 0.1, label.padding = unit(0.1, "lines"), size = 3, position = position_stack(vjust = 0.5), fill = alpha(c(boxcolor),boxalpha), col = euccoef_col)
  } else {p}      
}












