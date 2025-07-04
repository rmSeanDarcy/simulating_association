###########################################################################################
### SAMSARA - Figure 3 and plots for related Supporting information                     ###
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
##### Fig.3: Main manuscript                                                          #####
###########################################################################################
###########################################################################################

###########################################################################################
##### Fig.3 A, B) Species richness and dissimilarity                                  #####
###########################################################################################
# Richness
richness_radius <- function(pdc, plotdims, titltxt) {
  fd_habx <- fd_hab
  fd_habx$run_nm <- as.character(fd_habx$run_nm)
  fd_habx$radius <- str_split_i(fd_habx$run_nm, '__Xradius', -1)#, "[[", 2)
  fd_habx$radius <- as.numeric(str_split_i(fd_habx$radius, '__Z', 1))#, "[[", 2)
  fdhx <- fd_habx
  fdhx$rscdist <- factor(fdhx$rscdist, levels = dist_nms)
  
  ggplot(fdhx, aes(x = radius, y = mn_rch, col = rscdist, group = rscdist, fill = rscdist, shape = rscdist)) +
    stat_smooth(se = FALSE) + #geom_ribbon(aes(ymin = mn_rch-sd_rch, ymax = mn_rch+sd_rch), alpha = 0.1, colour = NA) +  
    geom_point(size = 1, stroke = 1) + theme_classic() + geom_vline(xintercept = 0.22, linetype = "dashed") + 
    xlab(TeX('Connectivity radius')) + ylab('Mean species\nrichness') +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = round(seq(min(fdhx$radius), max(fdhx$radius), by = 0.25),2)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Bray-Curtis dissimilarity (communities)
dissimilarity_radius <- function(fd_hab, plotdims, titltxt) {
  fd_habx <- fd_hab
  fd_habx$run_nm <- as.character(fd_habx$run_nm)
  fd_habx$radius <- str_split_i(fd_habx$run_nm, '__Xradius', -1)#, "[[", 2)
  fd_habx$radius <- as.numeric(str_split_i(fd_habx$radius, '__Z', 1))#, "[[", 2)
  fdhx <- fd_habx
  fdhx$rscdist <- factor(fdhx$rscdist, levels = dist_nms)
  
  ggplot(fdhx, aes(x = radius, y = mn_bc, col = rscdist, group = rscdist, fill = rscdist, shape = rscdist)) +
    stat_smooth(se = FALSE) + #geom_ribbon(aes(ymin = mn_rch-sd_rch, ymax = mn_rch+sd_rch), alpha = 0.1, colour = NA) +  
    geom_point(size = 1, stroke = 1) + theme_classic() + geom_vline(xintercept = 0.22, linetype = "dashed") + 
    xlab(TeX('Connectivity radius')) + ylab('Mean species\nrichness') +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = round(seq(min(fdhx$radius), max(fdhx$radius), by = 0.25),2)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}

###########################################################################################
##### Fig.3 C, D) Species richness and dissimilarity                                  #####
###########################################################################################
# Positive associations
cooc_radius_pos <- function(pd_hab, nhab, rscdistrib, plotdims, titltxt) {
  pd_hab$radius <- as.numeric(sapply(str_split(sapply(str_split(pd_hab$run_nm, '__Xradius'), "[[", 2), '__Z'), "[[", 1))
  pd_hab <- pd_hab[pd_hab$rscdist == rscdistrib,]
  pd <- as.data.frame(rbind(cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$eucp_corA,'eucp_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$intp_corA,'intp_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$eucintp_corA,'eucintp_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$unxpA,'unxpA')))
  colnames(pd) <- c('radius','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$n_assoc <- pd_hab[match(pd$radius,pd_hab$radius),]$totp
  pd$perc_n <- pd$perc*pd_hab[match(pd$radius,pd_hab$radius),]$totp
  sumsx <- c()
  for (i in unique(pd$radius)) {
    sumsx <- c(sumsx,sum(pd[pd$radius == i,]$perc_n))
  }
  pd$rscdist <- factor(pd$rscdist, levels = c('same0','same5','same10','diffL'))
  xlabels <- c(TeX('$\\Lambda_{0}$'),TeX('$\\Lambda_{5}$'),TeX('$\\Lambda_{10}$'), TeX('$\\Lambda_{mix}$'))
  pd$matches <- factor(pd$matches, levels = c('unxpA','eucp_corA','eucintp_corA','intp_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) +  #col = matches, 
    geom_area(position="stack", stat="identity") +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Positive\nassociations') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Association\nmatches:', labels = colabels) +
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_area_pattern(aes(pattern = matches),pattern_color = NA,
                      pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    scale_x_continuous(breaks = round(seq(min(pd$radius), max(pd$radius), by = 0.25),2)) + geom_vline(xintercept = 0.22, linetype = "dashed") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative associations
cooc_radius_neg <- function(pd_hab, nhab, rscdistrib, plotdims, titltxt) {
  pd_hab$radius <- as.numeric(sapply(str_split(sapply(str_split(pd_hab$run_nm, '__Xradius'), "[[", 2), '__Z'), "[[", 1))
  pd_hab <- pd_hab[pd_hab$rscdist == rscdistrib,]
  pd <- as.data.frame(rbind(cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$eucn_corA,'eucn_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$intn_corA,'intn_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$eucintn_corA,'eucintn_corA'),
                            cbind(pd_hab$radius,pd_hab$rscdist,pd_hab$unxnA,'unxnA')))
  colnames(pd) <- c('radius','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$n_assoc <- pd_hab[match(pd$radius,pd_hab$radius),]$totn
  pd$perc_n <- pd$perc*pd_hab[match(pd$radius,pd_hab$radius),]$totn
  sumsx <- c()
  for (i in unique(pd$radius)) {
    sumsx <- c(sumsx,sum(pd[pd$radius == i,]$perc_n))
  }
  pd$rscdist <- factor(pd$rscdist, levels = c('same0','same5','same10','diffL'))
  xlabels <- c(TeX('$\\Lambda_{0}$'),TeX('$\\Lambda_{5}$'),TeX('$\\Lambda_{10}$'), TeX('$\\Lambda_{mix}$'))
  pd$matches <- factor(pd$matches, levels = c('unxnA','eucn_corA','eucintn_corA','intn_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) + 
    geom_area(position="stack", stat="identity") +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Negative\nassociations') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Association\nmatches:', labels = colabels) +
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_area_pattern(aes(pattern = matches),pattern_color = NA,
                      pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    scale_x_continuous(breaks = round(seq(min(pd$radius), max(pd$radius), by = 0.25),2)) + geom_vline(xintercept = 0.22, linetype = "dashed") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}


###########################################################################################
###########################################################################################
##### Fig.3: Supporting information                                                   #####
###########################################################################################
###########################################################################################

###########################################################################################
##### SFig3_noise                                                                     #####
###########################################################################################


























