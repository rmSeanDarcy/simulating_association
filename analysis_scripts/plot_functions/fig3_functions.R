###########################################################################################
### SAMSARA - Figure 3 and plots for related Supporting information                     ###
###########################################################################################

###########################################################################################
##### Fig.3: Main manuscript                                                          #####
###########################################################################################

###########################################################################################
##### Fig.3 A, B) Species richness and dissimilarity                                  #####
###########################################################################################
# Richness
richness_radius <- function(fdh) {
  fdh$radius <- as.numeric(str_split_i(fdh$simulation, 'radius_', -1))
  ggplot(fdh, aes(x = radius, y = mn_rch, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    stat_smooth(se = FALSE) + #geom_ribbon(aes(ymin = mn_rch-sd_rch, ymax = mn_rch+sd_rch), alpha = 0.1, colour = NA) +  
    geom_point(size = 1, stroke = 1) + theme_classic() + geom_vline(xintercept = 0.22, linetype = "dashed") + 
    xlab(TeX('Connectivity radius')) + ylab('Mean species\nrichness') +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = round(seq(min(fdh$radius), max(fdh$radius), by = 0.25),2)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Bray-Curtis dissimilarity (communities)
dissimilarity_radius <- function(fdh) {
  fdh$radius <- as.numeric(str_split_i(fdh$simulation, 'radius_', -1))
  ggplot(fdh, aes(x = radius, y = mn_bc, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    stat_smooth(se = FALSE) + #geom_ribbon(aes(ymin = mn_rch-sd_rch, ymax = mn_rch+sd_rch), alpha = 0.1, colour = NA) +  
    geom_point(size = 1, stroke = 1) + theme_classic() + geom_vline(xintercept = 0.22, linetype = "dashed") + 
    xlab(TeX('Connectivity radius')) + ylab('Mean Bray-Curtis\ndissimilarity') +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = round(seq(min(fdh$radius), max(fdh$radius), by = 0.25),2)) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}

###########################################################################################
##### Fig.3 C, D) Species richness and dissimilarity                                  #####
###########################################################################################
# Positive associations
cooc_radius_pos <- function(pdh,titltxt) {
  pdh$radius <- as.numeric(str_split_i(pdh$simulation, 'radius_', -1))
  pd <- as.data.frame(rbind(cbind(pdh$radius,pdh$treatment,pdh$totp,pdh$eucp_corA,'eucp_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totp,pdh$intp_corA,'intp_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totp,pdh$eucintp_corA,'eucintp_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totp,pdh$unxpA,'unxpA')))
  colnames(pd) <- c('radius','treatment','totp','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$perc_n <- pd$perc*as.numeric(pd$totp)
  sumsx <- c()
  pd$treatment <- factor(pd$treatment)
  pd$matches <- factor(pd$matches, levels = c('unxpA','eucp_corA','eucintp_corA','intp_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) +  #col = matches, 
    geom_area(position="stack", stat="identity") +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Positive\nco-occurrences') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Co-occurrence\nmatches:', labels = colabels) +
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_area_pattern(aes(pattern = matches),pattern_color = NA,
                      pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    scale_x_continuous(breaks = round(seq(min(pd$radius), max(pd$radius), by = 0.25),2)) + geom_vline(xintercept = 0.22, linetype = "dashed") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative associations
cooc_radius_neg <- function(pdh,titltxt) {
  pdh$radius <- as.numeric(str_split_i(pdh$simulation, 'radius_', -1))
  pd <- as.data.frame(rbind(cbind(pdh$radius,pdh$treatment,pdh$totn,pdh$eucn_corA,'eucn_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totn,pdh$intn_corA,'intn_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totn,pdh$eucintn_corA,'eucintn_corA'),
                            cbind(pdh$radius,pdh$treatment,pdh$totn,pdh$unxnA,'unxnA')))
  colnames(pd) <- c('radius','treatment','totp','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$perc_n <- pd$perc*as.numeric(pd$totp)
  pd$treatment <- factor(pd$treatment)
  pd$matches <- factor(pd$matches, levels = c('unxnA','eucn_corA','eucintn_corA','intn_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) +  #col = matches, 
    geom_area(position="stack", stat="identity") +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Positive\nco-occurrences') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Co-occurrence\nmatches:', labels = colabels) +
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
### Positive co-occurrences
cooc_radius_noise_pos <- function(pdh) {
  pdh$radius <- as.numeric(unlist(strsplit(unlist(strsplit(pdh$simulation, '__'))[(2*(1:length(pdh$simulation)))-1], '_'))[2*(1:length(pdh$simulation))])
  pdh$noise <- as.numeric(unlist(strsplit(unlist(strsplit(pdh$simulation, '__'))[2*(1:length(pdh$simulation))], '_'))[2*(1:length(pdh$simulation))])
  pd <- c()
  for (i in 1:nrow(pdh)) {
    pdx <- as.data.frame(rbind(cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totp,pdh[i,]$noise,pdh[i,]$eucp_corA,'eucp_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totp,pdh[i,]$noise,pdh[i,]$intp_corA,'intp_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totp,pdh[i,]$noise,pdh[i,]$eucintp_corA,'eucintp_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totp,pdh[i,]$noise,pdh[i,]$unxpA,'unxpA')))
    pd <- rbind(pd,pdx)
  }
  colnames(pd) <- c('radius','treatment','totp','noise','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$perc_n <- pd$perc*as.numeric(pd$totp)
  pd$treatment <- factor(pd$treatment)
  pd$matches <- factor(pd$matches, levels = c('unxpA','eucp_corA','eucintp_corA','intp_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) +  #col = matches, 
    geom_area(position="stack", stat="identity") + facet_grid(~ noise) +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Positive\nco-occurrences') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Co-occurrence\nmatches:', labels = colabels) +
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_area_pattern(aes(pattern = matches),pattern_color = NA,
                      pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    scale_x_continuous(breaks = round(seq(min(pd$radius), max(pd$radius), by = 0.25),2)) + geom_vline(xintercept = 0.22, linetype = "dashed") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control
}
### Negative co-occurrences
cooc_radius_noise_neg <- function(pdht) {
  pdh$radius <- as.numeric(unlist(strsplit(unlist(strsplit(pdh$simulation, '__'))[(2*(1:length(pdh$simulation)))-1], '_'))[2*(1:length(pdh$simulation))])
  pdh$noise <- as.numeric(unlist(strsplit(unlist(strsplit(pdh$simulation, '__'))[2*(1:length(pdh$simulation))], '_'))[2*(1:length(pdh$simulation))])
  pd <- c()
  for (i in 1:nrow(pdh)) {
    pdx <- as.data.frame(rbind(cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totn,pdh[i,]$noise,pdh[i,]$eucn_corA,'eucn_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totn,pdh[i,]$noise,pdh[i,]$intn_corA,'intn_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totn,pdh[i,]$noise,pdh[i,]$eucintn_corA,'eucintn_corA'),
                               cbind(pdh[i,]$radius,pdh[i,]$treatment,pdh[i,]$totn,pdh[i,]$noise,pdh[i,]$unxnA,'unxnA')))
    pd <- rbind(pd,pdx)
  }
  colnames(pd) <- c('radius','treatment','totn','noise','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$perc_n <- pd$perc*as.numeric(pd$totn)
  pd$treatment <- factor(pd$treatment)
  pd$noise <- factor(pd$noise)
  pd$matches <- factor(pd$matches, levels = c('unxnA','eucn_corA','eucintn_corA','intn_corA'))
  colabels <- c('Rest','Similar pref.','Similar pref. &\ninteractions','Interactions')
  pd$radius <- as.numeric(pd$radius)
  pd[is.na(pd)] <- 0
  ggplot(pd, aes(y=perc_n, x=radius, fill = matches, goup = matches)) +  #col = matches, 
    geom_area(position="stack", stat="identity") + facet_grid(~ noise) +
    theme_classic() + xlab(TeX('Connectivity radius $\\d_{e}$')) + ylab('Negative\nco-occurrences') +
    scale_fill_manual(values = c(rest_col,euc_col,euc_col,int_col)
                      , name = 'Co-occurrence\nmatches:', labels = colabels) +
    scale_pattern_manual('stack', values = c('none','none','stripe','none'), guide = "none") +
    geom_area_pattern(aes(pattern = matches),pattern_color = NA,
                      pattern_size = 1,pattern_density = .5,pattern_spacing = .075,pattern_fill = int_col,color = NA) +
    guides(fill = guide_legend(override.aes = list(pattern = c("none", "none", "stripe", "none")))) + #guides(fill=guide_legend(override.aes=list(pattern="none"))) +
    scale_x_continuous(breaks = round(seq(min(pd$radius), max(pd$radius), by = 0.25),2)) + geom_vline(xintercept = 0.22, linetype = "dashed") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control
}


