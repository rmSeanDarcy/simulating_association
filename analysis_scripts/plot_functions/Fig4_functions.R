###########################################################################################
##### SAMSARA - Figure 4: Cubes and prediction metrics                                #####
###########################################################################################

###########################################################################################
##### Plot co-occurrence data
###########################################################################################
coocurrence_plots <- function(pd_cub, sim_nms, colpal, symbolvals) {
  pdc <- pd_cub[pd_cub$treatment %in% sim_nms,]
  p1 <- coocpos(pdc, colpal, symbolvals, titltxt = 'Positive co-occurrences')
  p2 <- coocneg(pdc, colpal, symbolvals, titltxt = 'Negative co-occurrences')
  p3 <- precision_intp(pdc, colpal, symbolvals, titltxt = 'Positive interactions', nintp = 4)
  p4 <- precision_intn(pdc, colpal, symbolvals, titltxt = 'Negative interactions', nintn = 18)
  p5 <- precision_eucp(pdc, colpal, symbolvals, titltxt = 'Similar env.', neuc = 11)
  p6 <- precision_eucn(pdc, colpal, symbolvals, titltxt = 'Dissimilar env.', neuc = 11)
  cowplot::plot_grid(p1,p3,p5,p2,p4,p6, nrow = 2)
}
###########################################################################################
### Positive co-occurrences with scale
coocpos <- function(pdc, colpal, symbolvals, titltxt) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = totp, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + 
    scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) + ylim(0,13) +
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Positive\nco-occurrences') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Negative co-occurrences with scale
coocneg <- function(pdc, colpal, symbolvals, titltxt) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = totn, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + 
    scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Negative\nco-occurrences') +# scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Precision for environmental preference similarity with scale - positive co-occurrences
precision_eucp <- function(pdc, colpal, symbolvals, titltxt, neuc) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = eucp, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + geom_hline(yintercept = neuc/45, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('Precision') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Precision for interactions with scale - positive co-occurrences
precision_intp <- function(pdc, colpal, symbolvals, titltxt, nintp) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = intp, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + geom_hline(yintercept = nintp/45, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('Precision') +# scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Precision for environmental preference dissimilarity with scale - negative co-occurrences
precision_eucn <- function(pdc, colpal, symbolvals, titltxt, neuc) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = eucn, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + geom_hline(yintercept = neuc/45, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('Precision') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Precision for interactions with scale - negative co-occurrences
precision_intn <- function(pdc, colpal, symbolvals, titltxt, nintn) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = intn, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() + geom_hline(yintercept = nintn/45, linetype = "dashed") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('Precision') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}

###########################################################################################
##### Community data plots                                                           ######
###########################################################################################
community_plots <- function(fd_cub, sim_nms, colpal, symbolvals) {
  pdc <- pd_cub[pd_cub$treatment %in% sim_nms,]
  fdc <- fd_cub[fd_cub$treatment %in% sim_nms,]
  p1 <- cub_rch(fdc, colpal, symbolvals, titltxt = 'Richness')
  p2 <- cub_bc(fdc, colpal, symbolvals, titltxt = 'Bray-Curtis dissim.')
  p3 <- cub_numzeroassocn(pdc, colpal, symbolvals, titltxt = 'Percent zero runs with zero negative co-occurrences')
  cowplot::plot_grid(p1,p2,p3, ncol = 3)
}
###########################################################################################
### Richness with scale
cub_rch <- function(fdc, colpal, symbolvals, titltxt) {
  fdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(fdc$d, '__'))[2*(1:length(fdc$d))], '_'))[2*(1:length(fdc$d))])
  ggplot(fdc, aes(x = dim, y = mn_rch, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() +
    scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Mean species\nrichness') + # (±sd)
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Bray-Curtis dissimilarity with scale
cub_bc <- function(fdc, colpal, symbolvals, titltxt) {
  fdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(fdc$d, '__'))[2*(1:length(fdc$d))], '_'))[2*(1:length(fdc$d))])
  ggplot(fdc, aes(x = dim, y = mn_bc, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('Mean Bray-Curtis\ndissimilarity') + #(±sd)
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Percentage of zero co-occurrences per run
cub_numzeroassocn <- function(pdc, colpal, symbolvals, titltxt) {
  pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
  ggplot(pdc, aes(x = dim, y = nzeron, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
    geom_line() +
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('% runs with zero\nnegative associations') +
    scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}


###########################################################################################
###########################################################################################
##### Supporting information Figures                                                  #####
###########################################################################################
###########################################################################################

###########################################################################################
##### Co-occurrence plots for noise                                                   #####
###########################################################################################
coocurrence_plots_noise <- function(pd_cub, sim_nms, noise_lvl, colpal, symbolvals) {
  pdc <- pd_cub[pd_cub$treatment %in% sim_nms & pd_cub$simulation %in% noise_lvl,]
  p1 <- coocpos(pdc, colpal, symbolvals, titltxt = 'Positive co-occurrences')
  p2 <- coocneg(pdc, colpal, symbolvals, titltxt = 'Negative co-occurrences')
  p3 <- precision_intp(pdc, colpal, symbolvals, titltxt = 'Positive interactions', nintp = 4)
  p4 <- precision_intn(pdc, colpal, symbolvals, titltxt = 'Negative interactions', nintn = 18)
  p5 <- precision_eucp(pdc, colpal, symbolvals, titltxt = 'Similar env.', neuc = 11)
  p6 <- precision_eucn(pdc, colpal, symbolvals, titltxt = 'Dissimilar env.', neuc = 11)
  cowplot::plot_grid(p1,p3,p5,p2,p4,p6, nrow = 2)
}


###########################################################################################
##### Preference effects on correlations and  co-occurrence                           #####
###########################################################################################
##### Master function
# subdir <- "./simulation_data/fig4/rc1/rc1" #list.dirs(path='./simulation_data/fig4')
get_pref_effect <- function(subdir) {
  ### Get data on species pairwise preferences and co-occurrence
  res <- get_preferences_and_cooccurrence(subdir)
  infmres <- res$infmres
  res_pos <- res$res_pos
  res_neg <- res$res_neg
  # Reorganise data
  res <- prep_preferences_and_cooccurrence_data(infmres)
  pres <- res$pres 
  nres <- res$nres
  # Prepare for plotting
  res <- get_pairwise_preference_correlation_data(pres, nres)
  pdf <- res$pdf
  ndf <- res$ndf
  # Plot Spearmans correlation for same resource preference values (gamma)
  plt_corcoeffs_pos_rs3_samepref(pdf)
}

###########################################################################################
##### Pairwise preferences and co-occurrence
###########################################################################################
### Get data on species pairwise preferences and co-occurrence
get_preferences_and_cooccurrence <- function(subdir) {
  infmres <- c()
  res_pos <- c()
  res_neg <- c()
  ###########################################################################################
  ### Add correlation coefficients to every pair
  sp_correl_log <- read.csv(paste0(subdir,'/sp_correl_log.csv'))[,-1]
  adjm <- c()
  #d <- unique(sp_correl_log$d)[1]
  for (d in unique(sp_correl_log$d)) {
    sp_correld <- sp_correl_log[sp_correl_log$d == d,]
    
    #r <- unique(sp_correl_log$repl)[1]
    for (r in unique(sp_correld$repl)) {
      adjmr <- sp_correld[sp_correld$repl == r,1:10]
      adjmr[is.na(adjmr)] <- 0
      diag(adjmr) <- NA
      adjmr[lower.tri(adjmr)] <- NA
      rownames(adjmr) <- colnames(adjmr)
      adjmr <- reshape2::melt(as.matrix(adjmr))
      adjmr <- adjmr[!is.na(adjmr$value),]
      adjm <- rbind(adjm, cbind(adjmr,d = d, repl = r))
    }
  }
  ###########################################################################################
  ### Load community data
  npop_cubd_log <- read.csv(paste0(subdir,'/npop_cubd_log.csv'))[,-1]
  kabs_cubd_log <- read.csv(paste0(subdir,'/kabs_cubd_log.csv'))[,-1]
  ##########################################################################################
  ### Main datasets on associations and species resource preferences
  infmat_res_comp <- read.csv(paste0(subdir,'/infmat_res_comp.csv'))[,-1]
  infm <- infmat_res_comp 
  infm$corcoeff <- adjm$value 
  ###########################################################################################
  ### Adding species preferences to association data
  sp_rsc_prf <- read.csv(paste0(subdir,'/sp_rsc_prf_log.csv'))[,-1]
  # Simplistic info on whether species are specialists (defined as preference for one resource >0.5)
  sp_rsc_prf$spgen <- NA
  sp_rsc_prf[sp_rsc_prf$R1 > 0.5,]$spgen <- 'R1'
  sp_rsc_prf[sp_rsc_prf$R2 > 0.5,]$spgen <- 'R2'
  sp_rsc_prf[sp_rsc_prf$R3 > 0.5,]$spgen <- 'R3'
  sp_rsc_prf[sp_rsc_prf$R1 < 0.5 & sp_rsc_prf$R2 < 0.5 & sp_rsc_prf$R3 < 0.5,]$spgen <- 'G'
  sp_rsc_prf$GR1 <- sqrt((sp_rsc_prf$R1 - 1)^2 + 
                           (sp_rsc_prf$R2 - 0)^2 + 
                           (sp_rsc_prf$R3 - 0)^2)
  sp_rsc_prf$GR2 <- sqrt((sp_rsc_prf$R1 - 0)^2 + 
                           (sp_rsc_prf$R2 - 1)^2 + 
                           (sp_rsc_prf$R3 - 0)^2)
  sp_rsc_prf$GR3 <- sqrt((sp_rsc_prf$R1 - 0)^2 + 
                           (sp_rsc_prf$R2 - 0)^2 + 
                           (sp_rsc_prf$R3 - 1)^2)
  sp_rsc_prf$GG <- sqrt((sp_rsc_prf$R1 - 0.333333)^2 + 
                          (sp_rsc_prf$R2 - 0.333333)^2 + 
                          (sp_rsc_prf$R3 - 0.333333)^2)
  
  #d <- unique(infm$d)[9]
  
  for (d in unique(sp_correl_log$d)) {
    infmd <- infm[infm$d == d,]
    res_posr <- c()
    res_negr <- c()
    #r <- 1
    for (r in unique(infm$repl)) {
      infmr <- infmd[infmd$repl == r,]
      srprfr <- sp_rsc_prf[sp_rsc_prf$repl == r,]
      infmr$R1i = NA
      infmr$R1j = NA
      infmr$R2i = NA
      infmr$R2j = NA
      infmr$R3i = NA
      infmr$R3j = NA
      #i <- 0
      for (i in 1:nrow(infmr)) {
        infmr[i,]$R1i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R1
        infmr[i,]$R1j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R1
        infmr[i,]$R2i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R2
        infmr[i,]$R2j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R2
        infmr[i,]$R3i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R3
        infmr[i,]$R3j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R3
      }
      
      infmr$GR1i = NA
      infmr$GR1j = NA
      infmr$GR2i = NA
      infmr$GR2j = NA
      infmr$GR3i = NA
      infmr$GR3j = NA
      for (i in 1:nrow(infmr)) {
        infmr[i,]$GR1i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR1
        infmr[i,]$GR1j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR1
        infmr[i,]$GR2i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR2
        infmr[i,]$GR2j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR2
        infmr[i,]$GR3i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR3
        infmr[i,]$GR3j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR3
      }
      
      infmres <- rbind(infmres, cbind(runm = subdir, infmr))
      
      ###########################################################################################
      infmrp <- infmr[infmr$cooc > 0,]
      R1i_R1jp <- cor(infmrp$R1i, infmrp$R1j, method = 'spearman')
      R2i_R2jp <- cor(infmrp$R2i, infmrp$R2j, method = 'spearman')
      R3i_R3jp <- cor(infmrp$R3i, infmrp$R3j, method = 'spearman')
      
      R1i_R2jp <- cor(infmrp$R1i, infmrp$R2j, method = 'spearman')
      R2i_R1jp <- cor(infmrp$R2i, infmrp$R1j, method = 'spearman')
      
      R1i_R3jp <- cor(infmrp$R1i, infmrp$R3j, method = 'spearman')
      R3i_R1jp <- cor(infmrp$R3i, infmrp$R1j, method = 'spearman')
      
      R2i_R3jp <- cor(infmrp$R2i, infmrp$R3j, method = 'spearman')
      R3i_R2jp <- cor(infmrp$R3i, infmrp$R2j, method = 'spearman')
      
      res_posr <- rbind(res_posr, c(cR1 = R1i_R1jp, cR2 = R2i_R2jp, cR3 = R3i_R3jp,
                                    cR12 = (R1i_R2jp+R2i_R1jp)/2, cR13 = (R1i_R3jp+R3i_R1jp)/2, cR23 = (R2i_R3jp+R3i_R2jp)/2))
      
      ###########################################################################################
      infmrn <- infmr[infmr$cooc < 0,]
      R1i_R1jn <- cor(infmrn$R1i, infmrn$R1j, method = 'spearman')
      R2i_R2jn <- cor(infmrn$R2i, infmrn$R2j, method = 'spearman')
      R3i_R3jn <- cor(infmrn$R3i, infmrn$R3j, method = 'spearman')
      
      R1i_R2jn <- cor(infmrn$R1i, infmrn$R2j, method = 'spearman')
      R2i_R1jn <- cor(infmrn$R2i, infmrn$R1j, method = 'spearman')
      
      R1i_R3jn <- cor(infmrn$R1i, infmrn$R3j, method = 'spearman')
      R3i_R1jn <- cor(infmrn$R3i, infmrn$R1j, method = 'spearman')
      
      R2i_R3jn <- cor(infmrn$R2i, infmrn$R3j, method = 'spearman')
      R3i_R2jn <- cor(infmrn$R3i, infmrn$R2j, method = 'spearman')
      
      res_negr <- rbind(res_negr, c(cR1 = R1i_R1jn, cR2 = R2i_R2jn, cR3 = R3i_R3jn,
                                    cR12 = (R1i_R2jn+R2i_R1jn)/2, cR13 = (R1i_R3jn+R3i_R1jn)/2, cR23 = (R2i_R3jn+R3i_R2jn)/2))
    }
    res_pos <- rbind(res_pos, c(runm = subdir, d = d, colMeans(res_posr, na.rm = TRUE)))
    res_neg <- rbind(res_neg, c(runm = subdir, d = d, colMeans(res_negr, na.rm = TRUE)))
    print(paste('Done for',d,subdir))
  }
  return(list(infmres = infmres, res_pos = res_pos, res_neg = res_neg))
}
###########################################################################################
### Prepare for plotting
prep_preferences_and_cooccurrence_data <- function(infmres) {
  res <- as.data.frame(infmres)
  res$dim <- as.numeric(str_split_fixed(res$d, pattern = '_', 5)[,5])
  res[,! colnames(res) %in% c('runm','N1','N2','d')] <- sapply(res[,! colnames(res) %in% c('runm','N1','N2','d')],as.numeric)
  res$intbin <- 'Other'
  res[res$int>0,]$intbin <- 'Competition'
  res[res$int<0,]$intbin <- 'Mutualism'
  pres <- res[res$corcoeff > 0,]
  nres <- res[res$corcoeff < 0,]
  return(list(pres = pres, nres = nres))  
}
###########################################################################################
prep_d <- function(data) {
  dat <- as.data.frame(data)
  dat$dim <- as.numeric(str_split_fixed(dat$d, pattern = '_', 5)[,5])
  dat$rscdist <- str_split_fixed(dat$runm, pattern = '_', 5)[,4]
  dat[,4:10] <- sapply(dat[,4:10],as.numeric)
  return(dat)
}
###########################################################################################
### Collect data
get_pairwise_preference_correlation_data <- function(pres, nres) {
  res_p_all <- c()
  res_p_mut <- c()
  res_p_comp <- c()
  res_p_nonint <- c()
  
  res_n_all <- c()
  res_n_comp <- c()
  res_n_nonint <- c()
  #rnm <- unique(pres$runm)[1]
  for (rnm in unique(pres$runm)) {
    #d <- unique(pres$d)[1]
    for (d in unique(pres$d)) {
      ##############################################################
      ### Positive
      # All pairs
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7,]
      resp <- c(R1_R1 = cor(presx$R1i, presx$R1j, method = 'spearman'),
                R2_R2 = cor(presx$R2i, presx$R2j, method = 'spearman'),
                R3_R3 = cor(presx$R3i, presx$R3j, method = 'spearman'),
                R1_R2 = (cor(presx$R1i, presx$R2j, method = 'spearman') + cor(presx$R2i, presx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(presx$R1i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(presx$R2i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R2j, method = 'spearman'))/2)
      res_p_all <- rbind(res_p_all, c(runm = rnm, d = d, type = 'All', resp))
      # Only mutualistic pairs
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7 & pres$int < 0,]
      resp <- c(R1_R1 = cor(presx$R1i, presx$R1j, method = 'spearman'),
                R2_R2 = cor(presx$R2i, presx$R2j, method = 'spearman'),
                R3_R3 = cor(presx$R3i, presx$R3j, method = 'spearman'),
                R1_R2 = (cor(presx$R1i, presx$R2j, method = 'spearman') + cor(presx$R2i, presx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(presx$R1i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(presx$R2i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R2j, method = 'spearman'))/2)
      res_p_mut <- rbind(res_p_mut, c(runm = rnm, d = d, type = 'Mutualism', resp))
      
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7 & pres$int > 0,]
      resp <- c(R1_R1 = cor(presx$R1i, presx$R1j, method = 'spearman'),
                R2_R2 = cor(presx$R2i, presx$R2j, method = 'spearman'),
                R3_R3 = cor(presx$R3i, presx$R3j, method = 'spearman'),
                R1_R2 = (cor(presx$R1i, presx$R2j, method = 'spearman') + cor(presx$R2i, presx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(presx$R1i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(presx$R2i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R2j, method = 'spearman'))/2)
      res_p_comp <- rbind(res_p_comp, c(runm = rnm, d = d, type = 'Competition', resp))
      
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7 & pres$int == 0,]
      resp <- c(R1_R1 = cor(presx$R1i, presx$R1j, method = 'spearman'),
                R2_R2 = cor(presx$R2i, presx$R2j, method = 'spearman'),
                R3_R3 = cor(presx$R3i, presx$R3j, method = 'spearman'),
                R1_R2 = (cor(presx$R1i, presx$R2j, method = 'spearman') + cor(presx$R2i, presx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(presx$R1i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(presx$R2i, presx$R3j, method = 'spearman') + cor(presx$R3i, presx$R2j, method = 'spearman'))/2)
      res_p_nonint <- rbind(res_p_nonint, c(runm = rnm, d = d, type = 'No interaction', resp))
      
      ##############################################################
      ### Positive
      # All pairs
      nresx <- nres[nres$runm == rnm & nres$d == d & nres$corcoeff < -0.7,]
      resn <- c(R1_R1 = cor(nresx$R1i, nresx$R1j, method = 'spearman'),
                R2_R2 = cor(nresx$R2i, nresx$R2j, method = 'spearman'),
                R3_R3 = cor(nresx$R3i, nresx$R3j, method = 'spearman'),
                R1_R2 = (cor(nresx$R1i, nresx$R2j, method = 'spearman') + cor(nresx$R2i, nresx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(nresx$R1i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(nresx$R2i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R2j, method = 'spearman'))/2)
      res_n_all <- rbind(res_n_all, c(runm = rnm, d = d, type = 'All', resn))
      # Only mutualistic pairs
      nresx <- nres[nres$runm == rnm & nres$d == d & nres$corcoeff < -0.7 & nres$int > 0,]
      resn <- c(R1_R1 = cor(nresx$R1i, nresx$R1j, method = 'spearman'),
                R2_R2 = cor(nresx$R2i, nresx$R2j, method = 'spearman'),
                R3_R3 = cor(nresx$R3i, nresx$R3j, method = 'spearman'),
                R1_R2 = (cor(nresx$R1i, nresx$R2j, method = 'spearman') + cor(nresx$R2i, nresx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(nresx$R1i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(nresx$R2i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R2j, method = 'spearman'))/2)
      res_n_comp <- rbind(res_n_comp, c(runm = rnm, d = d, type = 'Competition', resn))
      
      nresx <- nres[nres$runm == rnm & nres$d == d & nres$corcoeff < -0.7 & nres$int == 0,]
      resn <- c(R1_R1 = cor(nresx$R1i, nresx$R1j, method = 'spearman'),
                R2_R2 = cor(nresx$R2i, nresx$R2j, method = 'spearman'),
                R3_R3 = cor(nresx$R3i, nresx$R3j, method = 'spearman'),
                R1_R2 = (cor(nresx$R1i, nresx$R2j, method = 'spearman') + cor(nresx$R2i, nresx$R1j, method = 'spearman'))/2,
                R1_R3 = (cor(nresx$R1i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R1j, method = 'spearman'))/2,
                R2_R3 = (cor(nresx$R2i, nresx$R3j, method = 'spearman') + cor(nresx$R3i, nresx$R2j, method = 'spearman'))/2)
      res_n_nonint <- rbind(res_n_nonint, c(runm = rnm, d = d, type = 'No interaction', resn))
    }
  }
  pdf <- rbind(prep_d(res_p_all),prep_d(res_p_mut),prep_d(res_p_comp),prep_d(res_p_nonint))
  ndf <- rbind(prep_d(res_n_all),prep_d(res_n_comp),prep_d(res_n_nonint))
  return(list(pdf = pdf, ndf = ndf))  
}
###########################################################################################
### Plot preferences
plt_corcoeffs_pos_rs3_samepref <- function(pdf) {
  pdfx <- pdf
  pdfx <- pdfx[pdfx$runm == subdir,]
  R1_R1 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R1_R1')], 'R1 x R1')
  R2_R2 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R2_R2')], 'R2 x R2')
  R3_R3 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R3_R3')], 'R3 x R3')
  colnames(R1_R1) <- c('runm','type','comp','d','spec')
  colnames(R2_R2) <- c('runm','type','comp','d','spec')
  colnames(R3_R3) <- c('runm','type','comp','d','spec')
  df <- rbind(R1_R1,R2_R2,R3_R3)
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(df, aes(x = d, y = comp, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + geom_line() + #stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ spec) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}

