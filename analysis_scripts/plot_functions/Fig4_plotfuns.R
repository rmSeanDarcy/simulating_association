###########################################################################################
##### SAMSARA - Figure 4: Cubes and prediction metrics                                #####
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
##### Plots for Figure 4 main manuscript                                              #####
###########################################################################################
###########################################################################################
# Most plots from this section are used in the main manuscript, some are added to Figures in the Supporting information

###########################################################################################
##### Plot all co-occurrence data                                                     #####
###########################################################################################

###########################################################################################
### Master plot for better expoting
coocurrence_plots <- function(run_nms,dist_nms, plotdims, hab_num) {
  labs <- as.data.frame(cbind(run_nms,dist_nms))
  pdc <- pd_cub1[pd_cub1$run_nm %in% labs$run_nms,]
  fdc <- fd_cub1[fd_cub1$run_nm %in% labs$run_nms,]
  p1 <- coocpos(pdc, labs, plotdims, titltxt = 'Positive associations', hab_num)
  p2 <- coocneg(pdc, labs, plotdims, titltxt = 'Negative associations', hab_num)
  p3 <- precision_intp(pdc, labs, plotdims, nintp = 4, titltxt = 'Precision pos. interactions', hab_num)
  p4 <- precision_intn(pdc, labs, plotdims, nintn = 18, titltxt = 'Precision neg. interactions', hab_num)
  p5 <- precision_eucp(pdc, labs, neuc = 11,plotdims, titltxt = 'Precision preference sim.', hab_num)
  p6 <- precision_eucn(pdc, labs, neuc = 11, plotdims, titltxt = 'Precision preference dissim.', hab_num)
  cowplot::plot_grid(p1,p3,p5,p2,p4,p6, nrow = 2)
}
###########################################################################################
### Positive co-occurrences with scale
coocpos <- function(pdc, labs, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$mean_habnum <- as.numeric(as.character(pdcx$mean_habnum))
  pdcx$labs <- factor(pdcx$labs, levels = dist_nms)
  if(hab_num == TRUE) {
    ggplot(pdcx, aes(x = mean_habnum, y = totp, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Number of habitats')) + ylab('Positive\nassociations') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)    
  } else if (hab_num == FALSE) {
    ggplot(pdcx, aes(x = dims, y = totp, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Composite sample side length l')) + ylab('Positive\nassociations') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Negative co-occurrences with scale
coocneg <- function(pdc, labs, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$mean_habnum <- as.numeric(as.character(pdcx$mean_habnum))
  pdcx$labs <- factor(pdcx$labs, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pdcx, aes(x = mean_habnum, y = totn, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Number of habitats')) + ylab('Negative\nassociations') +# scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pdcx, aes(x = dims, y = totn, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Composite sample side length l')) + ylab('Negative\nassociations') +# scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
##### Precision plots                                                                 #####
###########################################################################################
### Precision for environmental preference similarity with scale - positive co-occurrences
precision_eucp <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pd <- as.data.frame(cbind(as.character(pdcx$mean_habnum),as.numeric(as.character(pdcx$dim)),pdcx$labs,pdcx$eucp,'eucp'))
  colnames(pd) <- c('mean_habnum','dim','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$mean_habnum <- as.numeric(pd$mean_habnum)
  pd$dim <- as.numeric(pd$dim)
  pd$rscdist <- factor(pd$rscdist, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pd, aes(x = mean_habnum, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.25, linetype = "dashed") +
      xlab(TeX('Number of habitats')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pd, aes(x = dim, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.25, linetype = "dashed") +
      xlab(TeX('Composite sample side length l')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Precision for interactions with scale - positive co-occurrences
precision_intp <- function(pdc, labs, nintp, plotdims, titltxt, hab_num) {
  
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pd <- as.data.frame(cbind(as.character(pdcx$mean_habnum),as.numeric(as.character(pdcx$dim)),pdcx$labs,pdcx$intp,'intp'))
  colnames(pd) <- c('mean_habnum','dim','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$mean_habnum <- as.numeric(pd$mean_habnum)
  pd$dim <- as.numeric(pd$dim)
  pd$rscdist <- factor(pd$rscdist, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pd, aes(x = mean_habnum, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.1, linetype = "dashed") +
      xlab(TeX('Number of habitats')) + ylab('Precision') +# scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pd, aes(x = dim, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.1, linetype = "dashed") +
      xlab(TeX('Composite sample side length l')) + ylab('Precision') +# scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Precision for environmental preference dissimilarity with scale - negative co-occurrences
precision_eucn <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pd <- as.data.frame(cbind(as.character(pdcx$mean_habnum),as.numeric(as.character(pdcx$dim)),pdcx$labs,pdcx$eucn,'eucn'))
  colnames(pd) <- c('mean_habnum','dim','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$mean_habnum <- as.numeric(pd$mean_habnum)
  pd$dim <- as.numeric(pd$dim)
  pd$rscdist <- factor(pd$rscdist, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pd, aes(x = mean_habnum, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.25, linetype = "dashed") +
      xlab(TeX('Number of habitats')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pd, aes(x = dim, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.25, linetype = "dashed") +
      xlab(TeX('Composite sample side length l')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Precision for interactions with scale - negative co-occurrences
precision_intn <- function(pdc, labs, nintn, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pd <- as.data.frame(cbind(as.character(pdcx$mean_habnum),as.numeric(as.character(pdcx$dim)),pdcx$labs,pdcx$intn,'intn'))
  colnames(pd) <- c('mean_habnum','dim','rscdist','perc','matches')
  pd$perc <- as.numeric(pd$perc) 
  pd$mean_habnum <- as.numeric(pd$mean_habnum)
  pd$dim <- as.numeric(pd$dim)
  pd$rscdist <- factor(pd$rscdist, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pd, aes(x = mean_habnum, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.4, linetype = "dashed") +
      xlab(TeX('Number of habitats')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pd, aes(x = dim, y = perc, col = rscdist, group = rscdist, shape = rscdist)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + geom_hline(yintercept = 0.4, linetype = "dashed") +
      xlab(TeX('Number of habitats')) + ylab('Precision') + #scale_x_log10() +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}

###########################################################################################
##### Community data relevant plots                                                  ######
###########################################################################################

###########################################################################################
### Master plot for easier exporting
community_plots <- function(run_nms,dist_nms, plotdims, hab_num) {
  labs <- as.data.frame(cbind(run_nms,dist_nms))
  pdc <- pd_cub1[pd_cub1$run_nm %in% labs$run_nms,]
  fdc <- fd_cub1[fd_cub1$run_nm %in% labs$run_nms,]
  p1 <- cub_rch(fdc, labs, plotdims, titltxt = 'Richness', hab_num)
  p2 <- cub_bc(fdc, labs, plotdims, titltxt = 'Dissimilarity', hab_num)
  p3 <- cub_numzeroassocn(pdc, labs, plotdims, titltxt = 'cub_numzeroassocn', hab_num)
  cowplot::plot_grid(p1,p2,p3, nrow = 3)
}
###########################################################################################
### Richness with scale
cub_rch <- function(fdc, labs, plotdims, titltxt, hab_num) {
  fdcx <- fdc
  pdcx <- pdc
  pdcx$xxx <- paste(pdcx$d,pdcx$run_nm)
  fdcx$xxx <- paste(fdcx$d,fdcx$run_nm)
  fdcx$mean_habnum <- pdcx[match(fdcx$xxx, pdcx$xxx),]$mean_habnum
  fdcx$dim <- pdcx[match(fdcx$xxx, pdcx$xxx),]$dims
  fdcx$labs <- labs[match(fdcx$run_nm,labs$run_nms),]$dist_nms
  fdcx$mean_habnum <- as.numeric(as.character(fdcx$mean_habnum))
  fdcx$labs <- factor(fdcx$labs, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(fdcx, aes(x = mean_habnum, y = mn_rch, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() +
      xlab(TeX('Number of habitats')) + ylab('Mean species\nrichness') + # (±sd)
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(fdcx, aes(x = dim, y = mn_rch, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + 
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Number of habitats')) + ylab('Mean species\nrichness') + # (±sd)
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Bray-Curtis dissimilarity with scale
cub_bc <- function(fdc, labs, plotdims, titltxt, hab_num) {
  fdcx <- fdc
  pdcx <- pdc
  pdcx$xxx <- paste(pdcx$d,pdcx$run_nm)
  fdcx$xxx <- paste(fdcx$d,fdcx$run_nm)
  fdcx$mean_habnum <- pdcx[match(fdcx$xxx, pdcx$xxx),]$mean_habnum
  fdcx$dim <- pdcx[match(fdcx$xxx, pdcx$xxx),]$dims  
  fdcx$labs <- labs[match(fdcx$run_nm,labs$run_nms),]$dist_nms
  fdcx$mean_habnum <- as.numeric(as.character(fdcx$mean_habnum))
  fdcx$labs <- factor(fdcx$labs, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(fdcx, aes(x = mean_habnum, y = mn_bc, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() + 
      xlab(TeX('Number of habitats')) + ylab('Mean Bray-Curtis\ndissimilarity') + #(±sd)
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(fdcx, aes(x = dim, y = mn_bc, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() +
      xlab(TeX('Number of habitats')) + ylab('Mean Bray-Curtis\ndissimilarity') + #(±sd)
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }
}
###########################################################################################
### Percentage of zero co-occurrences per run
cub_numzeroassocn <- function(pdc, labs, plotdims, titltxt, hab_num) {
  fdcx <- fdc
  pdcx <- pdc
  pdcx$xxx <- paste(pdcx$d,pdcx$run_nm)
  fdcx$xxx <- paste(fdcx$d,fdcx$run_nm)
  fdcx$mean_habnum <- pdcx[match(fdcx$xxx, pdcx$xxx),]$mean_habnum
  fdcx$dim <- pdcx[match(fdcx$xxx, pdcx$xxx),]$dims  
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$mean_habnum <- as.numeric(as.character(pdcx$mean_habnum))
  pdcx$labs <- factor(pdcx$labs, levels = dist_nms)
  if (hab_num == TRUE) {
    ggplot(pdcx, aes(x = mean_habnum, y = nzeron, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() +
      xlab(TeX('Number of habitats')) + ylab('% runs with zero\nnegative associations') +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  } else if (hab_num == FALSE) {
    ggplot(pdcx, aes(x = mean_habnum, y = nzeron, col = labs, group = labs, fill = labs, shape = labs)) +
      stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
      geom_point(size = 1, stroke = 1) + theme_classic() +
      xlab(TeX('Number of habitats')) + ylab('% runs with zero\nnegative associations') +
      scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
      #scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
      force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
  }  
}


###########################################################################################
###########################################################################################
##### Supporting information Figures                                                  #####
###########################################################################################
###########################################################################################

###########################################################################################
##### SFig4_pairs                                                                     #####
###########################################################################################
get_correlation_driver_map_cubes <- function(infmat_res_comp, dim, plotdims, titltxt) {
  dx <- infmat_res_comp[infmat_res_comp$d == dim,]
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
##### Density plots
### Density for all
get_euc_densities_overall <- function(infmat_res_comp, dim, plotdims, titltxt) {
  dx <- infmat_res_comp[infmat_res_comp$d == dim,]
  res <- c()
  sq <- seq(0,1.6,0.1)
  #i <- 0
  for (i in sq) {
    ndx <- dx[dx$euc >= i & dx$euc < i+0.1,]
    res <- rbind(res, c(nrow(ndx)))
  }
  res <- as.data.frame(res)
  colnames(res) <- c('all_pos_assoc')
  rownames(res) <- -1*sq
  
  res <- as.data.frame(reshape2::melt(as.matrix(res)))
  colnames(res) <- c('bin','type','val')
  res$bin <- as.numeric(res$bin)
  res$bin <- res$bin-0.05
  res$val <- as.numeric(res$val)
  res$val <- res$val/100
  levels(res$type)[levels(res$type)=="all_pos_assoc"] <- "All pairs"
  
  ggplot(res, aes(x = bin)) + coord_flip() + theme_classic() + xlim(-1.5,0) +
    geom_col(aes(y = val), fill = 'black',position="identity", lwd = 0) + ylab('Number of associations') + xlab(TeX('$E_{ij}$')) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
### Density for subgroups
#what_int <- 'noint'
#what_int <- 'justint'
#what_int <- 'all'
get_correlation_driver_densities_cubes <- function(infmat_res_comp, dim, what_int, plotdims, titltxt) {
  dx <- infmat_res_comp[infmat_res_comp$d == dim,]
  
  if (what_int == 'noint') {
    dx <- dx[dx$int == 0,]
  } else if (what_int == 'justint') {
    dx <- dx[dx$int != 0,]
  } else if (what_int == 'all') {
    print('all')
    #dx <- dx[dx$int == 0,]
  }
  
  res <- c()
  sq <- seq(0,1.6,0.1)
  #i <- 0
  for (i in sq) {
    ndx <- dx[dx$euc >= i & dx$euc < i+0.1,]
    res <- rbind(res, c(sum(ndx$coocpos == 1), sum(ndx$coocneg == 1)))
  }
  res <- as.data.frame(res)
  colnames(res) <- c('Positive','Negative')
  rownames(res) <- -1*sq
  res <- as.data.frame(reshape2::melt(as.matrix(res)))
  colnames(res) <- c('bin','type','val')
  res$bin <- as.numeric(res$bin)
  res$bin <- res$bin-0.05
  res$val <- as.numeric(res$val)
  res$val <- res$val/100
  res$typ <- factor(res$type, levels = c('Positive','Negative'))
  
  ggplot(res, aes(x = bin),) + coord_flip() + theme_classic() + xlim(-1.5,0) + ylim(0,1.2) +
    geom_col(aes(y = val, fill = typ, alpha = 0.7),position="identity", lwd = 0) + ylab('Number of associations') + xlab(TeX('$E_{ij}$')) +
    scale_fill_manual(values = c('red','blue')) + labs(color = 'Species pair', shape = 'Resource\nrelationship') +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}


###########################################################################################
##### SFig4_recall_F1                                                                 #####
###########################################################################################

###########################################################################################
### Master plot for better expoting
# Recovery
coocurrence_plots_recovery <- function(run_nms,dist_nms, plotdims, hab_num) {
  labs <- as.data.frame(cbind(run_nms,dist_nms))
  pdc <- pd_cub1[pd_cub1$run_nm %in% labs$run_nms,]
  fdc <- fd_cub1[fd_cub1$run_nm %in% labs$run_nms,]
  p1 <- recovery_intp(pdc, labs, plotdims, nintp = 4, titltxt = 'Recovery pos. interactions', hab_num)
  p2 <- recovery_intn(pdc, labs, plotdims, nintn = 18, titltxt = 'Recovery neg. interactions', hab_num)
  p3 <- recovery_eucp(pdc, labs, plotdims, neuc = 11, titltxt = 'Recovery preference sim.', hab_num)
  p4 <- recovery_eucn(pdc, labs, plotdims, neuc = 11, titltxt = 'Recovery preference dissim.', hab_num)
  cowplot::plot_grid(p1,p2,p3,p4, nrow = 2)
}
# F1
coocurrence_plots_f1 <- function(run_nms,dist_nms, plotdims, hab_num) {
  labs <- as.data.frame(cbind(run_nms,dist_nms))
  pdc <- pd_cub1[pd_cub1$run_nm %in% labs$run_nms,]
  fdc <- fd_cub1[fd_cub1$run_nm %in% labs$run_nms,]
  p1 <- f1_intp(pdc, labs, plotdims, nintp = 4, titltxt = 'Recovery pos. interactions', hab_num)
  p2 <- f1_intn(pdc, labs, plotdims, nintn = 18, titltxt = 'Recovery neg. interactions', hab_num)
  p3 <- f1_eucp(pdc, labs, plotdims, neuc = 11, titltxt = 'Recovery preference sim.', hab_num)
  p4 <- f1_eucn(pdc, labs, plotdims, neuc = 11, titltxt = 'Recovery preference dissim.', hab_num)
  cowplot::plot_grid(p1,p2,p3,p4, nrow = 2)
}


###########################################################################################
##### Recovery
# Positive co-occurrences and mutualsim
recovery_intp <- function(pdc, labs, nintp, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_p <- pdcx$intp*pdcx$cor_totp
  pdcx$recall <- pdcx$perc_p/nintp
  pdcx$f1 <- (2/(pdcx$perc_p+pdcx$recall))
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = recall, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Recovery') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Positive co-occurrences and environmental preference similarity
recovery_eucp <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_p <- pdcx$eucp*pdcx$cor_totp
  pdcx$recall <- pdcx$perc_p/neuc
  pdcx$f1 <- (2/(pdcx$perc_p+pdcx$recall))
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = recall, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Recovery') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative co-occurrences and competition
recovery_intn <- function(pdc, labs, nintn, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_n <- pdcx$intn*pdcx$cor_totn
  pdcx$recall <- pdcx$perc_n/nintn
  pdcx$f1 <- (2/(pdcx$perc_n+pdcx$recall))
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = recall, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Recovery') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative co-occurrences and environmental preference dissimilarity
recovery_eucn <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_n <- pdcx$eucn*pdcx$cor_totn
  pdcx$recall <- pdcx$perc_n/neuc
  pdcx$f1 <- (2/(pdcx$perc_n+pdcx$recall))
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = recall, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('Recovery') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
###########################################################################################
##### F1
# Positive co-occurrences and mutualsim
f1_intp <- function(pdc, labs, nintp, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_p <- pdcx$intp*pdcx$cor_totp
  pdcx$recall <- pdcx$perc_p/nintp
  pdcx$f1 <- 2* ( (pdcx$intp*pdcx$recall) /  (pdcx$intp+pdcx$recall) )
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = f1, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('F1') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Positive co-occurrences and environmental preference similarity
f1_eucp <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_p <- pdcx$eucp*pdcx$cor_totp
  pdcx$recall <- pdcx$perc_p/neuc
  pdcx$f1 <- 2* ( (pdcx$eucp*pdcx$recall) /  (pdcx$eucp+pdcx$recall) )
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = f1, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('F1') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative co-occurrences and competition
f1_intn <- function(pdc, labs, nintn, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_n <- pdcx$intn*pdcx$cor_totn
  pdcx$recall <- pdcx$perc_n/nintn
  pdcx$f1 <- 2* ( (pdcx$intn*pdcx$recall) /  (pdcx$intn+pdcx$recall) )
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = f1, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() +
    xlab(TeX('Composite sample side length l')) + ylab('F1') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}
# Negative co-occurrences and environmental preference dissimilarity
f1_eucn <- function(pdc, labs, neuc, plotdims, titltxt, hab_num) {
  pdcx <- pdc
  pdcx$labs <- labs[match(pdcx$run_nm,labs$run_nms),]$dist_nms
  pdcx$perc_n <- pdcx$eucn*pdcx$cor_totn
  pdcx$recall <- pdcx$perc_n/neuc
  pdcx$f1 <- 2* ( (pdcx$eucn*pdcx$recall) /  (pdcx$eucn+pdcx$recall) )
  pdcx$rscdist <- factor(pdcx$rscdist, levels = dist_nms)
  ggplot(pdcx, aes(x = dim, y = f1, col = labs, group = labs, shape = labs)) +
    stat_smooth(se = FALSE) + #, formula = y~poly(x,2)
    geom_point(size = 1, stroke = 1) + theme_classic() + 
    xlab(TeX('Composite sample side length l')) + ylab('F1') + #scale_x_log10() +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)
}






