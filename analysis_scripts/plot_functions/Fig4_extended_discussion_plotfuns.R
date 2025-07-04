###########################################################################################
##### SAMSARA - Figure 4: Extended discussion - Analysis and plot functions           #####
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
#### Ternary plots of resource abundance distributions
###########################################################################################
### Get resource abundance data for all scenarios
get_resource_abundance_data_of_all_scenarios <- function(parent_set_of_analyses,parentdirnm) {
  allsums <- c()
  #z <- 2
  for (z in 1:length(list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE))) {
    subdir <- list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE)[z]
    subdirnm <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''), full.names=FALSE)[z]
    kabs_cubd_log <- read.csv(paste0(subdir,'/kabs_cubd_log.csv'))[,-1]
    kabs_cubd_log$rowsum <- rowSums(kabs_cubd_log[,1:3])
    kabs_cubd_log$rowsd <- rowSds(as.matrix(kabs_cubd_log[,1:3]))
    kabs_cubd_log$rowmn <- rowMeans(kabs_cubd_log[,1:3])
    kabs_cubd_log$rownormsd <- rowSds(as.matrix(kabs_cubd_log[,1:3]))/rowMeans(as.matrix(kabs_cubd_log[,1:3]))
    kabs_cubd_log$rownormsd_sum <- rowSds(as.matrix(kabs_cubd_log[,1:3]))/rowSums(as.matrix(kabs_cubd_log[,1:3]))
    allsums <- rbind(allsums, cbind(kabs_cubd_log, subdirnm))
  }
  allsums$dim <- as.numeric(str_split_fixed(allsums$d, pattern = '_', 5)[,5])
  allsums$rscdist <- str_split_fixed(allsums$subdirnm, pattern = '_', 5)[,4]
  allsums[allsums$rscdist == 'diffL',]$rscdist <- 'rs3'
  allsums[allsums$rscdist == 'same0',]$rscdist <- 'rc1'
  allsums[allsums$rscdist == 'same5',]$rscdist <- 'rs1'
  allsums[allsums$rscdist == 'same10',]$rscdist <- 'rs2'
  return(allsums)  
}
###########################################################################################
### Ternary plots of all resource abundances per resource scenario -> Colours etc. were changed in inkscape
ternary_resource_abundances <- function(allsums, d, nm) {
  dat <- as.data.frame(allsums[allsums$d == d & allsums$subdirnm == nm,])
  dat[,1:3] <- dat[,1:3]/max(dat[,1:3])
  ggtern(data = dat, aes(x = R1, y = R2, z = R3)) +
    geom_point(color="black",size=1) + theme_rgbw() 
}

###########################################################################################
#### Extract variation in total species and resource abundances for the different resource distribution scenarios
###########################################################################################
get_distribution_of_total_resource_and_species_abundances <- function(parent_set_of_analyses, parentdirnm) {
  res <- c()
  #z <- 2
  for (z in 1:length(list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE))) {
    subdir <- list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE)[z]
    subdirnm <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''), full.names=FALSE)[z]
    npop_cubd_log <- read.csv(paste0(subdir,'/npop_cubd_log.csv'))[,-1]
    kabs_cubd_log <- read.csv(paste0(subdir,'/kabs_cubd_log.csv'))[,-1]
    ###########################################################################################
    for (d in unique(kabs_cubd_log$d)) {
      nuresr <- c()
      #d <- unique(kabs_cubd_log$d)[9]
      #r <- 0
      for (r in unique(kabs_cubd_log$repl)) {
        npop <- npop_cubd_log[npop_cubd_log$d == d & npop_cubd_log$repl == r,]
        kabs <- kabs_cubd_log[kabs_cubd_log$d == d & kabs_cubd_log$repl == r,]
        npop[,1:10] <- sapply(npop[,1:10], as.numeric)
        kabs[,1:3] <- sapply(kabs[,1:3], as.numeric)
        mnpop <- mean(rowSums(npop[,1:10]))
        mnkabs <- mean(rowSums(kabs[,1:3]))
        sdpop <- sd(rowSums(npop[,1:10]))
        sdkabs <- sd(rowSums(kabs[,1:3]))
        mnsdpop <- sd(rowSums(npop[,1:10]))/mean(rowSums(npop[,1:10]))
        mnsdkabs <- sd(rowSums(kabs[,1:3]))/mean(rowSums(kabs[,1:3]))
        nuresr <- rbind(nuresr, c(mnpop=mnpop,mnkabs=mnkabs,sdpop=sdpop,sdkabs=sdkabs,mnsdpop=mnsdpop,mnsdkabs=mnsdkabs))
      }
      res <- rbind(res, c(runm = subdirnm, d = d, colMeans(nuresr)))
      print(paste('Done for',d,z))
    }
  }  
  res <- as.data.frame(res)
  res[,3:ncol(res)] <- sapply(res[,3:ncol(res)], as.numeric)
  res$dim <- as.numeric(str_split_fixed(res$d, pattern = '_', 5)[,5])
  res$rscdist <- str_split_fixed(res$runm, pattern = '_', 5)[,4]
  res$rscdist <- factor(res$rscdist, levels = c(''))
  return(res)
}
###########################################################################################
### Corresponding plot function
spc_rsc_abund <- function(res, colnm, ylabl) {
  res$plot_this <- res[,colnames(res) == colnm]
  ggplot(res, aes(x = dim, y = plot_this, col = labs, group = labs, shape = labs)) +
    geom_point(size = 1, stroke = 1) + theme_classic() + stat_smooth(se = FALSE) +#geom_line() + 
    xlab(TeX('Number of habitats')) + ylab(ylabl) +
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_fill_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)  
}

###########################################################################################
##### Specialism/generalsim in associations                                           #####
###########################################################################################
get_species_specialsim_info <- function(parent_set_of_analyses,parentdirnm) {
  standard_error <- function(x) {
    sd(x) / sqrt(length(x))
  }
  Gamma_resl <- c()
  Gres_pos <- c()
  Gres_neg <- c()
  ###########################################################################################
  #z <- 1
  for (z in 1:length(list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE))) {
    subdir <- list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE)[z]
    subdirnm <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''), full.names=FALSE)[z]
    sp_correl_log <- read.csv(paste0(subdir,'/sp_correl_log.csv'))[,-1]
    infmat_res_comp <- read.csv(paste0(subdir,'/infmat_res_comp.csv'))[,-1]
    infm <- infmat_res_comp 
    sp_rsc_prf <- read.csv(paste0(subdir,'/sp_rsc_prf_log.csv'))[,-1]
    ###########################################################################################
    ### Adding species preferences to association data
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
      Gres_posd <- c()
      Gres_negd <- c()
      #r <- 1
      for (r in unique(infm$repl)) {
        Gamma_res <- c()
        infmr <- infmd[infmd$repl == r,]
        srprfr <- sp_rsc_prf[sp_rsc_prf$repl == r,]
        #i <- 1
        for (i in 1:nrow(infmr)) {
          ###########################################################################################
          ### Get euc distance between both species
          Gij <- sqrt((srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R1 - srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R1)^2 + 
                        (srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R2 - srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R2)^2 + 
                        (srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$R3 - srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$R3)^2)
          ###########################################################################################
          ### Check which combination both species have (all combos documented except two different types of specialists = R_R)
          if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1'){
            RGij <- 'both_R1'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R2' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R2'){
            RGij <- 'both_R2'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R3' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R3'){
            RGij <- 'both_R3'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'G' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'G'){
            RGij <- 'both_G'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'G' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'G'){
            RGij <- 'R1_G'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'G' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R2' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R2' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'G'){
            RGij <- 'R2_G'
          } else if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'G' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R3' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R3' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'G'){
            RGij <- 'R3_G'
          } else {
            RGij <- 'R_R'
          }
          
          #R1_vs_R2R3ij <- 'Nope'
          #if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R2' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R2' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R3' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R3' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1'){
          #  R1_vs_R2R3ij <- 'R1_vs_R2R3'
          #}  
          R_vs_R <- 'Nope'
          if(srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R2' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R3' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R3' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R2'){
            R_vs_R <- 'R2_vs_R3'
          } else if (srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R2' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R2' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1') {
            R_vs_R <- 'R1_vs_R2'
          } else if (srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R1' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R3' | srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen == 'R3' & srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen == 'R1') {
            R_vs_R <- 'R1_vs_R3'
          }
          ###########################################################################################
          ### Save 
          Gamma_res <- rbind(Gamma_res, c(Ni = infmr[i,]$N1, Nj = infmr[i,]$N2, RGij = RGij, R_vs_R = R_vs_R,
                                          SPGENi = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$spgen, SPGENj = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$spgen,
                                          GR1i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR1, GR1j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR1,
                                          GR2i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR2, GR2j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR2,
                                          GR3i = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GR3, GR3j = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GR3,
                                          GGi = srprfr[srprfr$rsc_nms == infmr[i,]$N1,]$GG, GGj = srprfr[srprfr$rsc_nms == infmr[i,]$N2,]$GG,
                                          Gij = Gij, cooc = infmr[i,]$cooc, int_i = infmr[i,]$int))
        }
        ### Log file for all repls
        #Gamma_resl <- rbind(Gamma_resl, cbind(Gamma_res, repl = r, d = d, run_nm = subdirnm))      
        ### Aggregate data for d and runs 
        Gresx <- as.data.frame(Gamma_res)
        Gresx[,7:ncol(Gresx)] <- sapply(Gresx[,7:ncol(Gresx)], as.numeric)
        Gresx$RGij <- factor(Gresx$RGij, levels = c('both_R1','both_R2','both_R3','both_G','R1_G','R2_G','R3_G','R_R'))
        Gresx$R_vs_R <- factor(Gresx$R_vs_R, levels = c('Nope','R2_vs_R3','R1_vs_R2','R1_vs_R3'))
        ### Results for positive associations
        Gpos <- Gresx[Gresx$cooc > 0,] 
        tRGij <- table(Gpos$RGij)/nrow(Gpos)#sum(tRGij)
        tR_vs_R <- table(Gpos$R_vs_R)/nrow(Gpos)
        dist_res <- cbind(mnGR1 = (Gpos$GR1i+Gpos$GR1j)/2, mnGR2 = (Gpos$GR2i+Gpos$GR2j)/2, 
                          mnGR3 = (Gpos$GR3i+Gpos$GR3j)/2, mnGG = (Gpos$GGi+Gpos$GGj)/2)
        sedist_res <- c(seGR1 = standard_error((Gpos$GR1i+Gpos$GR1j)/2), seGR2 = standard_error((Gpos$GR2i+Gpos$GR2j)/2), 
                        seGR3 = standard_error((Gpos$GR3i+Gpos$GR3j)/2), seGG = standard_error((Gpos$GGi+Gpos$GGj)/2))
        Gres_posd <- rbind(Gres_posd, c(nrow(Gpos), colSums(Gpos[,7:ncol(Gpos)])/nrow(Gpos), colSums(dist_res)/nrow(Gpos), sedist_res, tRGij, tR_vs_R))
        ### Results for negative associations
        Gneg <- Gresx[Gresx$cooc < 0,] 
        tRGij <- table(Gneg$RGij)/nrow(Gneg)#sum(tRGij)
        tR_vs_R <- table(Gneg$R_vs_R)/nrow(Gneg)
        dist_res <- cbind(mnGR1 = (Gneg$GR1i+Gneg$GR1j)/2, mnGR2 = (Gneg$GR2i+Gneg$GR2j)/2, 
                          mnGR3 = (Gneg$GR3i+Gneg$GR3j)/2, mnGG = (Gneg$GGi+Gneg$GGj)/2)
        sedist_res <- c(seGR1 = standard_error((Gneg$GR1i+Gneg$GR1j)/2), seGR2 = standard_error((Gneg$GR2i+Gneg$GR2j)/2), 
                        seGR3 = standard_error((Gneg$GR3i+Gneg$GR3j)/2), seGG = standard_error((Gneg$GGi+Gneg$GGj)/2))
        Gres_negd <- rbind(Gres_negd, c(nrow(Gneg), colSums(Gneg[,7:ncol(Gneg)])/nrow(Gneg), colSums(dist_res)/nrow(Gneg), sedist_res, tRGij, tR_vs_R))
      }
      ###########################################################################################
      ### Get means for all replicates
      Gres_posd <- Gres_posd[Gres_posd[,1] > 0,] 
      Gres_negd <- Gres_negd[Gres_negd[,1] > 0,] 
      
      Gres_pos <- rbind(Gres_pos, c(run_nm = subdirnm, d = d, colSums(Gres_posd)/nrow(Gres_posd)))
      Gres_neg <- rbind(Gres_neg, c(run_nm = subdirnm, d = d, colSums(Gres_negd)/nrow(Gres_negd)))
      print(paste('Done for',d,subdirnm))
    }
  }
  return(list(Gamma_res = Gamma_res, Gres_pos = Gres_pos, Gres_neg = Gres_neg))
}
###########################################################################################
### Prepare data for plotting
prep_specialsim_info <- function(Gres_pos,Gres_neg) {
  Gres_pos <- as.data.frame(Gres_pos)
  Gres_neg <- as.data.frame(Gres_neg)
  Gres_pos[,3:ncol(Gres_pos)] <- sapply(Gres_pos[,3:ncol(Gres_pos)],as.numeric)
  Gres_neg[,3:ncol(Gres_neg)] <- sapply(Gres_neg[,3:ncol(Gres_neg)],as.numeric)
  Gres_pos$dim <- as.numeric(str_split_fixed(Gres_pos$d, pattern = '_', 5)[,5])
  Gres_neg$dim <- as.numeric(str_split_fixed(Gres_neg$d, pattern = '_', 5)[,5])
  Gres_pos$rscdist <- str_split_fixed(Gres_pos$run_nm, pattern = '_', 5)[,4]
  Gres_neg$rscdist <- str_split_fixed(Gres_neg$run_nm, pattern = '_', 5)[,4]
  Gres_pos$rscdist <- factor(Gres_pos$rscdist, levels = dist_nms)
  Gres_neg$rscdist <- factor(Gres_neg$rscdist, levels = dist_nms)
  return(list(Gres_pos,Gres_neg))
}
###########################################################################################
### Plotting function
#Gres <- Gres_pos
#Gres <- Gres_neg
pairs_generalsim <- function(Gres) {
  ggplot(Gres, aes(x = dim, y = mnGG, col = rscdist, group = rscdist, shape = rscdist)) +
    geom_point(size = 1, stroke = 1) + theme_classic() + stat_smooth(method = 'lm', se = FALSE) +#geom_line() + 
    #geom_errorbar(aes(ymin = mnGG-seGG, ymax = mnGG+seGG), width = 0.01) + 
    xlab(TeX('Sample cube side length l')) + ylab('Pairs mean\ndistance to G*') + 
    scale_shape_manual(values = symbolvals, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_color_manual(values = colpal_rscscenarios, labels = xlabels, name = 'Resource\ndistribution\nscenarios') +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    theme(legend.position = "none") +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + ggtitle(titltxt)#+ font_size_control 
}


###########################################################################################
###########################################################################################
##### Preference effects on correlations and  co-occurrence                           #####
###########################################################################################
###########################################################################################

###########################################################################################
##### Pairwise preferences and co-occurrence
###########################################################################################
### Get data on species pairwise preferences and co-occurrence
get_preferences_and_cooccurrence <- function(parent_set_of_analyses,parentdirnm) {
  infmres <- c()
  res_pos <- c()
  res_neg <- c()
  for (z in 1:length(list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE))) {
    subdir <- list.dirs(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''),recursive = FALSE)[z]
    subdirnm <- list.files(path=paste(workdir,"/Result_master_dir/",parent_set_of_analyses,'/',parentdirnm,sep=''), full.names=FALSE)[z]
    
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
        
        infmres <- rbind(infmres, cbind(runm = subdirnm, infmr))
        
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
      res_pos <- rbind(res_pos, c(runm = subdirnm, d = d, colMeans(res_posr, na.rm = TRUE)))
      res_neg <- rbind(res_neg, c(runm = subdirnm, d = d, colMeans(res_negr, na.rm = TRUE)))
      print(paste('Done for',d,subdirnm))
    }
  }
  return(list(infmres = infmres, res_pos = res_pos, res_neg = res_neg))
}

###########################################################################################
##### Plotting pairwise preferences for a single scneario and sampling scale
###########################################################################################

###########################################################################################
### Prepare for plotting
prep_preferences_and_cooccurrence_data <- function(infmres) {
  res <- as.data.frame(infmres)
  res$dim <- as.numeric(str_split_fixed(res$d, pattern = '_', 5)[,5])
  res$rscdist <- str_split_fixed(res$runm, pattern = '_', 5)[,4]
  res[,c(4:15, 17:30)] <- sapply(res[,c(4:15, 17:30)],as.numeric)
  res[res$rscdist == 'same0',]$rscdist <- 'rc1' 
  res[res$rscdist == 'same5',]$rscdist <- 'rs1' 
  res[res$rscdist == 'same10',]$rscdist <- 'rs2' 
  res[res$rscdist == 'diffL',]$rscdist <- 'rs3' 
  #res$rscdist <- factor(res$rscdist, levels = c('rc1','rs1','rs2','rs3'))
  res$intbin <- 'Other'
  res[res$int>0,]$intbin <- 'Competition'
  res[res$int<0,]$intbin <- 'Mutualism'
  pres <- res[res$corcoeff > 0,]
  nres <- res[res$corcoeff < 0,]
  return(list(pres = pres, nres = nres))  
}
###########################################################################################
##### Plot pairwise preferences
### Postive co-occurrence
plt_pos <- function(pres, dim, corthr, x, y) {
  presx <- pres[pres$dim == dim & pres$corcoeff > corthr,]
  if (x != y) {
    presxy <- as.data.frame(rbind(cbind(X = presx[,colnames(presx) == paste0(x,'i')], Y = presx[,colnames(presx) == paste0(y,'j')],
                                        Interaction = presx$intbin, rscdist = presx$rscdist),
                                  cbind(X = presx[,colnames(presx) == paste0(x,'j')], Y = presx[,colnames(presx) == paste0(y,'i')],
                                        Interaction = presx$intbin, rscdist = presx$rscdist)))    
  } else {
    presxy <- as.data.frame(cbind(X = presx[,colnames(presx) == paste0(x,'i')], Y = presx[,colnames(presx) == paste0(y,'j')],
                                  Interaction = presx$intbin, rscdist = presx$rscdist))
  }
  presxy[,1:2] <- sapply(presxy[,1:2], as.numeric)
  labls <- c()
  #i <- unique(presxy$rscdist)[1]
  for (i in unique(presxy$rscdist)) {
    labls <- as.data.frame(rbind(labls, c(rscdist = i, corelc = cor(presxy[presxy$rscdist == i,]$X,presxy[presxy$rscdist == i,]$Y, method = 'spearman'))))
  }
  labls$corelc <- as.numeric(labls$corelc)
  labls$lab <- paste0(labls$rscdist,' (',round(labls$corelc,2),')')
  presxy <- merge(presxy, labls[, c("lab","rscdist")], by = "rscdist")
  color_palette <- c("Other" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(presxy, aes(x = X, y = Y, col = Interaction)) +
    geom_point(alpha = 0.2) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ lab) + theme_bw() + xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    scale_color_manual(values = color_palette) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold'))
}
### Negative co-occurrence
plt_neg <- function(nres, dim, corthr, x, y) {
  nresx <- nres[nres$dim == dim & nres$corcoeff < corthr,]
  nresx[nresx$int<0,]$intbin <- 'Other'
  if (x != y) {
    nresxy <- as.data.frame(rbind(cbind(X = nresx[,colnames(nresx) == paste0(x,'i')], Y = nresx[,colnames(nresx) == paste0(y,'j')],
                                        Interaction = nresx$intbin, rscdist = nresx$rscdist),
                                  cbind(X = nresx[,colnames(nresx) == paste0(x,'j')], Y = nresx[,colnames(nresx) == paste0(y,'i')],
                                        Interaction = nresx$intbin, rscdist = nresx$rscdist)))    
  } else {
    nresxy <- as.data.frame(cbind(X = nresx[,colnames(nresx) == paste0(x,'i')], Y = nresx[,colnames(nresx) == paste0(y,'j')],
                                  Interaction = nresx$intbin, rscdist = nresx$rscdist))
  }
  nresxy[,1:2] <- sapply(nresxy[,1:2], as.numeric)
  labls <- c()
  #i <- unique(nresxy$rscdist)[1]
  for (i in unique(nresxy$rscdist)) {
    labls <- as.data.frame(rbind(labls, c(rscdist = i, corelc = cor(nresxy[nresxy$rscdist == i,]$X,nresxy[nresxy$rscdist == i,]$Y, method = 'spearman'))))
  }
  labls$corelc <- as.numeric(labls$corelc)
  labls$lab <- paste0(labls$rscdist,' (',round(labls$corelc,2),')')
  nresxy <- merge(nresxy, labls[, c("lab","rscdist")], by = "rscdist")
  color_palette <- c("Other" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(nresxy, aes(x = X, y = Y, col = Interaction)) +
    geom_point(alpha = 0.2) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ lab) + theme_bw() + xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
    scale_color_manual(values = color_palette) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold'))
}

###########################################################################################
##### Pairwise preference correlation coefficients for different scenarios and at different scales
###########################################################################################

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
prep_d <- function(data) {
  dat <- as.data.frame(data)
  dat$dim <- as.numeric(str_split_fixed(dat$d, pattern = '_', 5)[,5])
  dat$rscdist <- str_split_fixed(dat$runm, pattern = '_', 5)[,4]
  dat[,4:10] <- sapply(dat[,4:10],as.numeric)
  return(dat)
}
###########################################################################################
plt_corcoeffs_pos <- function(pdf, pltcol) {
  pdfx <- pdf
  pdfx[pdfx$rscdist == 'same0',]$rscdist <- 'rc1' 
  pdfx[pdfx$rscdist == 'same5',]$rscdist <- 'rs1' 
  pdfx[pdfx$rscdist == 'same10',]$rscdist <- 'rs2' 
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  pdfx$pltcol <- pdfx[,colnames(pdfx) == pltcol]
  ggplot(pdfx, aes(x = dim, y = pltcol, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + geom_line() + #stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ rscdist) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans correlation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}
# 
plt_corcoeffs_pos_rs3_samepref <- function(pdf) {
  pdfx <- pdf
  pdfx <- pdfx[pdfx$rscdist == 'diffL',]
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
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
# 
plt_corcoeffs_pos_rs3_difpref <- function(pdf) {
  pdfx <- pdf
  pdfx <- pdfx[pdfx$rscdist == 'diffL',]
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1_R1 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R1_R2')], 'R1 x R2')
  R2_R2 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R1_R3')], 'R1 x R3')
  R3_R3 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R2_R3')], 'R2 x R3')
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

###########################################################################################
plt_corcoeffs_neg <- function(ndf, pltcol) {
  ndfx <- ndf
  ndfx[ndfx$rscdist == 'same0',]$rscdist <- 'rc1' 
  ndfx[ndfx$rscdist == 'same5',]$rscdist <- 'rs1' 
  ndfx[ndfx$rscdist == 'same10',]$rscdist <- 'rs2' 
  ndfx[ndfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  color_palette_n <- c("All" = "black", "No interaction" = "darkgrey", "Competition" = "blue")
  ndfx$pltcol <- ndfx[,colnames(ndfx) == pltcol]
  ggplot(ndfx, aes(x = dim, y = pltcol, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + geom_line() + #stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ rscdist) + theme_bw() + 
    xlab('Composite sample sidelength l') + ylab('Spearmans correlation coefficient') +
    scale_color_manual(values = color_palette_n) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}
# 
plt_corcoeffs_neg_rs3_samepref <- function(ndf) {
  ndfx <- ndf
  ndfx <- ndfx[ndfx$rscdist == 'diffL',]
  ndfx[ndfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1_R1 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R1_R1')], 'R1 x R1')
  R2_R2 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R2_R2')], 'R2 x R2')
  R3_R3 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R3_R3')], 'R3 x R3')
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
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}
# 
plt_corcoeffs_neg_rs3_difpref <- function(ndf) {
  ndfx <- ndf
  ndfx <- ndfx[ndfx$rscdist == 'diffL',]
  ndfx[ndfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1_R1 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R1_R2')], 'R1 x R2')
  R2_R2 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R1_R3')], 'R1 x R3')
  R3_R3 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R2_R3')], 'R2 x R3')
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
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}

###########################################################################################
### Analysis - Individual species preferences and correlation coefficients
###########################################################################################

###########################################################################################
### Collect all the data needed
get_data_for_individual_species_preference <- function(pres, nres) {
  res_p_all <- c()
  res_p_mut <- c()
  res_p_comp <- c()
  res_p_nonint <- c()
  
  res_n_all <- c()
  res_n_comp <- c()
  res_n_nonint <- c()
  #rnm <- unique(pres$runm)[1]
  for (rnm in unique(pres$runm)) {
    #d <- unique(pres$d)[9]
    for (d in unique(pres$d)) {
      ##############################################################
      ### Positive
      # All pairs
      presx <- pres[pres$runm == rnm & pres$d == d,] # & pres$corcoeff > 0.7
      R1i = cor(presx$R1i, presx$corcoeff, method = 'spearman')
      R1j = cor(presx$R1j, presx$corcoeff, method = 'spearman')
      R2i = cor(presx$R2i, presx$corcoeff, method = 'spearman')
      R2j = cor(presx$R2j, presx$corcoeff, method = 'spearman')
      R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman')
      R3j = cor(presx$R3j, presx$corcoeff, method = 'spearman')
      resp <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_p_all <- rbind(res_p_all, c(runm = rnm, d = d, type = 'All', resp))
      # Only mutualistic pairs
      presx <- pres[pres$runm == rnm & pres$d == d & pres$int < 0,] #R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman'),
      R1i = cor(presx$R1i, presx$corcoeff, method = 'spearman')
      R1j = cor(presx$R1j, presx$corcoeff, method = 'spearman')
      R2i = cor(presx$R2i, presx$corcoeff, method = 'spearman')
      R2j = cor(presx$R2j, presx$corcoeff, method = 'spearman')
      R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman')
      R3j = cor(presx$R3j, presx$corcoeff, method = 'spearman')
      resp <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_p_mut <- rbind(res_p_mut, c(runm = rnm, d = d, type = 'Mutualism', resp))
      
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7 & pres$int > 0,]
      R1i = cor(presx$R1i, presx$corcoeff, method = 'spearman')
      R1j = cor(presx$R1j, presx$corcoeff, method = 'spearman')
      R2i = cor(presx$R2i, presx$corcoeff, method = 'spearman')
      R2j = cor(presx$R2j, presx$corcoeff, method = 'spearman')
      R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman')
      R3j = cor(presx$R3j, presx$corcoeff, method = 'spearman')
      resp <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_p_comp <- rbind(res_p_comp, c(runm = rnm, d = d, type = 'Competition', resp))
      
      presx <- pres[pres$runm == rnm & pres$d == d & pres$corcoeff > 0.7 & pres$int == 0,]
      R1i = cor(presx$R1i, presx$corcoeff, method = 'spearman')
      R1j = cor(presx$R1j, presx$corcoeff, method = 'spearman')
      R2i = cor(presx$R2i, presx$corcoeff, method = 'spearman')
      R2j = cor(presx$R2j, presx$corcoeff, method = 'spearman')
      R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman')
      R3j = cor(presx$R3j, presx$corcoeff, method = 'spearman')
      resp <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_p_nonint <- rbind(res_p_nonint, c(runm = rnm, d = d, type = 'No interaction', resp))
      
      ##############################################################
      ### Positive
      # All pairs
      nresx <- nres[nres$runm == rnm & nres$d == d,] # & nres$corcoeff < -0.7
      
      #xx = nresx$R1i
      #yy = abs(nresx$corcoeff)
      #cor(xx, yy, method = 'spearman')
      #ggplot(nresx, aes(x = xx, y = yy)) +
      #  geom_point() + stat_smooth(method = 'lm')

      R1i = cor(nresx$R1i, abs(nresx$corcoeff), method = 'spearman')
      R1j = cor(nresx$R1j, abs(nresx$corcoeff), method = 'spearman')
      R2i = cor(nresx$R2i, abs(nresx$corcoeff), method = 'spearman')
      R2j = cor(nresx$R2j, abs(nresx$corcoeff), method = 'spearman')
      R3i = cor(nresx$R3i, abs(nresx$corcoeff), method = 'spearman')
      R3j = cor(nresx$R3j, abs(nresx$corcoeff), method = 'spearman')
      resn <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_n_all <- rbind(res_n_all, c(runm = rnm, d = d, type = 'All', resn))
      # Only mutualistic pairs
      nresx <- nres[nres$runm == rnm & nres$d == d & nres$corcoeff < -0.7 & nres$int > 0,]
      R1i = cor(nresx$R1i, abs(nresx$corcoeff), method = 'spearman')
      R1j = cor(nresx$R1j, abs(nresx$corcoeff), method = 'spearman')
      R2i = cor(nresx$R2i, abs(nresx$corcoeff), method = 'spearman')
      R2j = cor(nresx$R2j, abs(nresx$corcoeff), method = 'spearman')
      R3i = cor(nresx$R3i, abs(nresx$corcoeff), method = 'spearman')
      R3j = cor(nresx$R3j, abs(nresx$corcoeff), method = 'spearman')
      resn <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_n_comp <- rbind(res_n_comp, c(runm = rnm, d = d, type = 'Competition', resn))
      
      nresx <- nres[nres$runm == rnm & nres$d == d & nres$corcoeff < -0.7 & nres$int == 0,]
      R1i = cor(nresx$R1i, abs(nresx$corcoeff), method = 'spearman')
      R1j = cor(nresx$R1j, abs(nresx$corcoeff), method = 'spearman')
      R2i = cor(nresx$R2i, abs(nresx$corcoeff), method = 'spearman')
      R2j = cor(nresx$R2j, abs(nresx$corcoeff), method = 'spearman')
      R3i = cor(nresx$R3i, abs(nresx$corcoeff), method = 'spearman')
      R3j = cor(nresx$R3j, abs(nresx$corcoeff), method = 'spearman')
      resn <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_n_nonint <- rbind(res_n_nonint, c(runm = rnm, d = d, type = 'No interaction', resn))
      
      print(paste0('Done ',rnm, d))
    }
  }
  pdf <- rbind(prep_indiv_d(res_p_all))#,prep_indiv_d(res_p_mut),prep_indiv_d(res_p_comp),prep_indiv_d(res_p_nonint))
  ndf <- rbind(prep_indiv_d(res_n_all))#,prep_indiv_d(res_n_comp),prep_indiv_d(res_n_nonint))
  return(list(pdf = pdf, ndf = ndf))
}
###########################################################################################
### Prepare data for plotting
prep_indiv_d <- function(data) {
  dat <- as.data.frame(data)
  dat$dim <- as.numeric(str_split_fixed(dat$d, pattern = '_', 5)[,5])
  dat$rscdist <- str_split_fixed(dat$runm, pattern = '_', 5)[,4]
  dat[,4:7] <- sapply(dat[,4:7],as.numeric)
  return(dat)
}
###########################################################################################
### Plotting individual preference correlation with pairwise abundance correlaitons - positive
plt_corcoeffs_pos <- function(pdf, pltcol) {
  pdfx <- pdf
  pdfx[pdfx$rscdist == 'same0',]$rscdist <- 'rc1' 
  pdfx[pdfx$rscdist == 'same5',]$rscdist <- 'rs1' 
  pdfx[pdfx$rscdist == 'same10',]$rscdist <- 'rs2' 
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  pdfx$pltcol <- pdfx[,colnames(pdfx) == pltcol]
  ggplot(pdfx, aes(x = dim, y = pltcol, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ rscdist) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}
#
plt_corcoeffs_rs3_pos <- function(pdf, pltcol) {
  pdfx <- pdf
  pdfx <- pdfx[pdfx$rscdist == 'diffL',]
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R1')], 'R1')
  R2 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R2')], 'R2')
  R3 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R3')], 'R3')
  colnames(R1) <- c('runm','type','comp','d','spec')
  colnames(R2) <- c('runm','type','comp','d','spec')
  colnames(R3) <- c('runm','type','comp','d','spec')
  df <- rbind(R1,R2,R3)
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(df, aes(x = d, y = comp, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ spec) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}


###########################################################################################
### Plotting individual preference correlation with pairwise abundance correlaitons - negative
plt_corcoeffs_neg <- function(ndf, pltcol) {
  ndfx <- ndf
  ndfx[ndfx$rscdist == 'same0',]$rscdist <- 'rc1' 
  ndfx[ndfx$rscdist == 'same5',]$rscdist <- 'rs1' 
  ndfx[ndfx$rscdist == 'same10',]$rscdist <- 'rs2' 
  ndfx[ndfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  color_palette_n <- c("All" = "black", "No interaction" = "darkgrey", "Competition" = "blue")
  ndfx$pltcol <- ndfx[,colnames(ndfx) == pltcol]
  ggplot(ndfx, aes(x = dim, y = pltcol, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ rscdist) + theme_bw() + 
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_n) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}
plt_corcoeffs_rs3_neg <- function(pdf, pltcol) {
  ndfx <- ndf
  ndfx <- ndfx[ndfx$rscdist == 'diffL',]
  ndfx[ndfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R1')], 'R1')
  R2 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R2')], 'R2')
  R3 <- cbind(ndfx[,colnames(ndfx) %in% c('runm','dim','type','R3')], 'R3')
  colnames(R1) <- c('runm','type','comp','d','spec')
  colnames(R2) <- c('runm','type','comp','d','spec')
  colnames(R3) <- c('runm','type','comp','d','spec')
  df <- rbind(R1,R2,R3)
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(df, aes(x = d, y = comp, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ spec) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="blue"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}










###########################################################################################
### Collect all the data needed
get_data_for_individual_species_preference_all_correlations <- function(infmres) {
  res_all <- c()
  #rnm <- unique(infmres$runm)[1]
  for (rnm in unique(infmres$runm)) {
    #d <- unique(infmres$d)[9]
    for (d in unique(infmres$d)) {
      ##############################################################
      ### Positive
      # All pairs
      presx <- infmres[infmres$runm == rnm & infmres$d == d,] # & pres$corcoeff > 0.7
      
      xx = presx$R1i
      yy = presx$corcoeff
      cor(xx, yy, method = 'spearman')
      ggplot(presx, aes(x = xx, y = yy)) +
        geom_point() + stat_smooth(method = 'lm')
      
      R1i = cor(presx$R1i, presx$corcoeff, method = 'spearman')
      R1j = cor(presx$R1j, presx$corcoeff, method = 'spearman')
      R2i = cor(presx$R2i, presx$corcoeff, method = 'spearman')
      R2j = cor(presx$R2j, presx$corcoeff, method = 'spearman')
      R3i = cor(presx$R3i, presx$corcoeff, method = 'spearman')
      R3j = cor(presx$R3j, presx$corcoeff, method = 'spearman')
      resp <- c(R1 = (R1i+R1j)/2, R2 = (R2i+R2j)/2, R3 = (R3i+R3j)/2)
      res_all <- rbind(res_p_all, c(runm = rnm, d = d, type = 'All', resp))

      print(paste0('Done ',rnm, d))
    }
  }
  pdf <- rbind(prep_indiv_d(res_all))#,prep_indiv_d(res_p_mut),prep_indiv_d(res_p_comp),prep_indiv_d(res_p_nonint))
  return(list(pdf = pdf))
}

###########################################################################################
### Collect all the data needed
plt_corcoeffs_rs3_all <- function(pdf, pltcol) {
  pdfx <- pdf
  pdfx <- pdfx[pdfx$rscdist == 'diffL',]
  pdfx[pdfx$rscdist == 'diffL',]$rscdist <- 'rs3' 
  
  R1 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R1')], 'R1')
  R2 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R2')], 'R2')
  R3 <- cbind(pdfx[,colnames(pdfx) %in% c('runm','dim','type','R3')], 'R3')
  colnames(R1) <- c('runm','type','comp','d','spec')
  colnames(R2) <- c('runm','type','comp','d','spec')
  colnames(R3) <- c('runm','type','comp','d','spec')
  df <- rbind(R1,R2,R3)
  color_palette_p <- c("All" = "black", "No interaction" = "darkgrey", "Mutualism" = "red", "Competition" = "blue")
  ggplot(df, aes(x = d, y = comp, col = type)) +
    geom_point(alpha = 1) + geom_hline(yintercept = 0, linetype = "dashed") + stat_smooth(method = 'lm', formula = y~poly(x,1), se = FALSE) + #, formula = y~poly(x,2)
    facet_wrap(~ spec) + theme_bw() + #xlab(latex2exp::TeX(paste0('$gamma_{i,k=',x,'}$'))) + ylab(latex2exp::TeX(paste0('$gamma_{j,k=',y,'}$'))) +
    xlab('Composite sample sidelength l') + ylab('Spearmans\ncorrelation coefficient') +
    scale_color_manual(values = color_palette_p) +
    scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
    force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) +
    theme(strip.background=element_rect(colour="black",fill="red"), strip.text = element_text(colour = 'white', face = 'bold')) + font_size_control
}

