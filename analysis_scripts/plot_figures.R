###########################################################################################
### SAMSARA - All plot functions                                                        ###
###########################################################################################
### Read libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(colorspace)
library(shades)
library(ggh4x)
library(ggpattern)
library(sjPlot)
library(latex2exp)
library(eulerr)
library(stringr)
library(rlang)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_analysis_functions'
for (i in list.files(path=('./analysis_scripts/plot_functions'),full.names=TRUE)) {
  source(i)
}

###########################################################################################
##### Fig.2: Simple system                                                            #####
###########################################################################################
### Load data
infm <- read.csv(paste0('./simulation_data/fig2/s1/s1/infmat_res_comp.csv'))
pd_hab <- read.csv(paste0('./analysis_data/fig2/infm_res_compiled.csv'))
fd_hab <- read.csv(paste0('./analysis_data/fig2/full_res_compiled.csv'))
###########################################################################################
##### Main: Fig.2 A) Species pairs co-occurrences
plotdims = c(4,7)
p <- get_correlation_driver_map(infm, plotdims, titltxt)
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/A_all_cors.svg', fig = p, width = 30, height = 20)
###########################################################################################
##### Main: Fig.2 B,C) Matching co-occurrences                                              #####
### Matching (precision)
plotdims = c(4,4)
p <- plot_match_allruns_eucintp(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintp',hnum, plotdims = c(4,4))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/B_match_positive.svg', fig = p, width = 30, height = 20)
p <- plot_match_allruns_eucintn(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintn',hnum, plotdims = c(4,4))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/C_match_negative.svg', fig = p, width = 30, height = 20)
### Main: Recovery
plotdims = c(4,1)
p <- plot_recovery_int(pd_hab, total = c(4,0,4,11,4,4,4), params = c('intp','totp'))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/recovery_int_positive.svg', fig = p)
p <- plot_recovery_int(pd_hab, total = c(18,0,18,11,18,18,18), params = c('intn','totn'))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/recovery_int_negative.svg', fig = p)
p <- plot_recovery_euc(pd_hab, total = c(0,11,11,11,11,11,11), params = c('eucp','totp'))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/recovery_euc_positive.svg', fig = p)
p <- plot_recovery_euc(pd_hab, total = c(0,11,11,11,11,11,11), params = c('eucn','totn'))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/recovery_euc_negative.svg', fig = p)


###########################################################################################
##### Fig.3: Meta-community dynamics                                                  #####
###########################################################################################
pd_hab <- read.csv(paste0('./analysis_data/fig3/infm_res_compiled.csv'))
fd_hab <- read.csv(paste0('./analysis_data/fig3/full_res_compiled.csv'))
###########################################################################################
### Main: Fig. 3 A, B
plotdims = c(3,2.5)
fdh <- fd_hab[fd_hab$treatment %in% c('rc1','rs1','rs2','rs3'),]
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1', 'rs1', 'rs2', 'rs3')
symbolvals <- c(0,2,6,5)
p <- richness_radius(fdh)
sjPlot::save_plot('./figures/main_manuscript/fig3/plots/dispersal_richness.svg', fig = p)
p <- dissimilarity_radius(fdh)
sjPlot::save_plot('./figures/main_manuscript/fig3/plots/dispersal_dissim.svg', fig = p)
###########################################################################################
### Main: Fig. 3 C, D
pdh <- pd_hab[pd_hab$treatment %in% c('rs3'),]
p <- cooc_radius_pos(pdh,titltxt = 'Positive')
sjPlot::save_plot('./figures/main_manuscript/fig3/plots/matching_positive.svg', fig = p)
p <- cooc_radius_neg(pdh,titltxt = 'Negative')
sjPlot::save_plot('./figures/main_manuscript/fig3/plots/matching_negative.svg', fig = p)

###########################################################################################
### Supporting information: Noise effect on co-occurrence
pdh <- pd_hab[pd_hab$treatment %in% c('rc1'),]
p <- cooc_radius_pos(pdh,titltxt = 'Positive rc1')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_positive_rc1.svg', fig = p)
p <- cooc_radius_neg(pdh,titltxt = 'Negative rc1')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_negative_rc1.svg', fig = p)
pdh <- pd_hab[pd_hab$treatment %in% c('rs1'),]
p <- cooc_radius_pos(pdh,titltxt = 'Positive rs1')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_positive_rs1.svg', fig = p)
p <- cooc_radius_neg(pdh,titltxt = 'Negative rs1')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_negative_rs1.svg', fig = p)
pdh <- pd_hab[pd_hab$treatment %in% c('rs2'),]
p <- cooc_radius_pos(pdh,titltxt = 'Positive rs2')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_positive_rs2.svg', fig = p)
p <- cooc_radius_neg(pdh,titltxt = 'Negative rs2')
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_scenarios/plots/matching_negative_rs2.svg', fig = p)
###########################################################################################
### Supporting information: Noise effect on co-occurrence
sim_vec <- c('rs3_noise')
pdh <- pd_hab[pd_hab$treatment %in% c('rs3_noise'),]
p <- cooc_radius_noise_pos(pdh)
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_noise/plots/matching_positive_noise.svg', fig = p, width = 30, height = 20)
p <- cooc_radius_noise_neg(pdh)
sjPlot::save_plot('./figures/supporting_information/sup_dispersal_noise/plots/matching_negative_noise.svg', fig = p, width = 30, height = 20)


###########################################################################################
##### Fig.4: Sampling volume                                                          #####
###########################################################################################
pd_cub <- read.csv(paste0('./analysis_data/fig4/infm_res_compiled.csv'))
fd_cub <- read.csv(paste0('./analysis_data/fig4/full_res_compiled.csv'))

#getOption("bitmapType")
#capabilities("cairo")
#grDevices::dev.cur()

###########################################################################################
### Main: Resource distributions (Plots A, B, C, D)
plotdims <- c(4.,3)
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1', 'rs1', 'rs2', 'rs3')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots(pd_cub, sim_nms, colpal, symbolvals)
save_svg('./figures/main_manuscript/fig4/plots/cooccurrence_rscdist.svg', fig = p, width = 30, height = 20, base_font_size = 10)
p <- community_plots(fd_cub, sim_nms, colpal, symbolvals)
save_svg('./figures/main_manuscript/fig4/plots/community_rscdist.svg', fig = p, width = 30, height = 20)
###########################################################################################
### Main: Alternative scenarios (Plots A, B, C, D)
colpal <- c('#22394AFF', '#AE2565FF','#33C6BBFF', '#666666FF')
sim_nms <- c('as1','as2','as3','rs3')
symbolvals <- c(5,1,3,4) #6,2,5,4,0,3)
p <- coocurrence_plots(pd_cub, sim_nms, colpal, symbolvals)
save_svg('./figures/main_manuscript/fig4/plots/cooccurrence_alternative.svg', fig = p, width = 30, height = 20)
p <- community_plots(fd_cub, sim_nms, colpal, symbolvals)
save_svg('./figures/main_manuscript/fig4/plots/community_alternative.svg', fig = p, width = 30, height = 20)





###########################################################################################
### Test with data variation
pdc <- pd_cub[pd_cub$treatment %in% sim_nms,]
pdc$dim <- as.numeric(unlist(strsplit(unlist(strsplit(pdc$d, '__'))[2*(1:length(pdc$d))], '_'))[2*(1:length(pdc$d))])
ggplot(pdc, aes(x = dim, y = totp, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
  geom_ribbon(aes(ymin = totp - totp_sd, ymax = totp + totp_sd), alpha = 0.1) +
  geom_line() + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1, big.mark = "")) + 
  geom_point(size = 1, stroke = 1) + theme_classic() + 
  xlab(TeX('Composite sample side length l')) + ylab('Positive\nco-occurrences') + #scale_x_log10() +
  scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
  scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
  scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
  scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
  force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)

ggplot(pdc, aes(x = dim, y = intp, col = treatment, group = treatment, fill = treatment, shape = treatment)) +
  geom_ribbon(aes(ymin = intp - intp_sd, ymax = intp + intp_sd), alpha = 0.1) +
  geom_line() + geom_hline(yintercept = nintp/45, linetype = "dashed") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, big.mark = "")) +
  geom_point(size = 1, stroke = 1) + theme_classic() +
  xlab(TeX('Composite sample side length l')) + ylab('Precision') +# scale_x_log10() +
  scale_shape_manual(values = symbolvals, name = 'Resource\ndistribution\nscenarios') +
  scale_color_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
  scale_fill_manual(values = colpal, name = 'Resource\ndistribution\nscenarios') +
  scale_x_continuous(breaks = c(0, seq(0.1, 0.3, by = 0.1)), labels = c("Habitat", seq(0.1, 0.3, by = 0.1))) +
  force_panelsizes(rows = unit(plotdims[1], "cm"),cols = unit(plotdims[2], "cm")) + font_size_control + ggtitle(titltxt)




###########################################################################################
### Supporting information: Simpsons paradox
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1_unequal', 'rs1_unequal', 'rs2_unequal', 'rs3_unequal')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots(pd_cub, sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/supporting_information/sup_unequal_habitats/plots/cooccurrence_rscdist_unequal.svg', fig = p, width = 30, height = 20)
###########################################################################################
### Supporting information: Co-occurrence plots for noise
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1_noise', 'rs1_noise', 'rs2_noise', 'rs3_noise')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots_noise(pd_cub, sim_nms, noise_lvl = 'noise_0.05', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist_noise0.05.svg', fig = p, width = 30, height = 20)
p <- coocurrence_plots_noise(pd_cub, sim_nms, noise_lvl = 'noise_0.1', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist_noise0.1.svg', fig = p, width = 30, height = 20)
p <- coocurrence_plots_noise(pd_cub, sim_nms, noise_lvl = 'noise_0.25', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist_noise0.25.svg', fig = p, width = 30, height = 20)
###########################################################################################
### Supporting information: Species pairwise preference correlations (rs3)
subdir <- "./simulation_data/fig4/rs3/rs3"
p <- get_pref_effect(subdir)
sjPlot::save_plot('./figures/supporting_information/sup_pref_r1/plots/cooccurrence_rscdist.svg', fig = p, width = 30, height = 20)


###########################################################################################
##### Fig.5: Network plots                                                            #####
###########################################################################################
subdir <- "./simulation_data/fig5/nw"
### Get interaction matrices and environmental preference dissimilarities (the same for all three data sets)
intm <- read.csv(paste0(subdir,'/s1/sp_int_mat_log.csv'))[,-1]
eucm <- read.csv(paste0(subdir,'/s1/euc_rsc_prf_sim_mat_log.csv'))[,-1]
pos_thr <- get_Qs(eucm, 0)[1] # Similar environmental preference threshold
neg_thr <- get_Qs(eucm, 0)[2] # Dissimilar environmental preference threshold
### Load co-occurrence data
basic_adj <- read.csv(paste0(subdir,'/s1/sign_ajd_res_log.csv'))[,-1]
basic_sprsc <- read.csv(paste0(subdir,'/s1/sprsc_cor_res_log.csv'))[,-1]
nullint_adj <- read.csv(paste0(subdir,'/c1/sign_ajd_res_log.csv'))[,-1]
nullint_sprsc <- read.csv(paste0(subdir,'/c1/sprsc_cor_res_log.csv'))[,-1]
nullrsc_adj <- read.csv(paste0(subdir,'/c2/sign_ajd_res_log.csv'))[,-1]
nullrsc_sprsc <- read.csv(paste0(subdir,'/c2/sprsc_cor_res_log.csv'))[,-1]
### Select a replicate
repl <- 95
p <- plot_nw_pos_allthree(basic_adj,basic_sprsc,nullint_adj,nullint_sprsc,nullrsc_adj,nullrsc_sprsc,eucm,intm,dset='HabSubsmpl_25',repl)
sjPlot::save_plot('./figures/main_manuscript/fig5/plots/networks.svg', fig = p, width = 30, height = 20)
# Context infor on specialisation (preference for an environmental factor > 0.5)
sprscprf <- read.csv(paste0(subdir,'/s1/sp_rsc_prf_log.csv'))[,-1]
sp_nms <- unique(sprscprf$rsc_nms)
srp <- sprscprf[sprscprf$repl == repl,1:3]
srp[srp < 0.5] <- NA
rownames(srp) <- sp_nms

