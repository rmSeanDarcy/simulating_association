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
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/A_all_cors.svg', fig = p)
###########################################################################################
##### Main: Fig.2 B,C) Matching co-occurrences                                              #####
### Matching (precision)
plotdims = c(4,4)
p <- plot_match_allruns_eucintp(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintp',hnum, plotdims = c(4,4))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/B_match_positive.svg', fig = p)
p <- plot_match_allruns_eucintn(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintn',hnum, plotdims = c(4,4))
sjPlot::save_plot('./figures/main_manuscript/fig2/plots/C_match_negative.svg', fig = p)
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
### Supporting information: Simpsons paradox
sim_vec <- c('rs3_noise')
pdh <- pd_hab[pd_hab$treatment %in% c('rs3_noise') & pd_hab$simulation  ,]






###########################################################################################
##### Fig.4: Sampling volume                                                          #####
###########################################################################################
pd_cub <- read.csv(paste0('./analysis_data/fig4/infm_res_compiled.csv'))
fd_cub <- read.csv(paste0('./analysis_data/fig4/full_res_compiled.csv'))

###########################################################################################
### Main: Resource distributions (Plots A, B, C, D)
plotdims <- c(3.,2.5)
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1', 'rs1', 'rs2', 'rs3')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots(sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/main_manuscript/fig4/plots/cooccurrence_rscdist.svg', fig = p)
p <- community_plots(sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/main_manuscript/fig4/plots/community_rscdist.svg', fig = p)
###########################################################################################
### Main: Alternative scenarios (Plots A, B, C, D)
colpal <- c('#22394AFF', '#AE2565FF','#33C6BBFF', '#666666FF')
sim_nms <- c('as1','as2','as3','rs3')
symbolvals <- c(5,1,3,4) #6,2,5,4,0,3)
p <- coocurrence_plots(sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/main_manuscript/fig4/plots/cooccurrence_alternative.svg', fig = p)
p <- community_plots(sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/main_manuscript/fig4/plots/community_alternative.svg', fig = p)

###########################################################################################
### Supporting information: Simpsons paradox
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1_unequal', 'rs1_unequal', 'rs2_unequal', 'rs3_unequal')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots(sim_nms, colpal, symbolvals)
sjPlot::save_plot('./figures/main_manuscript/sup_unequal_habitats/plots/cooccurrence_rscdist.svg', fig = p)
###########################################################################################
### Supporting information: Co-occurrence plots for noise
colpal <- c('#000000FF', '#999999FF', '#CCCCCCFF', '#666666FF')
sim_nms <- c('rc1_noise', 'rs1_noise', 'rs2_noise', 'rs3_noise')
symbolvals <- c(0,2,6,5)
p <- coocurrence_plots_noise(sim_nms, noise_lvl = 'noise_0.05', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist.svg', fig = p)
p <- coocurrence_plots_noise(sim_nms, noise_lvl = 'noise_0.1', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist.svg', fig = p)
p <- coocurrence_plots_noise(sim_nms, noise_lvl = 'noise_0.25', colpal, symbolvals) 
sjPlot::save_plot('./figures/supporting_information/sup_volume_noise/plots/cooccurrence_rscdist.svg', fig = p)
###########################################################################################
### Supporting information: Species pairwise preference correlations (rs3)
subdir <- "./simulation_data/fig4/rs3/rs3"
p <- get_pref_effect(subdir)
sjPlot::save_plot('./figures/supporting_information/sup_pref_r1/plots/cooccurrence_rscdist.svg', fig = p)


###########################################################################################
##### Fig.5: Network plots                                                            #####
###########################################################################################


