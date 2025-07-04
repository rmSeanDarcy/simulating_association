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
experiment <- 'fig2_test'

###########################################################################################
##### Load output data
fd_hab <- read.csv(paste0('./analysis_data/', experiment,'/full_res_compiled.csv'))
pd_hab <- read.csv(paste0('./analysis_data/', experiment,'/infm_res_compiled.csv'))

###########################################################################################
### Plotting settings:
int_col <- '#ffa500ff'    #'#F1DAFFFF'
eucint_col <- '#008080ff'
euc_col <- '#008080ff'
rest_col <- '#D1D3D2FF'
#
multipl_q = (((10*10)-10)/2)*0.25
multipl_posint = 4
multipl_negint = 16
boxalpha = 0.8
boxcolor = 'white'
sml_font <- 8
med_font <- 10
big_font <- 11
titltxt <- 'hi'
hnum <- 25
plotdims = c(4,4)

###########################################################################################
##### Fig.2: Simple system                                                            #####
###########################################################################################

###########################################################################################
##### Fig.2 B,C) Matching co-occurrences                                              #####
p <- plot_match_allruns_eucintp(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintp',hnum, plotdims = c(4,4))
sjPlot::save_plot(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses,'/Fig2_B.png'), fig = p)
p <- plot_match_allruns_eucintn(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintn',hnum, plotdims = c(4,4))
sjPlot::save_plot(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses,'/Fig2_C.png'), fig = p)


###########################################################################################
##### Fig.3: Meta-community dynamics                                                  #####
###########################################################################################


###########################################################################################
##### Fig.4: Sampling volume                                                          #####
###########################################################################################


###########################################################################################
##### Fig.5: Network plots                                                            #####
###########################################################################################


