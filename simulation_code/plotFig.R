###########################################################################################
### SAMSARA - Figure 2 and plots for related Supporting information                     ###
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
#install.packages("rlang", dependencies = TRUE)

### Initiate inputs
args <- commandArgs(trailingOnly=TRUE)
workdir <- args[1]
parent_set_of_analyses <- args[2]
### For testing
#workdir <- '/home/swani/Documents/computational_research_tools/homework4/SAMSARA'
#parent_set_of_analyses <- 'Fig2'


###########################################################################################
### Here we load all of the data analysis functions that are contained in scripts within the folder './Load_collectfunctions'
setwd(workdir)
for (i in list.files(path= paste(workdir,c('/simulation_code/Load_plot_functions'),sep=''), full.names=TRUE)) {
  source(i)
}

###########################################################################################
##### Load output data
pd_hab <- read.csv(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses,'/exp_results.csv'))

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
##### Fig.2: Main manuscript                                                          #####
###########################################################################################

###########################################################################################
##### Fig.2 B,C) Matching co-occurrences                                              #####
p <- plot_match_allruns_eucintp(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintp',hnum, plotdims = c(4,4))
sjPlot::save_plot(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses,'/Fig2_B.svg'), fig = p)
p <- plot_match_allruns_eucintn(pd_hab,label_mode = 'none',titltxt = 'plot_match_allruns_eucintn',hnum, plotdims = c(4,4))
sjPlot::save_plot(paste0(workdir,'/Result_master_dir/',parent_set_of_analyses,'/Fig2_C.svg'), fig = p)





