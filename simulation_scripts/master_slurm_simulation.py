###############################################################################################
#### SAMSARA Master Script - Simulations                                                    ###
###############################################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### This script is intended to run in the conda environment, which we provide in 
### 'SAMSARA_repository/Conda_environment/' and will likely not work if not activated
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Load libraries & setwd for loading functions                                            ###
import numpy as np
import random
import networkx as nx
import datetime
import pandas as pd
import os
import getopt
import sys

##############################################################################################################################
##### 1. Load function scripts
##############################################################################################################################
sys.path.append(os.path.abspath('./simulation_scripts/simulation_functions'))
os.sys.path.append('./simulation_scripts/simulation_functions')
from set_species import *
from set_NW import *
from set_dynamics import *
from set_standard_initial_conditions import *

##############################################################################################################################
##### 2. Set working directory, define initial set of conditions and create output directories
##############################################################################################################################
## Unhash if inputs given in .bash script
#experiment = sys.argv[1]
#group = sys.argv[2]
#simulation = sys.argv[3]
#std_init_cond = sys.argv[4]

experiment = 'fig1'
treatment = 'r1'
simulation = 'r1'
std_init_cond = 'basic'

## If the parent directory does not exist yet it is generated here, otherwise just the subdirectory is generated here
os.makedirs(os.path.join('./simulation_data', experiment), exist_ok=True)
os.makedirs(os.path.join('./simulation_data', experiment, treatment), exist_ok=True)
os.makedirs(os.path.join('./simulation_data', experiment, treatment, simulation), exist_ok=True)
os.chdir(os.path.join('./simulation_data', experiment, treatment, simulation))

##############################################################################################################################
##### 3. Initialise parameters from bash input
##############################################################################################################################
### All functions called in this section can be found in './Load_simulation_functions/set_standard_initial_conditions.py'
# Here the standard initial conditions are called -> These are based on the parameter 'std_init_cond' defined above (here = 'basic')
n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions(std_init_cond)
## This next function was added to overwrite what is defined in the standard initial conditions by adding additional arguments to the .bash command
#n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_cmd_conditions(n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod)

## If you are interested in changing certain settings this is where you can do so manually:
nrep = 2            # We here only simulate 2 runs (the standard number of replicates is 100)
n_num = 25          # We only simulate 25 habitats (the standard number of habitats being 300)
# See the function 'set_initial_conditions_basic' in './Load_simulation_functions/set_standard_initial_conditions.py'
# to find a list of all parameters and what they mean

## Here we need to calculate some some additional parameters from our initial inputs
# What these individual parameters mean is described in the function 'set_additional_conditions'
tot_perc_int, perc_comp, perc_facil, evt_dur, sp_ini_num, seed_pop, neg_grw_thr, rsc_cube_gen_vals, t_dg, t_g, D = set_additional_conditions(sp_num,per_sp_int,ratio_cf,evt_num,evt_dur_var,sp_ini_perc,sp_ini_abd,n_num,kcoeff,neg_frac_kcoeff,rsc_cube_gen_vals,rsc_cube_gen_rule,run_time,sim_steps,D_set)


##############################################################################################################################
##### 4. Save log files of either inputs or simulated data
##############################################################################################################################
resource_cube_log = [[]]*nrep                   # A numpy data file that saves resource cubes data for the first run generated
sp_hab_carcap_log = pd.DataFrame()              # Species habitat carrying capacities(K_{i,h})
sp_grwrt_log = pd.DataFrame()                   # Species growth rates
n_knorm_log = pd.DataFrame()                    # Relative resource abundances per habitat
n_kabs_log = pd.DataFrame()                     # Absolute resource abundances per habitat
node_xyz_log = pd.DataFrame()                   # Coordinates of all habitats
sp_int_mat_log = pd.DataFrame()                 # Log of Interaction matrices
sp_rsc_prf_log = pd.DataFrame()                 # Log of all species resource preferences (theta_{i,k})
euc_rsc_prf_sim_mat_log = pd.DataFrame()        # Species euclidean distances in resource preferences matrices

cooc_cor_adj_mat_log = pd.DataFrame()           # Species correlation matrices
signif_cooc_cor_adj_mat_log = pd.DataFrame()    # Species co-ocurrence matrices (-> trimmed for correlation threshold and significance testing)

n_pop_log_log = pd.DataFrame()                  # Main output: Species abundances at final timepoint
n_pop_log_log_first = pd.DataFrame()            # Log file to check species abundances over time of first simulated run
# -> At 100 timepoints across the simulated period abundances are saved. This allows us to plot abundance curves over time to see whether
# given our settings the populations have equilibrated.

param_log = pd.DataFrame()                      # Log file of settings
# Simple outputs of names, indices and errors -> Not used outside this script
sp_indiv_data_log = pd.DataFrame()
hab_indiv_data_log = pd.DataFrame()
rsc_indiv_data_log = pd.DataFrame()


##############################################################################################################################
##### 5. Run simulations
##############################################################################################################################

##############################################################################################################################
### Run simulations for all n replicates (nrep)
#q = 0
for q in range(nrep):
    
    print("nrep = "+str(q))

    ##########################################################################################################################
    ### Initialise data for simulation
    ## See specific description of used functions and their inputs in respective functions scripts (given in brackets)
    # Generates resource cubes
    resource_cubes_x = get_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_rule, rsc_cube_gen_vals) #(function in './Load_simulation_functions/set_NW.py')
    # Sets species resource preference
    sp_rsc_prf = set_prefernce(sp_num, rsc_num, pref_type)                                              #(function in './Load_simulation_functions/set_species.py')
    # Sets species interaction matrix
    sp_int_mat = rdm_sym_sp_int_mat(sp_num, perc_comp, perc_facil, int_minmax, int_mod)                 #(function in './Load_simulation_functions/set_species.py')
    # Generates the habitat network
    G, pos, n_num, n_nms, node_xyz, edge_xyz, n_nbr_nms, n_nbr_idx = set_NW_random_geometric_graph(n_num, radius) #(function in './Load_simulation_functions/set_NW.py')
    # Sets habitats individual resource values
    n_knorm, n_kabs = set_node_resource_vals(n_num, node_xyz, shape, resource_cubes_x, rsc_num)         #(function in './Load_simulation_functions/set_NW.py')
    # Calculates species habitat carrying capacity
    sp_hab_carcap = set_sp_hab_carcapx(n_num,sp_rsc_prf,n_knorm,n_kabs,neg_grw_thr,sp_ini_abd,kcoeff)   #(function in './Load_simulation_functions/set_species.py') 
    # Modifies or overrides species habitat carrying capacity, f.ex. for runs with random K's 
    sp_hab_carcap = mod_sp_hab_carcap(sp_hab_carcap, kmod, kcoeff)                                      #(function in './Load_simulation_functions/set_species.py') 
    # Sets species growth rates based on their habitat carrying capacity (K < 0.1 -> r = 0 ; K < 0.1 -> r = 1)  
    sp_grwrt = set_sp_grwrt_bin(n_num,sp_hab_carcap, neg_grw_thr, killem = False)                       #(function in './Load_simulation_functions/set_species.py')
    # Get number of neighbors a habitat has in the habitat network
    n_nbr_num = [len(n_nbr_idx[i]) for i in range(len(n_nbr_idx))]
        
    ##########################################################################################################################
    ### Main simulation
    # Main simulation function
    n_pop_log, error_count = full_sim_no_error_stop(t_g, t_dg, evt_num, evt_dur, D, n_num, sp_int_mat, seed_pop, n_nbr_idx, sp_hab_carcap, sp_grwrt) #(function in './Load_simulation_functions/set_NW.py')
    # Get species, habitat and resource identifyers 'names' for saving
    sp_nms, rsc_nms, n_nms = get_nms(n_num,rsc_num,sp_num)                                              #(function in './Load_simulation_functions/set_NW.py')
    # Generate euclidean distance matrix for species resource preferences
    eudimat = eudi_rscprf(sp_rsc_prf)                                                                   #(function in './Load_simulation_functions/set_species.py')

    ##########################################################################################################################
    ### Speces, habitat and resource identifyers (+ replicate index and number of error runs) to be added to respective save data
    sp_xtra =  pd.DataFrame({"repl":np.repeat(q, sp_num).tolist(),"err":np.repeat(error_count, sp_num).tolist(),"sp_nms":sp_nms, "Ridx":np.arange(1,(sp_num+1),1)})
    hab_xtra = pd.DataFrame({"repl":np.repeat(q, n_num).tolist(),"err":np.repeat(error_count, n_num).tolist(),"n_nms":n_nms, "Ridx":np.arange(1,(n_num+1),1)})        
    rsc_xtra = pd.DataFrame({"repl":np.repeat(q, rsc_num).tolist(),"err":np.repeat(error_count, rsc_num).tolist(), "rsc_nms":rsc_nms, "Ridx":np.arange(1,(rsc_num+1),1)})    
    ### Save data (-> Find content of each data set in 4.)
    mean_sp_int = np.mean(sp_int_mat)
    var_sp_int = np.var(sp_int_mat)
    param = pd.DataFrame([q,error_count,n_num,rsc_num,radius,rsc_cube_gen_rule,str(rsc_cube_gen_vals),sp_num,perc_comp,perc_facil,int_minmax,sp_ini_perc,D,var_sp_int,mean_sp_int,shape[0]])
    param_log = pd.concat([param_log,param.T])
    sp_int_mat_log = pd.concat([sp_int_mat_log,pd.concat([pd.DataFrame(sp_int_mat), sp_xtra], axis = 1)])
    euc_rsc_prf_sim_mat_log = pd.concat([euc_rsc_prf_sim_mat_log,pd.concat([pd.DataFrame(eudimat), sp_xtra], axis = 1)])
    sp_rsc_prf_log = pd.concat([sp_rsc_prf_log,pd.concat([pd.DataFrame(sp_rsc_prf), sp_xtra], axis = 1)])
    sp_hab_carcap_log = pd.concat([sp_hab_carcap_log,pd.concat([pd.DataFrame(sp_hab_carcap), hab_xtra], axis = 1)])
    sp_grwrt_log = pd.concat([sp_grwrt_log,pd.concat([pd.DataFrame(sp_grwrt), hab_xtra], axis = 1)])
    n_knorm_log = pd.concat([n_knorm_log,pd.concat([pd.DataFrame(n_knorm), hab_xtra], axis = 1)])
    n_kabs_log = pd.concat([n_kabs_log,pd.concat([pd.DataFrame(n_kabs), hab_xtra], axis = 1)])
    node_xyz_log = pd.concat([node_xyz_log,pd.concat([pd.DataFrame(node_xyz), hab_xtra], axis = 1)])
    sp_indiv_data_log = pd.concat([sp_indiv_data_log,sp_xtra])
    hab_indiv_data_log = pd.concat([hab_indiv_data_log,pd.concat([pd.DataFrame(n_nbr_num), hab_xtra], axis = 1)])
    rsc_indiv_data_log = pd.concat([rsc_indiv_data_log,rsc_xtra])
    n_pop_log_log = pd.concat([n_pop_log_log,pd.concat([pd.DataFrame(n_pop_log[-1]), hab_xtra], axis = 1)])
    ### In the first run both generated resource landscape data is saved, as well as species abundance over 100 timepoints is saved
    if q == 0:
        savesteps = np.round(np.linspace(0, evt_dur, num = 101))
        n_pop_log_log_first = pd.DataFrame(n_pop_log[:,0,:])
        savesteps = np.append(savesteps,np.linspace(evt_dur+1, evt_dur+3, num = 3))
        n_pop_log_log_first = pd.concat([n_pop_log_log_first.reset_index(drop=True), pd.DataFrame(savesteps)], axis=1)
        copyspnms = np.copy(sp_nms).tolist()
        copyspnms.append('savesteps')
        n_pop_log_log_first.columns = copyspnms
        n_pop_log_log_first.to_csv("n_pop_log_log_first.csv")
        np.save('resource_cube_log_first', resource_cubes_x)


##############################################################################################################################
##### 6. Save data
##############################################################################################################################

##############################################################################################################################
### Export data to .csv for further analysis in .R (open './Slurm_stuff/master_slurm_analysis.R')

## Settings and parameters saved in every single run (not used in further analysis)
param_log.columns = ["nrep","errors","n_num","rsc_num","radius","rsc_cube_gen_rule","rsc_cube_gen_vals","sp_num","perc_comp","perc_facil","int_minmax","sp_ini_perc", "D_set","var_sp_int","mean_sp_int","shape"]
param_log.to_csv("param_log.csv")
#param = pd.DataFrame([q,error_count,n_num,rsc_num,radius,rsc_cube_gen_rule,rsc_cube_gen_vals,sp_num,perc_comp,perc_facil,int_minmax,sp_ini_perc,D_set,var_sp_int,mean_sp_int,shape[0]])

## Settings and parameters: One values for all runs
# This is used! But keep in mind that this only saves the parameter values at this final point, so none that vary between runs
settings_log = pd.DataFrame([['experiment',experiment],
                            ['treatment',treatment],
                            ['simulation',simulation],
                            ['std_init_cond',std_init_cond],
                            ['n_num',n_num],
                            ['rsc_num',rsc_num],
                            ['rsc_cube_gen_rule',rsc_cube_gen_rule], 
                            ['rsc_cube_gen_vals',rsc_cube_gen_vals], 
                            ['shape',shape[0]],
                            ['radius',radius],
                            ['kcoeff',kcoeff],
                            ['neg_frac_kcoeff',neg_frac_kcoeff],
                            ['kmod',kmod],
                            ['squared_rsc',squared_rsc],
                            ['pref_type',pref_type],
                            ### Species
                            ['sp_num',sp_num],
                            ['per_sp_int',per_sp_int],
                            ['ratio_cf',ratio_cf],
                            ['int_minmax',int_minmax],
                            ### Simulation
                            ['D_set',D_set],
                            ['nrep',nrep],
                            ['t_g',t_g],
                            ['t_dg',t_dg],
                            ['evt_num',evt_num],
                            ['evt_dur',evt_dur],
                            ['evt_dur_var',evt_dur_var],
                            ['sp_ini_perc',sp_ini_perc],
                            ['sp_ini_abd',sp_ini_abd],
                            ### Calculated metrics
                            ['neg_grw_thr',neg_grw_thr],
                            ['tot_perc_int',tot_perc_int],
                            ['perc_comp',perc_comp],
                            ['perc_facil',perc_facil],
                            ['sp_ini_num',sp_ini_num]])
settings_log.columns = ["setting","value"]
settings_log.to_csv("settings_log.csv")

## All other output files
sp_int_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
sp_int_mat_log.to_csv("sp_int_mat_log.csv")
euc_rsc_prf_sim_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
euc_rsc_prf_sim_mat_log.to_csv("euc_rsc_prf_sim_mat_log.csv")
sp_rsc_prf_log.columns = rsc_nms + rsc_xtra.columns.tolist()
sp_rsc_prf_log.to_csv("sp_rsc_prf_log.csv")
sp_hab_carcap_log.columns = sp_nms + hab_xtra.columns.tolist()
sp_hab_carcap_log.to_csv("sp_hab_carcap_log.csv")
sp_grwrt_log.columns = sp_nms + hab_xtra.columns.tolist()
sp_grwrt_log.to_csv("sp_grwrt_log.csv")
n_knorm_log.columns = rsc_nms + hab_xtra.columns.tolist()
n_knorm_log.to_csv("n_knorm_log.csv")
n_kabs_log.columns = rsc_nms + hab_xtra.columns.tolist()
n_kabs_log.to_csv("n_kabs_log.csv")
node_xyz_log.columns = ['x_coord','y_coord','z_coord'] + hab_xtra.columns.tolist()
node_xyz_log.to_csv("node_xyz_log.csv")
## Main output: Populations per habitat
n_pop_log_log.columns = sp_nms + sp_xtra.columns.tolist()
n_pop_log_log.to_csv("n_pop_log_log.csv")


