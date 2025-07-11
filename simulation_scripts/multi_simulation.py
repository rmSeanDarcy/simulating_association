###############################################################################################
#### SAMSARA Master Script - Simulations                                                    ###
###############################################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### This script is intended to run in the provided conda environment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Load necessary modules                                                                  ###
import numpy as np
import random
import networkx as nx
import datetime
import pandas as pd
import os
import getopt
import sys
import argparse

##############################################################################################################################
##### 1. Load custom function scripts
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
experiment = 'fig5'
#treatment = sys.argv[2]        # Specified when saving!
#simulation = sys.argv[3]
std_init_cond = 'basic'
#experiment = 'fig1x'
#treatment = 'r1'
#simulation = 'r1'
#std_init_cond = 'basic'

## If the parent directory does not exist yet it is generated here, otherwise just the subdirectory is generated here
os.makedirs(os.path.join('./simulation_data/xfig5'), exist_ok=True)

os.makedirs(os.path.join('./simulation_data/xfig5/nw_c1'), exist_ok=True)
os.makedirs(os.path.join('./simulation_data/xfig5/nw_c2'), exist_ok=True)
os.makedirs(os.path.join('./simulation_data/xfig5/nw_s1'), exist_ok=True)

os.makedirs(os.path.join('./simulation_data/xfig5/nw_c1/nw_c1'), exist_ok=True)
os.makedirs(os.path.join('./simulation_data/xfig5/nw_c2/nw_c2'), exist_ok=True)
os.makedirs(os.path.join('./simulation_data/xfig5/nw_s1/nw_s1'), exist_ok=True)

##############################################################################################################################
##### 3. Initialise parameters from bash input
##############################################################################################################################
### All functions called in this section can be found in './Load_simulation_functions/set_standard_initial_conditions.py'
# Here the standard initial conditions are called -> These are based on the parameter 'std_init_cond' defined above (here = 'basic')
n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions(std_init_cond)
## This next function was added to overwrite what is defined in the standard initial conditions by adding additional arguments to the .bash command
n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_cmd_conditions(n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod)

## If you are interested in changing certain settings this is where you can do so manually:
nrep = 2            # We here only simulate 2 runs (the standard number of replicates is 100)
#n_num = 25          # We only simulate 25 habitats (the standard number of habitats being 300)
# See the function 'set_initial_conditions_basic' in './Load_simulation_functions/set_standard_initial_conditions.py'
# to find a list of all parameters and what they mean

## Here we need to calculate some some additional parameters from our initial inputs
# What these individual parameters mean is described in the function 'set_additional_conditions'
tot_perc_int, perc_comp, perc_facil, evt_dur, sp_ini_num, seed_pop, neg_grw_thr, rsc_cube_gen_vals, t_dg, t_g, D = set_additional_conditions(sp_num,per_sp_int,ratio_cf,evt_num,evt_dur_var,sp_ini_perc,sp_ini_abd,n_num,kcoeff,neg_frac_kcoeff,rsc_cube_gen_vals,rsc_cube_gen_rule,run_time,sim_steps,D_set)


##############################################################################################################################
##### 4. Save log files of either inputs or simulated data
##############################################################################################################################

##############################################################################################################################
### Once for c1
c1_resource_cube_log = [[]]*nrep                   # A numpy data file that saves resource cubes data for the first run generated
c1_sp_hab_carcap_log = []              # Species habitat carrying capacities(K_{i,h})
c1_sp_grwrt_log = []                   # Species growth rates
c1_n_knorm_log = []                    # Relative resource abundances per habitat
c1_n_kabs_log = []                     # Absolute resource abundances per habitat
c1_node_xyz_log = []                   # Coordinates of all habitats
c1_sp_int_mat_log = []                 # Log of Interaction matrices
c1_sp_rsc_prf_log = []                 # Log of all species resource preferences (theta_{i,k})
c1_euc_rsc_prf_sim_mat_log = []        # Species euclidean distances in resource preferences matrices
c1_cooc_cor_adj_mat_log = []           # Species correlation matrices
c1_signif_cooc_cor_adj_mat_log = []    # Species co-ocurrence matrices (-> trimmed for correlation threshold and significance testing)
c1_n_pop_log_log = []                  # Main output: Species abundances at final timepoint
c1_n_pop_log_log_first = []            # Log file to check species abundances over time of first simulated run
c1_param_log = []                      # Log file of settings
c1_sp_indiv_data_log = []
c1_hab_indiv_data_log = []
c1_rsc_indiv_data_log = []

##############################################################################################################################
### Once for c2
c2_resource_cube_log = [[]]*nrep                   # A numpy data file that saves resource cubes data for the first run generated
c2_sp_hab_carcap_log = []              # Species habitat carrying capacities(K_{i,h})
c2_sp_grwrt_log = []                   # Species growth rates
c2_n_knorm_log = []                    # Relative resource abundances per habitat
c2_n_kabs_log = []                     # Absolute resource abundances per habitat
c2_node_xyz_log = []                   # Coordinates of all habitats
c2_sp_int_mat_log = []                 # Log of Interaction matrices
c2_sp_rsc_prf_log = []                 # Log of all species resource preferences (theta_{i,k})
c2_euc_rsc_prf_sim_mat_log = []        # Species euclidean distances in resource preferences matrices
c2_cooc_cor_adj_mat_log = []           # Species correlation matrices
c2_signif_cooc_cor_adj_mat_log = []    # Species co-ocurrence matrices (-> trimmed for correlation threshold and significance testing)
c2_n_pop_log_log = []                  # Main output: Species abundances at final timepoint
c2_n_pop_log_log_first = []            # Log file to check species abundances over time of first simulated run
c2_param_log = []                      # Log file of settings
c2_sp_indiv_data_log = []
c2_hab_indiv_data_log = []
c2_rsc_indiv_data_log = []

##############################################################################################################################
### Once for s1
s1_resource_cube_log = [[]]*nrep                   # A numpy data file that saves resource cubes data for the first run generated
s1_sp_hab_carcap_log = []              # Species habitat carrying capacities(K_{i,h})
s1_sp_grwrt_log = []                   # Species growth rates
s1_n_knorm_log = []                    # Relative resource abundances per habitat
s1_n_kabs_log = []                     # Absolute resource abundances per habitat
s1_node_xyz_log = []                   # Coordinates of all habitats
s1_sp_int_mat_log = []                 # Log of Interaction matrices
s1_sp_rsc_prf_log = []                 # Log of all species resource preferences (theta_{i,k})
s1_euc_rsc_prf_sim_mat_log = []        # Species euclidean distances in resource preferences matrices
s1_cooc_cor_adj_mat_log = []           # Species correlation matrices
s1_signif_cooc_cor_adj_mat_log = []    # Species co-ocurrence matrices (-> trimmed for correlation threshold and significance testing)
s1_n_pop_log_log = []                  # Main output: Species abundances at final timepoint
s1_n_pop_log_log_first = []            # Log file to check species abundances over time of first simulated run
s1_param_log = []                      # Log file of settings
s1_sp_indiv_data_log = []
s1_hab_indiv_data_log = []
s1_rsc_indiv_data_log = []

##############################################################################################################################
##### 5. Run simulations
##############################################################################################################################

### Backup variables when changing in different scenarios
Xperc_comp = np.copy(perc_comp)
Xperc_facil = np.copy(perc_facil)
Xkmod = np.copy(kmod)

##############################################################################################################################
### Run simulations for all n replicates (nrep)
#q = 0
for q in range(nrep):
    
    print("nrep = "+str(q))

    ##########################################################################################################################
    ### All of these parameters will be set once before all runs and overwritten within a run
    perc_comp = np.copy(Xperc_comp)
    perc_facil = np.copy(Xperc_facil)
    kmod = np.copy(Xkmod)
    sp_rsc_prf = set_prefernce(sp_num, rsc_num, pref_type)
    eudimat = eudi_rscprf(sp_rsc_prf)
    sp_int_mat = rdm_sym_sp_int_mat(sp_num, perc_comp, perc_facil, int_minmax, int_mod)
    G, pos, n_num, n_nms, node_xyz, edge_xyz, n_nbr_nms, n_nbr_idx = set_NW_random_geometric_graph(n_num, radius)
    mean_sp_int = np.mean(sp_int_mat)
    var_sp_int = np.var(sp_int_mat)
    resource_cubes_x = get_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_rule, rsc_cube_gen_vals)    
    n_knorm, n_kabs = set_node_resource_vals(n_num, node_xyz, shape, resource_cubes_x, rsc_num)
    sp_hab_carcap = set_sp_hab_carcapx(n_num,sp_rsc_prf,n_knorm,n_kabs,neg_grw_thr,sp_ini_abd,kcoeff)    
    sp_hab_carcap = mod_sp_hab_carcap(sp_hab_carcap, kmod, kcoeff)
    sp_grwrt = set_sp_grwrt_bin(n_num,sp_hab_carcap, neg_grw_thr, killem = False)        
    sp_nms, rsc_nms, n_nms = get_nms(n_num,rsc_num,sp_num)    
    ### Save these initial conditions!
    Xsp_int_mat = np.copy(sp_int_mat)
    Xsp_hab_carcap = np.copy(sp_hab_carcap)
    Xsp_grwrt = np.copy(sp_grwrt)
    
    ##############################################################################################################################
    ### Simulation c1: Interactions only
    ##############################################################################################################################
    # Load saved data
    sp_int_mat = np.copy(Xsp_int_mat)
    sp_hab_carcap = np.copy(Xsp_hab_carcap)
    sp_grwrt = np.copy(Xsp_grwrt)
    # Modify simulation settings!!!
    kmod = 'rdm_normal'                         #'keep' 'rdm_uniform' 'rdm_normal'
    sp_hab_carcap = mod_sp_hab_carcap(sp_hab_carcap, kmod, kcoeff)
    sp_grwrt = set_sp_grwrt_bin(n_num,sp_hab_carcap, neg_grw_thr, killem = False)
    # Run simulation and save numeric outputs of run into log pd.DataFrames
    n_pop_log, error_count = full_sim_no_error_stop(t_g, t_dg, evt_num, evt_dur, D, n_num, sp_int_mat, seed_pop, n_nbr_idx, sp_hab_carcap, sp_grwrt)
    ##########################################################################################################################
    ### Fill log data
    # Speces, habitat and resource identifyers (+ replicate index and number of error runs) to be added to respective save data
    sp_xtra =  pd.DataFrame({"repl":np.repeat(q, sp_num).tolist(),"err":np.repeat(error_count, sp_num).tolist(),"sp_nms":sp_nms, "Ridx":np.arange(1,(sp_num+1),1)})
    hab_xtra = pd.DataFrame({"repl":np.repeat(q, n_num).tolist(),"err":np.repeat(error_count, n_num).tolist(),"n_nms":n_nms, "Ridx":np.arange(1,(n_num+1),1)})        
    rsc_xtra = pd.DataFrame({"repl":np.repeat(q, rsc_num).tolist(),"err":np.repeat(error_count, rsc_num).tolist(), "rsc_nms":rsc_nms, "Ridx":np.arange(1,(rsc_num+1),1)})    
    # Save data (-> Find content of each data set in 4.)
    mean_sp_int = np.mean(sp_int_mat)
    var_sp_int = np.var(sp_int_mat)
    param = [q,error_count,n_num,rsc_num,radius,rsc_cube_gen_rule,str(rsc_cube_gen_vals),sp_num,perc_comp,perc_facil,int_minmax,sp_ini_perc,D,var_sp_int,mean_sp_int,shape[0]]
    c1_param_log.append(param)
    c1_sp_int_mat_log.append(pd.concat([pd.DataFrame(sp_int_mat), sp_xtra], axis = 1))
    c1_euc_rsc_prf_sim_mat_log.append(pd.concat([pd.DataFrame(eudimat), sp_xtra], axis = 1))
    c1_sp_rsc_prf_log.append(pd.concat([pd.DataFrame(sp_rsc_prf), sp_xtra], axis = 1))
    c1_sp_hab_carcap_log.append(pd.concat([pd.DataFrame(sp_hab_carcap), hab_xtra], axis = 1))
    c1_sp_grwrt_log.append(pd.concat([pd.DataFrame(sp_grwrt), hab_xtra], axis = 1))
    c1_n_knorm_log.append(pd.concat([pd.DataFrame(n_knorm), hab_xtra], axis = 1))
    c1_n_kabs_log.append(pd.concat([pd.DataFrame(n_kabs), hab_xtra], axis = 1))
    c1_node_xyz_log.append(pd.concat([pd.DataFrame(node_xyz), hab_xtra], axis = 1))
    c1_n_pop_log_log.append(pd.concat([pd.DataFrame(n_pop_log[-1]), hab_xtra], axis = 1))
    
    ##############################################################################################################################
    ### Simulation c1: Environment only
    ##############################################################################################################################
    # Load saved data
    sp_int_mat = np.copy(Xsp_int_mat)
    sp_hab_carcap = np.copy(Xsp_hab_carcap)
    sp_grwrt = np.copy(Xsp_grwrt)
    ### Modify simulation settings!!!
    perc_comp = 0
    perc_facil = 0
    sp_int_mat = rdm_sym_sp_int_mat(sp_num, perc_comp, perc_facil, int_minmax, int_mod)
    # Run simulation and save numeric outputs of run into log pd.DataFrames
    n_pop_log, error_count = full_sim_no_error_stop(t_g, t_dg, evt_num, evt_dur, D, n_num, sp_int_mat, seed_pop, n_nbr_idx, sp_hab_carcap, sp_grwrt)
    ##########################################################################################################################
    ### Fill log data
    # Speces, habitat and resource identifyers (+ replicate index and number of error runs) to be added to respective save data
    sp_xtra =  pd.DataFrame({"repl":np.repeat(q, sp_num).tolist(),"err":np.repeat(error_count, sp_num).tolist(),"sp_nms":sp_nms, "Ridx":np.arange(1,(sp_num+1),1)})
    hab_xtra = pd.DataFrame({"repl":np.repeat(q, n_num).tolist(),"err":np.repeat(error_count, n_num).tolist(),"n_nms":n_nms, "Ridx":np.arange(1,(n_num+1),1)})        
    rsc_xtra = pd.DataFrame({"repl":np.repeat(q, rsc_num).tolist(),"err":np.repeat(error_count, rsc_num).tolist(), "rsc_nms":rsc_nms, "Ridx":np.arange(1,(rsc_num+1),1)})    
    # Save data (-> Find content of each data set in 4.)
    mean_sp_int = np.mean(sp_int_mat)
    var_sp_int = np.var(sp_int_mat)
    param = [q,error_count,n_num,rsc_num,radius,rsc_cube_gen_rule,str(rsc_cube_gen_vals),sp_num,perc_comp,perc_facil,int_minmax,sp_ini_perc,D,var_sp_int,mean_sp_int,shape[0]]
    c2_param_log.append(param)
    c2_sp_int_mat_log.append(pd.concat([pd.DataFrame(sp_int_mat), sp_xtra], axis = 1))
    c2_euc_rsc_prf_sim_mat_log.append(pd.concat([pd.DataFrame(eudimat), sp_xtra], axis = 1))
    c2_sp_rsc_prf_log.append(pd.concat([pd.DataFrame(sp_rsc_prf), sp_xtra], axis = 1))
    c2_sp_hab_carcap_log.append(pd.concat([pd.DataFrame(sp_hab_carcap), hab_xtra], axis = 1))
    c2_sp_grwrt_log.append(pd.concat([pd.DataFrame(sp_grwrt), hab_xtra], axis = 1))
    c2_n_knorm_log.append(pd.concat([pd.DataFrame(n_knorm), hab_xtra], axis = 1))
    c2_n_kabs_log.append(pd.concat([pd.DataFrame(n_kabs), hab_xtra], axis = 1))
    c2_node_xyz_log.append(pd.concat([pd.DataFrame(node_xyz), hab_xtra], axis = 1))
    c2_n_pop_log_log.append(pd.concat([pd.DataFrame(n_pop_log[-1]), hab_xtra], axis = 1))

    ##############################################################################################################################
    ### Simulation s1: Both drivers
    ##############################################################################################################################
    # Load saved data -> 'basic' simulation parameters
    sp_int_mat = np.copy(Xsp_int_mat)
    sp_hab_carcap = np.copy(Xsp_hab_carcap)
    sp_grwrt = np.copy(Xsp_grwrt)
    # Run simulation and save numeric outputs of run into log pd.DataFrames
    n_pop_log, error_count = full_sim_no_error_stop(t_g, t_dg, evt_num, evt_dur, D, n_num, sp_int_mat, seed_pop, n_nbr_idx, sp_hab_carcap, sp_grwrt)
    ##########################################################################################################################
    ### Fill log data
    # Speces, habitat and resource identifyers (+ replicate index and number of error runs) to be added to respective save data
    sp_xtra =  pd.DataFrame({"repl":np.repeat(q, sp_num).tolist(),"err":np.repeat(error_count, sp_num).tolist(),"sp_nms":sp_nms, "Ridx":np.arange(1,(sp_num+1),1)})
    hab_xtra = pd.DataFrame({"repl":np.repeat(q, n_num).tolist(),"err":np.repeat(error_count, n_num).tolist(),"n_nms":n_nms, "Ridx":np.arange(1,(n_num+1),1)})        
    rsc_xtra = pd.DataFrame({"repl":np.repeat(q, rsc_num).tolist(),"err":np.repeat(error_count, rsc_num).tolist(), "rsc_nms":rsc_nms, "Ridx":np.arange(1,(rsc_num+1),1)})    
    # Save data (-> Find content of each data set in 4.)
    mean_sp_int = np.mean(sp_int_mat)
    var_sp_int = np.var(sp_int_mat)
    param = [q,error_count,n_num,rsc_num,radius,rsc_cube_gen_rule,str(rsc_cube_gen_vals),sp_num,perc_comp,perc_facil,int_minmax,sp_ini_perc,D,var_sp_int,mean_sp_int,shape[0]]
    s1_param_log.append(param)
    s1_sp_int_mat_log.append(pd.concat([pd.DataFrame(sp_int_mat), sp_xtra], axis = 1))
    s1_euc_rsc_prf_sim_mat_log.append(pd.concat([pd.DataFrame(eudimat), sp_xtra], axis = 1))
    s1_sp_rsc_prf_log.append(pd.concat([pd.DataFrame(sp_rsc_prf), sp_xtra], axis = 1))
    s1_sp_hab_carcap_log.append(pd.concat([pd.DataFrame(sp_hab_carcap), hab_xtra], axis = 1))
    s1_sp_grwrt_log.append(pd.concat([pd.DataFrame(sp_grwrt), hab_xtra], axis = 1))
    s1_n_knorm_log.append(pd.concat([pd.DataFrame(n_knorm), hab_xtra], axis = 1))
    s1_n_kabs_log.append(pd.concat([pd.DataFrame(n_kabs), hab_xtra], axis = 1))
    s1_node_xyz_log.append(pd.concat([pd.DataFrame(node_xyz), hab_xtra], axis = 1))
    s1_n_pop_log_log.append(pd.concat([pd.DataFrame(n_pop_log[-1]), hab_xtra], axis = 1))
    

##############################################################################################################################
##### 6. Save data
##############################################################################################################################

##############################################################################################################################
### c1 
headers = ["nrep","errors","n_num","rsc_num","radius","rsc_cube_gen_rule","rsc_cube_gen_vals","sp_num","perc_comp","perc_facil","int_minmax","sp_ini_perc", "D_set","var_sp_int","mean_sp_int","shape"]
c1_param_log = pd.DataFrame(c1_param_log, columns = headers)
c1_param_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/param_log.csv")
## All other output files
c1_sp_int_mat_log = pd.concat(c1_sp_int_mat_log)
c1_sp_int_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
c1_sp_int_mat_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/sp_int_mat_log.csv")

c1_euc_rsc_prf_sim_mat_log = pd.concat(c1_euc_rsc_prf_sim_mat_log)
c1_euc_rsc_prf_sim_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
c1_euc_rsc_prf_sim_mat_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/euc_rsc_prf_sim_mat_log.csv")

c1_sp_rsc_prf_log = pd.concat(c1_sp_rsc_prf_log)
c1_sp_rsc_prf_log.columns = rsc_nms + rsc_xtra.columns.tolist()
c1_sp_rsc_prf_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/sp_rsc_prf_log.csv")

c1_sp_hab_carcap_log = pd.concat(c1_sp_hab_carcap_log)
c1_sp_hab_carcap_log.columns = sp_nms + hab_xtra.columns.tolist()
c1_sp_hab_carcap_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/sp_hab_carcap_log.csv")

c1_sp_grwrt_log = pd.concat(c1_sp_grwrt_log)
c1_sp_grwrt_log.columns = sp_nms + hab_xtra.columns.tolist()
c1_sp_grwrt_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/sp_grwrt_log.csv")

c1_n_knorm_log = pd.concat(c1_n_knorm_log)
c1_n_knorm_log.columns = rsc_nms + hab_xtra.columns.tolist()
c1_n_knorm_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/n_knorm_log.csv")

c1_n_kabs_log = pd.concat(c1_n_kabs_log)
c1_n_kabs_log.columns = rsc_nms + hab_xtra.columns.tolist()
c1_n_kabs_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/n_kabs_log.csv")

c1_node_xyz_log = pd.concat(c1_node_xyz_log)
c1_node_xyz_log.columns = ['x_coord','y_coord','z_coord'] + hab_xtra.columns.tolist()
c1_node_xyz_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/node_xyz_log.csv")

## Main output: Populations per habitat
c1_n_pop_log_log = pd.concat(c1_n_pop_log_log)
c1_n_pop_log_log.columns = sp_nms + sp_xtra.columns.tolist()
c1_n_pop_log_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/n_pop_log_log.csv")

##############################################################################################################################
### c2 
headers = ["nrep","errors","n_num","rsc_num","radius","rsc_cube_gen_rule","rsc_cube_gen_vals","sp_num","perc_comp","perc_facil","int_minmax","sp_ini_perc", "D_set","var_sp_int","mean_sp_int","shape"]
c2_param_log = pd.DataFrame(c2_param_log, columns = headers)
c2_param_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/param_log.csv")
## All other output files
c2_sp_int_mat_log = pd.concat(c2_sp_int_mat_log)
c2_sp_int_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
c2_sp_int_mat_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/sp_int_mat_log.csv")

c2_euc_rsc_prf_sim_mat_log = pd.concat(c2_euc_rsc_prf_sim_mat_log)
c2_euc_rsc_prf_sim_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
c2_euc_rsc_prf_sim_mat_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/euc_rsc_prf_sim_mat_log.csv")

c2_sp_rsc_prf_log = pd.concat(c2_sp_rsc_prf_log)
c2_sp_rsc_prf_log.columns = rsc_nms + rsc_xtra.columns.tolist()
c2_sp_rsc_prf_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/sp_rsc_prf_log.csv")

c2_sp_hab_carcap_log = pd.concat(c2_sp_hab_carcap_log)
c2_sp_hab_carcap_log.columns = sp_nms + hab_xtra.columns.tolist()
c2_sp_hab_carcap_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/sp_hab_carcap_log.csv")

c2_sp_grwrt_log = pd.concat(c2_sp_grwrt_log)
c2_sp_grwrt_log.columns = sp_nms + hab_xtra.columns.tolist()
c2_sp_grwrt_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/sp_grwrt_log.csv")

c2_n_knorm_log = pd.concat(c2_n_knorm_log)
c2_n_knorm_log.columns = rsc_nms + hab_xtra.columns.tolist()
c2_n_knorm_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/n_knorm_log.csv")

c2_n_kabs_log = pd.concat(c2_n_kabs_log)
c2_n_kabs_log.columns = rsc_nms + hab_xtra.columns.tolist()
c2_n_kabs_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/n_kabs_log.csv")

c2_node_xyz_log = pd.concat(c2_node_xyz_log)
c2_node_xyz_log.columns = ['x_coord','y_coord','z_coord'] + hab_xtra.columns.tolist()
c2_node_xyz_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/node_xyz_log.csv")

## Main output: Populations per habitat
c2_n_pop_log_log = pd.concat(c2_n_pop_log_log)
c2_n_pop_log_log.columns = sp_nms + sp_xtra.columns.tolist()
c2_n_pop_log_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/n_pop_log_log.csv")


##############################################################################################################################
### s1 
headers = ["nrep","errors","n_num","rsc_num","radius","rsc_cube_gen_rule","rsc_cube_gen_vals","sp_num","perc_comp","perc_facil","int_minmax","sp_ini_perc", "D_set","var_sp_int","mean_sp_int","shape"]
s1_param_log = pd.DataFrame(s1_param_log, columns = headers)
s1_param_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/param_log.csv")
## All other output files
s1_sp_int_mat_log = pd.concat(s1_sp_int_mat_log)
s1_sp_int_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
s1_sp_int_mat_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/sp_int_mat_log.csv")

s1_euc_rsc_prf_sim_mat_log = pd.concat(s1_euc_rsc_prf_sim_mat_log)
s1_euc_rsc_prf_sim_mat_log.columns = sp_nms + sp_xtra.columns.tolist()
s1_euc_rsc_prf_sim_mat_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/euc_rsc_prf_sim_mat_log.csv")

s1_sp_rsc_prf_log = pd.concat(s1_sp_rsc_prf_log)
s1_sp_rsc_prf_log.columns = rsc_nms + rsc_xtra.columns.tolist()
s1_sp_rsc_prf_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/sp_rsc_prf_log.csv")

s1_sp_hab_carcap_log = pd.concat(s1_sp_hab_carcap_log)
s1_sp_hab_carcap_log.columns = sp_nms + hab_xtra.columns.tolist()
s1_sp_hab_carcap_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/sp_hab_carcap_log.csv")

s1_sp_grwrt_log = pd.concat(s1_sp_grwrt_log)
s1_sp_grwrt_log.columns = sp_nms + hab_xtra.columns.tolist()
s1_sp_grwrt_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/sp_grwrt_log.csv")

s1_n_knorm_log = pd.concat(s1_n_knorm_log)
s1_n_knorm_log.columns = rsc_nms + hab_xtra.columns.tolist()
s1_n_knorm_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/n_knorm_log.csv")

s1_n_kabs_log = pd.concat(s1_n_kabs_log)
s1_n_kabs_log.columns = rsc_nms + hab_xtra.columns.tolist()
s1_n_kabs_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/n_kabs_log.csv")

s1_node_xyz_log = pd.concat(s1_node_xyz_log)
s1_node_xyz_log.columns = ['x_coord','y_coord','z_coord'] + hab_xtra.columns.tolist()
s1_node_xyz_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/node_xyz_log.csv")

## Main output: Populations per habitat
s1_n_pop_log_log = pd.concat(s1_n_pop_log_log)
s1_n_pop_log_log.columns = sp_nms + sp_xtra.columns.tolist()
s1_n_pop_log_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/n_pop_log_log.csv")

##############################################################################################################################
### While slightly different in cases of kmod or so, we just export a single file for all of the simulations
## Settings and parameters: One values for all runs
# This is used! But keep in mind that this only saves the parameter values at this final point, so none that vary between runs
c1_settings_log = pd.DataFrame([['experiment',experiment],
                            ['treatment','nw_c1'],
                            ['simulation','nw_c1'],
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
c1_settings_log.columns = ["setting","value"]
c1_settings_log.to_csv("./simulation_data/xfig5/nw_c1/nw_c1/settings_log.csv")

##############################################################################################################################
### While slightly different in cases of kmod or so, we just export a single file for all of the simulations
## Settings and parameters: One values for all runs
# This is used! But keep in mind that this only saves the parameter values at this final point, so none that vary between runs
c2_settings_log = pd.DataFrame([['experiment',experiment],
                            ['treatment','nw_c2'],
                            ['simulation','nw_c2'],
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
c2_settings_log.columns = ["setting","value"]
c2_settings_log.to_csv("./simulation_data/xfig5/nw_c2/nw_c2/settings_log.csv")

##############################################################################################################################
### While only slightly different we need to specify the experiment, treatment and simulation label
s1_settings_log = pd.DataFrame([['experiment',experiment],
                            ['treatment','nw_s1'],
                            ['simulation','nw_s1'],
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
s1_settings_log.columns = ["setting","value"]
s1_settings_log.to_csv("./simulation_data/xfig5/nw_s1/nw_s1/settings_log.csv")



