###############################################################################################
#### SAMSARA Master Script - Fixed initial conditions                                       ###
###############################################################################################
import getopt
import sys
import numpy as np
import random
import os


###############################################################################################
### Defining which set of basic initial conditions a script will be working with            ###
###############################################################################################

###############################################################################################
### Combining both interaction and resource effects
## Call 'basic' in function 'set_initial_conditions(std_init_cond = basic)'
def set_initial_conditions_basic():

    ### Basic environment settings:
    n_num = 300                     		# Number of habitats
    rsc_num = 3             				# Number of resources
    radius = 0              		  		# Connectivity radius d_{e} for habitats to become connected
    shape = (100, 100, 100)             	# Dimension of resource landscapes -> Gives number of voxels that are generated 

    ### Settings for resource distribution scenarios
    ## 1. Simple input  
    # 'fix_same' -> All three (depending on rsc_num) resource landscapes have the same lambda input (prevalence of small scale spatial structures, see 'Fyeld_generator' documentation)
    rsc_cube_gen_rule = 'fix_same'          # 'fix_same', 'specific_dist'
    # If rsc_cube_gen_rule = 'fix_same', then the 'rsc_cube_gen_vals' parameter can only take a single value as input
    rsc_cube_gen_vals = '0'                 # [5], [10] (or any other value)
    # If rsc_cube_gen_rule = 'fix_same', then this 'squared_rsc' parameter can be set to either take the basic outputs from 'Fyeld_generator' ('normal') or to subsequently square the values ('squared')
    squared_rsc = 'squared'                 # 'squared', 'normal'

    ## 2. Manual input
    # 'specific_dist' -> Tells the function that every resource landscape will be generated individually 
    #rsc_cube_gen_rule = 'specific_dist'         # 'fix_same', 'specific_dist'
    # If rsc_cube_gen_rule = 'specific_dist', then the rsc_cube_gen_vals contains a vector of strings for every number of resources (rsc_num).
    #rsc_cube_gen_vals = 'pwl_spcXsqrX0,pwl_spcXsqrX5,pwl_spcXsqrX9.99'   # 'pwl_spcXsqrX10','gradXsqrXtd','gradXnmlXlr'
    
    # The logic of these strings goes as follows:    
    # A. The first part, so before the first 'X' of the string, contains an identifyer what kind of resource landscape should be generated
    # pwl_spc 	-> landscape via Fyeld_generator
    # grad 	-> landscape via simple gradient
    # B. The second part, between the 'X's sets whether outputs should be squared or normal values (goes for both the gradient as well as the Fyeld_generator outpurs)
    # C. The last part, after the second 'X' depends on the first input
    # If the first input was 'pwl_spc' then this gives the lambda input
    # If the first input was 'grad' then here one can define from which side to which a gradient should be drawn (f.ex. 'td' = top -> down, or 'lr' = left -> right)
    
    ### Species settings:
    sp_num = 10                             # Number of species
    per_sp_int = 2.25 	         			# Number of interactions per species 
    # -> This is not a percentage. So if there are 10 species, there are 45 potential interactions. If 50 percent of these interactions should be realsied then there are 10*2.25 (=22.5) interactions  
    ratio_cf = 0.8 		               		# The ratio of competition to facilitation in the realised interactions
    int_minmax = 0.5    	      			# This gives the maximal interaction strength (not altered in the manuscript), minimum values is 0.1 
    int_mod = 'rdm_unif'                    # Mode for which interaction coefficients are generated (not altered in the manuscript) -> Random uniform value between 0.1 and 0.5 (or -0.1 and -0.5) 
    
    kmod = 'keep'                           # Modifies K for the control scenario c1 (interactions only) 
    # If 'rdm_uniform' or 'rdm_normal' -> Random K's are drawn (from uniform distribution or normal distribution); If 'equal' -> Keep generated K's 
    
    kcoeff = 1 				             	# Value for modifying species K's (not altered in manuscript) -> Keep fixed at 1!	
    neg_frac_kcoeff = 0.1 	         		# Threshold for K beneath which, species growth rates are set to zero
    # Multiplied with kcoeff in 'set_additional_conditions' function below, but as kcoeff fixed at 1 -> neg_frac_kcoeff = threshold
    
    pref_type = 'sequential'                # 'sequential', 'independent'    
    # Method of assigning species preference. Only 'sequential' is used in paper, 'independent' favours more generalist preferences
    
    ### Simulation parameters:
    evt_num = 1                             # Allows time periods of dispersal ('events') and periods without. Fixed at 1 in all analyses in paper! 
    evt_dur_var = 100,100                   # Related to setting above and gives duration and variability of events. Fixed at 100,100!
    run_time = 100                          # This is the total runtime of simulation
    sim_steps = 5000                        # These are the simulation step sizes input 
    # See 'set_additional_conditions' function below for how this is further modified as an input for the integrator
    
    sp_ini_perc = 1                         # This gives the fraction of species that are randomly selected to colonise a habitat
    # In scenario s5 sp_ini_perc = 0.8
    sp_ini_abd = 0.001                      # The abundance of species when starting simulations

    D_set = 0.0                             # The diffusion coefficient
    # This is further modified in the 'set_additional_conditions' function below as the strength of the diffusion also depends on time steps

    nrep = 100                              # Number of simulations or replicates
    
    ### Deprecated and not being used in functions anymore, so set to any value
    neg_grw_thr = 1                         # Now re-calculated in 'set_additional_conditions', still requires input here, but will be overwritten
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod 

###############################################################################################
### Control c1: Null condition for interactions -> No resource based stuff active -> For testing interaction based assembly
## Call 'null_int' in function 'set_initial_conditions(std_init_cond = null_int)'
def set_initial_conditions_null_int():
    ### Basic environment settings:
    n_num = 300                     	
    rsc_num = 3 			
    radius = 0  			
    shape = (100, 100, 100)             
    rsc_cube_gen_rule = 'fix_same'            
    rsc_cube_gen_vals = '0'              
    squared_rsc = 'squared'                 
    ### Species settings:
    sp_num = 10
    per_sp_int = 2.25
    ratio_cf = 0.8
    int_minmax = 0.5
    int_mod = 'rdm_unif'                     
    kcoeff = 1
    neg_frac_kcoeff = 0.1
    kmod = 'rdm_normal'                     
    neg_grw_thr = 1
    pref_type = 'sequential'                 
    ### Simulation parameters:
    evt_num = 1                     
    evt_dur_var = 100,100   
    run_time = 100
    sim_steps = 5000
    sp_ini_perc = 1   
    sp_ini_abd = 0.001                
    D_set = 0.0
    ### Number of simulations
    nrep = 100
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod 


###############################################################################################
### Control c2: Null conditions for resource based assembly -> No interactions are active
## Call 'null_rsc' in function 'set_initial_conditions(std_init_cond = null_rsc)'
def set_initial_conditions_null_rsc():
    ### Habitat network settings:
    n_num = 300                    
    ### Resource distribution settings:
    rsc_num = 3
    radius = 0
    shape = (100, 100, 100)            
    squared_rsc = 'squared'                   
    rsc_cube_gen_rule = 'fix_same'       
    rsc_cube_gen_vals = '0'              
    ### Species settings:
    sp_num = 10
    per_sp_int = 0
    ratio_cf = 0.8
    int_minmax = 0.
    int_mod = 'rdm_unif'                     
    kcoeff = 1
    neg_frac_kcoeff = 0.1
    kmod = 'keep'                          
    neg_grw_thr = 1
    pref_type = 'sequential'                  
    ### Simulation parameters:
    evt_num = 0                     
    evt_dur_var = 100,100           
    run_time = 100
    sim_steps = 5000
    sp_ini_perc = 1   
    sp_ini_abd = 0.001                
    D_set = 0.0
    ### Number of simulations
    nrep = 100
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod 


###############################################################################################
### Adding metacommunity dynamics
###############################################################################################
### Settings used for figure 3 
# -> Sets default connectivity radius d_e{e} (radius = 0.2) and diffusion coefficient (D_set = 0.1)
# We also setz the 'sim_steps' argument to 10000! This is because more integration steps improve the accurracy of solved differential equation (here via Runge Kutta method)
## Call 'metacom' in function 'set_initial_conditions(std_init_cond = metacom)'
def set_initial_conditions_metacom():
    ### Habitat network settings:
    n_num = 300                    
    ### Resource distribution settings:
    rsc_num = 3
    radius = 0.2
    shape = (100, 100, 100)            
    squared_rsc = 'squared'                      
    rsc_cube_gen_rule = 'fix_same'            
    rsc_cube_gen_vals = '0'                    
    ### Species settings:
    sp_num = 10
    per_sp_int = 2.25
    ratio_cf = 0.8
    int_minmax = 0.5
    int_mod = 'rdm_unif'                  
    kcoeff = 1
    neg_frac_kcoeff = 0.1
    kmod = 'keep'                           
    neg_grw_thr = 1
    pref_type = 'sequential'                
    ### Simulation parameters:
    evt_num = 1                     
    evt_dur_var = 100,100           
    run_time = 100
    sim_steps = 10000
    sp_ini_perc = 1   
    sp_ini_abd = 0.001                
    D_set = 0.1
    ### Number of simulations
    nrep = 100
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod 


###############################################################################################
### Calling standard initial conditions based on bash input                                 ###
###############################################################################################
### Main function that sets the standard initial conditions
def set_initial_conditions(std_init_cond):
    if std_init_cond == 'null_int':
        n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions_null_int()
    if std_init_cond == 'null_rsc':
        n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions_null_rsc()        
    elif std_init_cond == 'basic':
        n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions_basic()
    elif std_init_cond == 'metacom':
        n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod = set_initial_conditions_metacom()
    else:
        print('Error: Input for standard initial conditions incorrect')
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod


###############################################################################################
### Change initial conditions via bash script                                               ###
###############################################################################################
### Function that overrides the standard initial conditions if passed to .bash script
## In a bash script calling 'master_slurm_simulation.py' the first four additional arguments give 
# 1. working directory
# 2 and 3. Parent and subfolder name
# 4. Setting for standard initial conditions (see above)
# 5+ taken by the following function -> Add f.ex. 'python3 1. 2. 3. 4. sp_num=20 rsc_num=5 ...' to change 
def set_cmd_conditions(n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod): #pwl_spc,
    argv = sys.argv[5:]
    try:
        short_option = "nn:rn:sn:psi:rcf:imm:en:rtm:stp:sip:sia:d:ngt:rd:nr:wd:kc:nfkc:sr:pt:km:rgr:rgv:imd"   #ps:
        long_option = ["n_num=","rsc_num=","sp_num=","per_sp_int=","ratio_cf=","int_minmax=","evt_num=",  #"pwl_spc=",
                       "run_time=","sim_steps=","sp_ini_perc=","sp_ini_abd=","D_set=","neg_grw_thr=","radius=","nrep=", "kcoeff=", 
                       "neg_frac_kcoeff=", "squared_rsc=", "pref_type=", 'kmod=', 'rsc_cube_gen_rule=', 'rsc_cube_gen_vals=','int_mod=']
        opts,args = getopt.getopt(argv, short_option, long_option)
    except getopt.GetoptError as err:
        print(err)
        opts = []        
    for opt, arg in opts:
        if opt in ["-nn","--n_num"]:
            n_num = int(arg)
        elif opt in ["-rn","--rsc_num"]:
            rsc_num = int(arg)
        elif opt in ["-sn","--sp_num"]:
            sp_num = int(arg)
        elif opt in ["-psi","--per_sp_int"]:
            per_sp_int = float(arg)
        elif opt in ["-rcf","--ratio_cf"]:
            ratio_cf = float(arg)
        elif opt in ["-imm","--int_minmax"]:
            int_minmax = float(arg)
        elif opt in ["-en","--evt_num"]:
            evt_num = int(arg)
        elif opt in ["-rtm","--run_time"]:
            run_time = int(arg)
        elif opt in ["-stp","--sim_steps"]:
            sim_steps = int(arg)
        elif opt in ["-sip","--sp_ini_perc"]:
            sp_ini_perc = float(arg)
        elif opt in ["-sia","--sp_ini_abd"]:
            sp_ini_abd = float(arg)
        elif opt in ["-d","--D_set"]:
            D_set = float(arg)
        elif opt in ["-ngt","--neg_grw_thr"]:
            neg_grw_thr = float(arg)
        elif opt in ["-rd","--radius"]:
            radius = float(arg)
        elif opt in ["-nr","--nrep"]:
            nrep = int(arg)
        elif opt in ["-kc","--kcoeff"]:
            kcoeff = float(arg)
        elif opt in ["-nfkc","--neg_frac_kcoeff"]:
            neg_frac_kcoeff = float(arg)
        elif opt in ["-sr","--squared_rsc"]:
            squared_rsc = str(arg)
        elif opt in ["-pt","--pref_type"]:
            pref_type = str(arg)
        elif opt in ["-km","--kmod"]:
            kmod = str(arg)            
        elif opt in ["-rgr","--rsc_cube_gen_rule"]:
            rsc_cube_gen_rule = str(arg)            
        elif opt in ["-rgv","--rsc_cube_gen_vals"]:
            rsc_cube_gen_vals = str(arg)                        
        elif opt in ["-imd","--int_mod"]:
            int_mod = str(arg)
    return n_num,rsc_num,shape,sp_num,per_sp_int,ratio_cf,int_minmax,evt_num,evt_dur_var,run_time,sim_steps,sp_ini_perc,sp_ini_abd,D_set,neg_grw_thr,radius,nrep,kcoeff,neg_frac_kcoeff,squared_rsc,pref_type,kmod,rsc_cube_gen_rule,rsc_cube_gen_vals,int_mod
#pwl_spc,


###############################################################################################
### Additional code for calculating or generating required parameters/data                  ###
###############################################################################################

### Initializing seed populations
## Used in 'set_additional_conditions' function below
def set_seed_pop(n_num, sp_ini_num, sp_ini_abd, sp_num):
    n_pop = [] #np.empty(shape=(n_num))
    for i in range(n_num): 
        n_pop_fill = np.random.choice(sp_num, sp_ini_num, replace=False)
        n_pop_empty = np.repeat(0.0, sp_num)
        n_pop_empty[n_pop_fill] = sp_ini_abd
        n_pop.append(n_pop_empty)  # ini_pop
    seed_pop = np.array(n_pop)
    return seed_pop

### Get additional settings calculated from initial inputs
def set_additional_conditions(sp_num,per_sp_int,ratio_cf,evt_num,evt_dur_var,sp_ini_perc,sp_ini_abd,n_num,kcoeff,neg_frac_kcoeff,rsc_cube_gen_vals,rsc_cube_gen_rule,run_time,sim_steps,D_set):
    tot_perc_int = sp_num*per_sp_int/(((sp_num*sp_num)-sp_num)/2)       # Gives the absolute total number of species interactions      
    if tot_perc_int > 1:
        print('More interactions to input than can possibly filled')
        tot_perc_int = 1
    perc_comp = ratio_cf*tot_perc_int                                   # Percent competition
    perc_facil = (1-ratio_cf)*tot_perc_int                              # Percent mutualsim
    sp_ini_num = int(np.floor(sp_num * sp_ini_perc))                    # Total number of species initialised in habitats
    seed_pop = set_seed_pop(n_num, sp_ini_num, sp_ini_abd, sp_num)      # Set seep species abundances in habitats -> First abundance table for gLV
    neg_grw_thr = kcoeff*neg_frac_kcoeff                                # Threshold for K_{i,h} -> Below which species growth rates are set zero
    if rsc_cube_gen_rule != 'fix_same':
        rsc_cube_gen_vals = rsc_cube_gen_vals.split(',')                # If rsc_cube_gen_rule != 'fix_same' this creates a vector from the string input (separated by commas)
    t_dg = run_time/sim_steps                                           # Integration timesteps for gLV are set by dividing the total run time by the number steps for simulation
    t_g = run_time/sim_steps                                            # This can be ignored as only 't_dg' is used for simulation in the paper
    evt_dur = int(sim_steps)                                            # Creates integer for the number of simulation steps
    D = D_set*t_dg                                                      # The strength of diffusion per time step -> Needs to be modified by intergration timestep size!
    return tot_perc_int, perc_comp, perc_facil, evt_dur, sp_ini_num, seed_pop, neg_grw_thr, rsc_cube_gen_vals, t_dg, t_g, D

            


