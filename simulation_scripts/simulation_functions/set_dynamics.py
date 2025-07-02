###############################################################################################
#### SAMSARA - Seeting up population dynamics & diffusion models                            ###
###############################################################################################
import numpy as np
from scipy.integrate import ode
#from scipy.integrate import solve_ivp



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Full simulation function                                                                ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### We start with the full simulation function -> Below you can find the individual components called in this function
def full_sim_no_error_stop(t_g, t_dg, evt_num, evt_dur, D, n_num, sp_int_mat, seed_pop, n_nbr_idx, sp_hab_carcap, sp_grwrt):
    
    ### For the first run we want to document multiple timepoints of the simulations to see the trajectories of species abundances and to visually estimate whether abundances have reached a point of equillibration
    # Depending on the number of simulation steps we set (sim_steps = evt_dur) we divide select 101 timepoints in equal distances from start to finish
    # On each of these steps species abundances is saved
    savesteps = np.round(np.linspace(0, evt_dur, num = 101))
    n_pop_log   = []                        # Empty result array
    error_count = 0                         # Error at first timepoint (t0) 
    n_pop = seed_pop                        # Initial population density
    n_pop_log.append(n_pop.tolist())        # Becomes first entry in result array
    
    ## Here we run the very first simulation step, before we start the dispersal x gLV simulation loop
    t = t_dg                                # This is the step size for the gLV simulations
    # Main gLV simulation function
    n_pop, error_count = g_sim(t, n_num, n_pop, sp_int_mat, sp_hab_carcap, sp_grwrt, error_count)     
    n_pop_log.append(n_pop)                 # Saves first step
    
    countrr = 0
    ## Here we begin with the main simulation loop -> We alternate diffusion and gLV simulation with very narrow timesteps
    # effectively simulating both processes in parallel
    
    # For every simulation step (evt_dur minus the first)
    #e = 0
    for e in range(evt_dur-1):
        
        ## Diffusion simulation     -> See below
        n_pop = d_sim(n_num, n_nbr_idx, n_pop, D)
        
        ## gLV simulation           -> See below
        n_pop, error_count = g_sim(t, n_num, n_pop, sp_int_mat, sp_hab_carcap, sp_grwrt, error_count) 

        # If the current simulation step is in the vector of steps I want to save (e in savesteps) -> then save this result
        if (e in savesteps):                    # and (countrr % 2) == 0    ; countrr += 1
            n_pop_log.append(n_pop)
    
    ## We save both final diffusion and gLV steps to see how each process effects abundances in the narrow timesteps chosen and when species are equlibrated
    n_pop = d_sim(n_num, n_nbr_idx, n_pop, D)                                
    n_pop_log.append(n_pop)
    n_pop, error_count = g_sim(t, n_num, n_pop, sp_int_mat, sp_hab_carcap, sp_grwrt, error_count) 
    n_pop_log.append(n_pop)

    return np.array(n_pop_log), error_count
# The error count gives us an idea whether the gLV simulation has produced unrealistic population 'explosions', 
# a phenomenon that is sometimes encountered with gLV simulations -> see more below



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Diffusion simulation                                                                    ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### Diffusion or what we refer to in the paper as dispersal is very simply implemented (equation can be found in the Supporting material)
# A species abundance in a given habitat is compared to the same species abundances in a neighbouring (connected) habitat
# -> The difference is assessed and multiplied with the diffusion coefficient
# The sum of all differences to all neighbouring (connected) habitats is then added to a local habitats species abundance
# This means if this species abundance is highest in the assessed habitat it will lose abundance to its neighbouring habitats (source population)
# If a habitat is lower in this species abundance than its neighbours then it is a sink habitat
def d_sim(n_num, n_nbr_idx, n_pop, D):
    n_pop = np.array(n_pop)
    n_pop_new = []
    for i in range(n_num):
        if len(n_nbr_idx[i]) > 0:

            add_sub = []                
            for j in range(len(n_nbr_idx[i])):                    
                add_sub.append(D * ( n_pop[n_nbr_idx[i][j]] - n_pop[i] ))    
            n_pop_new.append(list(n_pop[i] + (np.array(add_sub).sum(axis=0))))            
        else:
            n_pop_new.append(list(n_pop[i]))
    n_pop_new = np.array(n_pop_new)
    n_pop_new[n_pop_new<0.001] = 0
    n_pop_new = n_pop_new.tolist()      
    return(n_pop_new)



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### gLV function                                                                            ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### We solve the gLV differential equations using the ode.set_integrator function from scipy.integrate
# The method we use is 'dopri5' which applies the runge-kutta method

### The main gLV integration function
def glv(t, x, r, sp_int_mat, k):    
    dxdt = ((r*x)/k) * (k - x - (sp_int_mat @ x) )
    return dxdt
# This function is analogous to Eq.1 paper, but has been slightly rearranged (k outside of fraction)

# Keep in mind that in the paper interaction coefficients signs are switched! 
# Here in the code negative interaction coefficients are mutualistic, while postiive coefficients are competitive!

### Integration solver main function
def solve_glv(r, sp_int_mat, k, x0, t_end):
    def solout(t, y):
        sol.append([*y])
    solver = ode(glv).set_integrator('dopri5')
    #def int_glv(x0, sp_int_mat, r, k, t_end): 
    solver.set_initial_value(x0, 0).set_f_params(r, sp_int_mat, k)
    sol = []
    solver.set_solout(solout)
    solver.integrate(t_end)
    sol = np.array(sol)
    return sol

### Main function for simulation population dynamics with gLV
def g_sim(t, n_num, n_pop, sp_int_mat, sp_hab_carcap, sp_grwrt, error_count):   
    sol = []
    for i in range(n_num): 
        k = sp_hab_carcap[i]
        r = sp_grwrt[i]
        x0 = n_pop[i]
        
        solx = solve_glv(r, sp_int_mat, k, x0, t)
        
        ## We have added this 'warning' here when species abundances far exceed what are reasonably expected values
        # This is an issue inherent to gLV simulations with randomised interaction matrices
        # Depending on the number of positive interactions and interaction strengths positive feedbacks can produce 
        # 'explosive' growth and unrealisitcally high abundances
        if (solx.max() > 1000) or (solx.min() < 0):
            print('ERROR')
            error_count += 1      
        # Runs with errors will be excluded from subsequent analysis!    
        
        # Species with abundances below 0.001 will be set to zero after every simulation time step
        solx[-1][solx[-1]<0.001] = 0
        sol.append(list(solx[-1]))
    return sol, error_count


