#######################################################################################
### SAMSARA - Set species characteristics                                           ###
#######################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Load libraries & setwd for loading functions                                            ###
import numpy as np
import random
from scipy.stats import truncnorm
from scipy.spatial.distance import braycurtis



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Set random species interaction matrices                                                 ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### Set random species interaction matrices 
# Generates a fully symmetric matrix with 0 for non interacting speces and draws random values for interacting species
# These values are between 0.1 and 0.5 (input fixed in paper -> max_int_vals = 0.5) and -0.1 and -0.5 for competing and mutualistic pairs respectively
def rdm_sym_sp_int_mat(sp_num, perc_comp, perc_facil, max_int_vals, int_mod):
    int_idx = np.tril_indices(sp_num,k=-1)
    repl_comp = int(np.round( perc_comp * (len(int_idx[0])) ))
    repl_facil = int(np.round( perc_facil * (len(int_idx[0])) ))
    repl_idx = np.random.choice(len(int_idx[0]), repl_comp+repl_facil, replace = False)
    repl_mat_idx = []
    for i in range(repl_comp+repl_facil):
        repl_mat_idx.append((int_idx[0][repl_idx[i]],int_idx[1][repl_idx[i]]))
    if max_int_vals == 0.:
        sp_int_mat = np.zeros(shape=(sp_num,sp_num))
        for i in range(repl_comp):
            sp_int_mat[repl_mat_idx[i]] = np.random.uniform(0., max_int_vals)
        for i in range(repl_comp, repl_comp+repl_facil):
            sp_int_mat[repl_mat_idx[i]] = np.random.uniform(0., -max_int_vals)
    else:                                   # Other tested options deleted from here as this is the only method used in the paper
        if int_mod == 'rdm_unif':
            sp_int_mat = np.zeros(shape=(sp_num,sp_num))
            for i in range(repl_comp):
                sp_int_mat[repl_mat_idx[i]] = np.random.uniform(0.1, max_int_vals)
            for i in range(repl_comp, repl_comp+repl_facil):
                sp_int_mat[repl_mat_idx[i]] = np.random.uniform(-0.1, -max_int_vals)
    sp_int_mat = sp_int_mat + sp_int_mat.T - np.diag(np.diag(sp_int_mat))
    return sp_int_mat
# Keep in mind that in the paper interaction coefficients signs are switched! 
# Here in the code negative interaction coefficients are mutualistic, while postiive coefficients are competitive!



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Species resource preferences                                                            ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### Main function for generating species resource preferences
def set_prefernce(sp_num, rsc_num, pref_type):
    if pref_type == 'identical':
        sp_rsc_prf = set_identical_preference(sp_num, rsc_num)
    elif pref_type == 'independent':
        sp_rsc_prf = set_independent_preference(sp_num, rsc_num)
    elif pref_type == 'sequential':
        sp_rsc_prf = set_sequential_preference(sp_num, rsc_num)
    else:
        print("ERROR in set_prefernce input")
        return
    return sp_rsc_prf
# This function can generate species resource preferences in three different based on the 'pref_type' input -> See next

## Set species resource preference completely identical -> Preference is an identical fraction of 1 
# If there are three resources (rsc_num = 3) each speces preference for each resource is 0.333333 (not used in paper)
def set_identical_preference(sp_num, rsc_num):
    sp_rsc_prf = []
    for i in range(sp_num):
        sp_i_prf = []
        for j in range(rsc_num):
            sp_i_prf.append(1/rsc_num)
        sp_rsc_prf.append(list(sp_i_prf))
    return sp_rsc_prf

## Set species resource preference independently -> Preferences are randomly drawn from 0 to 1, then normalised so the total sum is 1
# This produces more similar, generalist resource preferences between species (not used in paper)
def set_independent_preference(sp_num, rsc_num):
    sp_rsc_prf = []
    for i in range(sp_num):
        sp_i_prf = []
        for j in range(rsc_num):
            sp_i_prf.append(random.uniform(0,1))
        sp_rsc_prf.append(list(np.array(sp_i_prf)/sum(sp_i_prf)))
    return sp_rsc_prf

## Set species resource preference sequentially -> Method described in paper
def set_sequential_preference(sp_num, rsc_num):
    sp_rsc_prf = []
    for i in range(sp_num):
        sp_i_prf = []
        sp_i_prf.append(random.uniform(0, 1))
        for j in range(rsc_num-2):
            sp_i_prf.append(random.uniform(0, 1 - sum(sp_i_prf)))   
        sp_i_prf.append(1 - sum(sp_i_prf))    
        random.shuffle(sp_i_prf)
        sp_rsc_prf.append(sp_i_prf)
    return sp_rsc_prf

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Get Euclidean distance matrix -> Gives environmental preference similarity values E_{ij}
def eudi_rscprf(sp_rsc_prf):
    eudimat = np.zeros(shape = (len(sp_rsc_prf),len(sp_rsc_prf)))
    for i in range(len(sp_rsc_prf)):
        for j in range(len(sp_rsc_prf)):
            eudimat[i,j] = np.linalg.norm(np.array(sp_rsc_prf[i]) - np.array(sp_rsc_prf[j]))        
    return eudimat



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Species habitat carrying capacity K_{i,h}                                               ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### Main function to set species habitat carrying capacity
## Calculation can be found in Eq.2 in the paper
def set_sp_hab_carcapx(n_num,sp_rsc_prf,n_knorm,n_kabs,neg_grw_thr,sp_ini_abd, kcoeff):
    sp_hab_carcap = []
    for i in range(n_num): 
        calc_sp_hab_carcap = np.array((np.array(sp_rsc_prf) @ np.array(n_kabs[i])))
        calc_sp_hab_carcap[calc_sp_hab_carcap <= neg_grw_thr] = sp_ini_abd
        sp_hab_carcap.append(list(calc_sp_hab_carcap))
    return sp_hab_carcap

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Here we add the ability to override species habitat carrying capacities for the control scenario c2
# By randomly generating K_{i,h}'s for every species we allow some variability in species abundances
# and therefore impacts from interactions between samples

## First we need a function to get values from a truncated normal distribution
# As we dont want negative K_{i,h} we truncate the normally distributed values we generate. 
# The distributions parameters are based of the main reference scenario's distribution (s1) 
def get_truncnorm_val(mean, sd, low, upp):
    trunc_val = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
    return trunc_val

## Main function for augumenting species K_{i,h}'s
# In the paper we only use kmod = 'rdm_normal' for c1 -> Values for the distribution from main reference scenario's distribution (s1)
# And in all other cases we simply set kmod = 'keep' which keeps the generated K_{i,h}'s unaltered
def mod_sp_hab_carcap(sp_hab_carcap, kmod, kcoeff):
    if kmod == 'equal':
        sp_hab_carcap = np.full(np.shape(sp_hab_carcap), kcoeff).tolist()
    elif kmod == 'rdm_uniform':
        sp_hab_carcap = np.random.uniform(size = np.shape(sp_hab_carcap), low = 0.19, high = 0.35).tolist()
    elif kmod == 'rdm_normal':
        tval = get_truncnorm_val(mean=0.27, sd=0.13, low=0, upp=kcoeff) # Informed by distribution of diffL = L(mix)
        sp_hab_carcap = tval.rvs(np.shape(sp_hab_carcap)).tolist()
    elif kmod == 'keep':
        print('K is kept as resource landscape based')
        sp_hab_carcap = sp_hab_carcap
    else:
        print('Error: Wrong input for modifying species habitat carrying capacity')
    return sp_hab_carcap

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Set species growth rates
# When species are in a very unfavourable habitat -> So when K_{i,h} is below our threshold value (neg_grw_thr = 0.1)
# we want them not to grow in this habitat. Species with K_{i,h}'s above this threshold have the growth rate r = 1, those below get r = 0
def set_sp_grwrt_bin(n_num,sp_hab_carcap, neg_grw_thr, killem):
    sp_grwrt = []
    for i in range(n_num):
        sp_hab_carcap_x = np.copy(sp_hab_carcap[i])
        if killem == True:
            sp_hab_carcap_x[sp_hab_carcap_x < (neg_grw_thr)] = -1
        else:
            sp_hab_carcap_x[sp_hab_carcap_x < (neg_grw_thr)] = 0
        sp_hab_carcap_x[sp_hab_carcap_x >= (neg_grw_thr)] = 1
        sp_grwrt.append(sp_hab_carcap_x.tolist())
    return sp_grwrt


