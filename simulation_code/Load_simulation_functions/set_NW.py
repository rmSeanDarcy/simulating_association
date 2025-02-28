###############################################################################################
#### SAMSARA - Seeting up habitat network and node property calculations                    ###
###############################################################################################
import numpy as np
import random
import networkx as nx
from Fyld_generator import *
#from FyeldGenerator import *


###############################################################################################
#### Set up names
# Set habitat names
def get_simple_n_nms(n_num):
    n_nms = []
    for i in range(1,n_num+1):
        n_nms.append('H'+str(i))   
    return n_nms
# Set resource names
def get_simple_rsc_nms(rsc_num):
    rsc_nms = []
    for i in range(1,rsc_num+1):
        rsc_nms.append('R'+str(i))   
    return rsc_nms
# Set species names
def get_simple_sp_nms(sp_num):
    sp_nms = []
    for i in range(1,sp_num+1):
        sp_nms.append('S'+str(i))   
    return sp_nms
# Get all names
def get_nms(n_num,rsc_num,sp_num):
    sp_nms  = get_simple_sp_nms(sp_num)
    rsc_nms = get_simple_rsc_nms(rsc_num)
    n_nms   = get_simple_n_nms(n_num)
    return sp_nms, rsc_nms, n_nms  



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Network construction                                                                    ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##### 3D networks
### Random network based of of networkx's random_geometric_graph() function
# Random points are generated in within a space (L³ = 1)
def set_NW_random_geometric_graph(n_num, radius):
    
    # Generate a dict of random node positions
    pos = {i: (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(n_num)}
    # Create random 3D network -> radius dictates how close nodes are needed to be for edges to be established
    G = nx.random_geometric_graph(n_num, radius, pos=pos)
    n_nms = list(G.nodes) # Get vector (set) of node names (here simply indices)
    node_xyz = np.array([pos[v] for v in G])
    edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])
    n_nbr_nms = np.empty(shape=(n_num), dtype=object)
    n_nbr_idx = np.empty(shape=(n_num), dtype=object)  
    for i in range(n_num):
        n = list(n_nms)[i]
        n_nbr_nms[i] = list(G.adj[n])
        n_nbr_number = np.zeros(len(n_nbr_nms[i]), dtype=int)
        for j in range(len(n_nbr_number)):
            m = n_nbr_nms[i][j]
            n_nbr_number[j] = n_nms.index(m)
        n_nbr_idx[i] = n_nbr_number

    return G, pos, n_num, n_nms, node_xyz, edge_xyz, n_nbr_nms, n_nbr_idx
# Outputs are:
# G         =  the nx. graph object
# pos       = the positions of the nodes (habitats)
# n_num     = The number of habitats
# n_nms     = Names of habitats
# node_xyz  = Node positions
# edge_xyz  = Edge positions (not further needed)
# n_nbr_nms = List of neighboring habitat (i.e. connected to habitat in network) indices
# n_nbr_idx = List of neighboring habitat (i.e. connected to habitat in network) indices

### Cloned from aboove -> Allows to fixed the positions of habitats first (outside of this function) and then generate networks with different connectivity radii
def set_NW_random_geometric_graph_fixed_pos(n_num, radius, pos):
#    pos = {i: (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(n_num)}
    G = nx.random_geometric_graph(n_num, radius, pos=pos)
    n_nms = list(G.nodes) # Get vector (set) of node names (here simply indices)
    node_xyz = np.array([pos[v] for v in G])
    edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])
    n_nbr_nms = np.empty(shape=(n_num), dtype=object)
    n_nbr_idx = np.empty(shape=(n_num), dtype=object)  
    for i in range(n_num):
        n = list(n_nms)[i]
        n_nbr_nms[i] = list(G.adj[n])
        n_nbr_number = np.zeros(len(n_nbr_nms[i]), dtype=int)
        for j in range(len(n_nbr_number)):
            m = n_nbr_nms[i][j]
            n_nbr_number[j] = n_nms.index(m)
        n_nbr_idx[i] = n_nbr_number

    return G, pos, n_num, n_nms, node_xyz, edge_xyz, n_nbr_nms, n_nbr_idx



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Resource landscape generation                                                           ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

###############################################################################################
### Generating random gaussian 3D field -> Function taken from library: FyldGenerator
## from Fyld_generator import *
# See documentation on https://github.com/cphyc/FyeldGenerator

# Inputs are: 1. the power law spectrum (pwl_spc), which gives the prevalence of small scale spatial structures
# 2. the shape argument, which gives the dimension of voxels that are simulated. In our simulations this is always [100,100,100]
def set_rdm_gaussian_field_3d(pwl_spc, shape):
    n = pwl_spc
#    shape = dim_3D
    # Helper that generates power-law power spectrum
    def Pkgen(n):
        def Pk(k):
            return np.power(k, -n)
        return Pk
    # Draw samples from a normal distribution
    def distrib(shape):
        a = np.random.normal(loc=0, scale=1, size=shape)
        b = np.random.normal(loc=0, scale=1, size=shape)
        return a + 1j * b
    resource_cube = generate_field(distrib, Pkgen(n), shape)
    resource_cube = resource_cube[:] + abs(resource_cube.min())
    resource_cube = resource_cube / resource_cube.max()

    return resource_cube
# Basic output of FyldGenerator library

###############################################################################################
### Generate a simple cube of the same dimensions with a gradient -> Values from 0 to 1 from one side to another
def set_rdm_grad_cube(shape, val_set):
    mat = np.zeros(shape = shape)
    filldim = np.linspace(0, 1, shape[0])
    
    # The input val_set can take 6 values -> 'bf' or 'fb' generate these gradients from the back to front of the cubes and front to back respectively
    # Inputs are analogous for top-down and left-right
    if val_set == 'bf':
        mat[:] = filldim[:, np.newaxis]
    if val_set == 'fb':
        mat[:] = (1-filldim[:, np.newaxis])
    
    if val_set == 'td':
        mat[:] = filldim[:]
    if val_set == 'dt':
        mat[:] = (1 - filldim[:])

    for m in range(len(mat)):
        if val_set == 'rl':
            mat[m] = filldim[m]
        if val_set == 'lr':
            n = len(mat) - 1 - m
            mat[m] = filldim[n]
    resource_cube = mat
    
    return resource_cube
# Output data matches output format of FyldGenerator library

###############################################################################################
### Main function for processing the input to then generate resource cubes
def get_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_rule, rsc_cube_gen_vals):    
    if rsc_cube_gen_rule == 'fix_same':
        resource_cubes_x = fix_same_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_vals)    
    if rsc_cube_gen_rule == 'specific_dist':
        resource_cubes_x = specific_dist_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_vals)
    return resource_cubes_x  
## We have two possible inputs for 'rsc_cube_gen_rule' -> Related functions that are called here are explained next

# Set multiple resource cubes independently with the same spatial input (rsc_cube_gen_rule == 'fix_same') 
def fix_same_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_vals):
    resource_cubes_x = []
    for i in range(rsc_num):
        # rsc_cube_gen_vals -> Contains the power law spectrum input for generating the resource landscapes (pwl_spc)
        resource_cube = set_rdm_gaussian_field_3d(float(rsc_cube_gen_vals[0]), shape)
        # We added the input to keep the normal generated values (squared_rsc = 'normal') or to square values (squared_rsc = 'squared')
        if squared_rsc == 'normal':        
            resource_cubes_x.append(np.array(resource_cube))
        elif squared_rsc == 'squared':      
            resource_cubes_x.append(np.array(resource_cube**2))
    resource_cubes_x = np.array(resource_cubes_x)
    return resource_cubes_x

# If inputs were selected to differ between resource landscapes generated (rsc_cube_gen_rule == 'specific_dist') some logic is required
def specific_dist_resource_cubes(rsc_num, shape, squared_rsc, rsc_cube_gen_vals):
    resource_cubes_x = []
    
    # In the script './Load_simulation_functions/set_standard_initial_conditions.py' 
    # the function 'set_initial_conditions_basic' explains the logic of this 'rsc_cube_gen_vals' when rsc_cube_gen_rule == 'specific_dist'
    # Depending on the specific inputs cubes are generated for each resource
    for k in range(rsc_num):
        typ_set = rsc_cube_gen_vals[k].split('X')[0]
        sq_set = rsc_cube_gen_vals[k].split('X')[1]
        val_set = rsc_cube_gen_vals[k].split('X')[2]
        
        if typ_set == 'pwl_spc':        
            resource_cube = set_rdm_gaussian_field_3d(float(val_set), shape)        
            if sq_set == 'nml':
                resource_cubes_x.append(np.array(resource_cube))
            elif sq_set == 'sqr':
                resource_cubes_x.append(np.array(resource_cube**2))

        if typ_set == 'grad':        
            resource_cube = set_rdm_grad_cube(shape, val_set)        
            if sq_set == 'nml':
                resource_cubes_x.append(np.array(resource_cube))
            elif sq_set == 'sqr':
                resource_cubes_x.append(np.array(resource_cube**2))

    return resource_cubes_x



###############################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Assigning NW nodes voxel specific resource abundances                                   ###
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################

### Based on the position of habitats (in node_xyz) their corresponding voxel in each resource landscape is found
# Every habitat gets a corresponding resource abundance for every resource landscape
def set_node_resource_vals(n_num, node_xyz, shape, resource_cubes_x, rsc_num): #resource_cubes_i, 
    n_kabs = []
    for i in range(n_num):
        n_kj = []
        for j in range(rsc_num):
            vx = int(np.floor(node_xyz[i][0]*shape[0]))
            vy = int(np.floor(node_xyz[i][1]*shape[0]))
            vz = int(np.floor(node_xyz[i][2]*shape[0]))
            n_kj.append(resource_cubes_x[j][vx][vy][vz])
        n_kabs.append(n_kj)
    n_knorm = []
    for i in range(n_num):
        n_knorm.append(list(n_kabs[i]/sum(n_kabs[i])))

    return n_knorm, n_kabs
# Outputs are both the absolute resource abundacnes (k_abs) and the relative (k_norm)

