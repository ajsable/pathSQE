import numpy as np

########################################################################################################
# Define dictionary input file for running pathSQE
# Author: Aiden Sable. Sept 2024.
########################################################################################################

# To run pathSQE_driver.py, you need the following files in addition to this script:
#       pathSQE_driver.py
#       pathSQE_functions.py
#       define_data.py
#       reduce_data_to_MDE.py
#       slice_utils.py
#       POSCAR (if pathSQE_params['use_seeKpath_path']=True)
#
# Reach out to Aiden for copies of these files that are known to work properly with pathSQE


def define_pathSQE_params(**kwargs):
    pathSQE_params={
    # a few required inputs
    'sample':'SnGeTe',
    'spacegroup':'Fm-3m', # don't worry about spaces - code can theoretically figure it out 
    'Temp':300,
    'Ei':40, 
    'use_seeKpath_path':False, # current only works with False
    'which_seeKpath_pathSegments':None, # seeKpath segments to use
    'user_defined_path':{'path':[ ('X', 'Gamma'), ('Gamma', 'W'), ('W', 'L'), ('L', 'Gamma'), ('Gamma', 'K') ], 
                         # defined in terms of the conv cell (not primitive that pathSQE defaults to)
                         # original prim point coords: 'GAMMA': [0.0, 0.0, 0.0], 'X': [0.5, 0.0, 0.5], 'L': [0.5, 0.5, 0.5], 'W': [0.5, 0.25, 0.75], 'K': [0.375, 0.375, 0.75], 'U': [0.625, 0.25, 0.625]
                         'point_coords':{ 'X': [0.0, 1.0, 0.0], 'Gamma': [0.0, 0.0, 0.0], 'W': [0.5, 1.0, 0.0],'L': [0.5, 0.5, 0.5],'Gamma': [0.0, 0.0, 0.0],'K': [0.75, 0.75, 0.0]}},
    
    # bins and integration ranges
    'qdim0_step_size':0.025, # fraction of qdim0
    'E_bins':'0,0.25,17', # min, step, max
    'qdim1_int_range':0.1, # fraction of qdim1
    'qdim2_int_range':0.1, # fraction of qdim2
    
    # saving outputs
    'BZ_report':True,
    'saveDir':'/SNS/CNCS/IPTS-25269/shared/Aiden/comprehensive/BZfold_pathSQE/300K_thickInt_exclude2BZ_fullPath/',
    
    # parameters relevant to plotting
    'vmin':2e-3, # 
    'vmax':2e-2, # 
    'cmap':'viridis',
    
    # for automatic slice quality evaluation 
    'filters':[], # under development
    
    # folding within 1 BZ based on point group (intra-BZ folding)
    'fold_one_zone':False,
    'BZ_offset':np.array([3,3,1]),
    
    # folding multiple BZ's together (inter-BZ folding - each BZ is folded separately then all are combined)
    'fold_all_zones':True, # temp hard code to only do BZ_list = 3,3,1 and 5,1,1
    'H_bound':[-2,6], # range of H to check for data
    'K_bound':[-2,6], #   " "    K         " "
    'L_bound':[-2,6], #   " "    L         " " 
    'E_bound':[0,17]}  #   " "    E         " "

    return pathSQE_params




