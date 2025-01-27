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
    'sample':'Bi',
    'spacegroup':'R -3 m :r', # using rhombohedral 
    'Temp':2,
    'Ei':20, 
    'use_seeKpath_path':False, # current only works with False
    'which_seeKpath_pathSegments':None, # seeKpath segments to use
    'user_defined_path':{'path':[ ('Gamma','T'), ('T','X'), ('X','L'), ('L','Gamma'), ('Gamma','X') ], 
                         # defined in terms of the primitive rhombohedral basis (i.e. not the basis being used in mantid which is often the conventional basis)
                         'point_coords':{'Gamma': [0.0, 0.0, 0.0], 'T': [0.5, 0.5, 0.5], 'X': [0 , 0.5, 0.5], 'L': [0, 0.5, 0] }},
    'prim2mantid_transMatrix': np.linalg.inv(np.array([[2/3, -1/3, -1/3], [1/3, 1/3, -2/3], [1/3, 1/3, 1/3]])), # bilbao space group R-3m (166) to hexagonal conv used in mantid UB matrix
    
    
    # bins and integration ranges
    'qdim0_step_size':0.025, # fraction of qdim0
    'E_bins':'-15,0.2,15', # min, step, max
    'qdim1_int_range':0.05, # fraction of qdim1
    'qdim2_int_range':0.05, # fraction of qdim2
    
    # saving outputs
    'BZ_report':True,
    'saveDir':'/SNS/CNCS/IPTS-31965/shared/Aiden/comprehensive/BZfold_pathSQE/2K_20meV/',
    
    # parameters relevant to plotting
    #'vmin':1e-5, # 
    #'vmax':7e-4, # 
    'cmap':'viridis',
    
    # for automatic slice quality evaluation 
    'filters':[], # under development
    
    # folding within 1 BZ based on point group (intra-BZ folding)
    'fold_one_zone':False,
    'BZ_offset':np.array([3,3,1]),
    
    # folding multiple BZ's together (inter-BZ folding - each BZ is folded separately then all are combined)
    'fold_all_zones':True, 
    'H_bound':[-10,10], # range of H to check for data
    'K_bound':[-10,10], #   " "    K         " "
    'L_bound':[-10,10], #   " "    L         " " 
    'E_bound':[0,15]}  #   " "    E         " "

    return pathSQE_params




