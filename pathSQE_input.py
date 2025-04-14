import numpy as np

########################################################################################################
# Define dictionary input file for running pathSQE
# Author: Aiden Sable. April 2025.
########################################################################################################

# To run pathSQE_driver.py, you need the following files in addition to this script:
#       pathSQE_driver.py
#       pathSQE_functions.py
#       define_data.py
#       reduce_data_to_MDE.py
#       slice_utils.py
#       POSCAR (if 'use seeKpath path'=True or 'perform simulations=True)
#       FORCE_CONSTANTS (if 'perform simulations=True)


def define_pathSQE_params(**kwargs):
    pathSQE_params={
    # Sample and Q point info
    'sample':'Ge', # separate elements with spaces; don't include stoichiometric subscripts, e.g. for Sr3As2, do 'Sr As'; if desired, some specific isotopes can be specified, e.g. instead of 'Fe Si', could do '56Fe 28Si'
    'space group':'Fd-3m', # don't worry about spaces
    'T and Ei conditions':[(5,40)], # (Temperature, Indicent energy) for each dataset in a list. MDEs must follow consistent naming scheme you specify in define_data.py
    'use seeKpath path':False, # POSCAR required
    'user defined Qpoints':{'path':[ ('Gamma', 'X'), ('X', 'U'), ('K', 'Gamma'), ('Gamma', 'L'), ('L', 'W'), ('W', 'X') ], # e.g. [('Gamma', 'X'), ('X', 'U')], total piece-wise path along which to process 2D inelastic slices, names must match strings in 'point_coords'
                            '1d_points':[], # e.g. ['Gamma', 'X'], all 1d points at which to make 1d S(E) cuts, names must match strings in 'point_coords'
                            'point_coords':{ 'Gamma': [0.0, 0.0, 0.0], 'X': [0, 1, 0], 'L': [0.5, 0.5, 0.5], 'W': [0.5, 1, 0], 'K': [0.75, 0.75, 0], 'U': [0.25, 1, 0.25] }}, # Q point definitions in Q basis associated with MDE's UB matrix
    'primitive to mantid transformation matrix': np.array([[-1,1,1],[1,-1,1],[1,1,-1]]), # matrix defining transformation from primtive reciprocal lattice basis to Q basis associated with MDE's UB matrix
    
    # bins and integration ranges
    'E bins':'0,0.5,40', # min, step, max for energy binning
    'qdim0 step size':0.025, # if 2d 'path', this is step along qdim0 in rlu; if 1d 'point_coords', this is +/- integration range for qdim0 in rlu
    'qdim1 integration range':0.1, # +/- integration range for qdim1 in rlu
    'qdim2 integration range':0.1, # +/- integration range for qdim2 in rlu

    # more advanced / comprehensive processing options
    'all symmetric 2d slices':True, # whether to process all symmetrically equivalent 2d inelastic slices (also includes bz folding)
    'all symmetric 1d cuts':False, # whether to process all symmetrically equivalent 1d cuts
    'BZ to process':0.89, # float or list of lists; if float between [0,1], then BZ with fractional qE coverage above provided threshold are processed (e.g. 0.5); if list of lists, data in each listed BZ will be processed (e.g. [[1,1,1], [2,0,0], [4,2,0]])
    'slice filter functions':[], # list of string(s), for automatic slice evaluation and filtering, string(s) must match filter function(s) specified in filter_functions.py
    
    # simulation details - POSCAR and FORCE_CONSTANTS required
    'perform simulations':True, # whether to perform analogous S(Q,E) simulations
    'supercell dimensions':[4,4,4], # supercell dimensions used to generate FORCE_CONSTANTS
    'use experimental coverage mask':True, # whether to impose experimental Q,E coverage mask on simulations
    'energy blurring sigma':0.75, # resolution blurring stdev in E, unit of E step (i.e. setting 2 here if E_bins[1] is 0.5 meV would give blurring with 1 meV stdev)
    
    # saving and outputs
    'output directory':'/SNS/ARCS/IPTS-13861/shared/Aiden/comprehensive/pathSQE_testing_Ge/5K_40meV_testRefactored/', # where to save everything (will be created if doesn't exist)
    'cmap':'viridis', # colormap for 2d plots

    # temporary - need to figure out how these can be accessed from loaded MDE workspace - ask Andrei
    'u_vec':np.array([1,0,0]),
    'v_vec':np.array([0,1,0])
    }

    return pathSQE_params




