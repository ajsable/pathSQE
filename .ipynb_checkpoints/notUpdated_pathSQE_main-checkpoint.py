import os
import imp
import sys
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mantid.simpleapi import CreateMDHistoWorkspace
from mantid import plots
from pathSQE_functions import *

sys.path.append('/SNS/groups/dgs/DGS_SC_scripts')
from reduce_data_to_MDE import *
from slice_utils import *

sys.path.append('manual_mantid_reduceAndViz')
import define_data

#############################################
# Find path and make all slice descriptions
#############################################
dsl = make_all_slice_descs(path_to_poscar='files_for_nb/POSCAR')


#############################################
# Taking slices w/ Mantid
#############################################
imp.reload(define_data)

mde_data=define_data.define_data_set()
reduce_data_to_MDE(mde_data)

for slice_description in dsl:
    make_slice(mde_data[0],slice_description, ASCII_slice_folder='', MD_slice_folder='')

