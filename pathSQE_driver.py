import numpy as np
import seekpath
import os
import sys
import imp
from mantid.simpleapi import *
import time
import shutil

import pathSQE_input
imp.reload(pathSQE_input)
import pathSQE_functions
imp.reload(pathSQE_functions)
import define_data
imp.reload(define_data)
from reduce_data_to_MDE_07142023 import *
from slice_utils_07142023 import *


########################################################################################################
# Driver script for running pathSQE main program in conjuction with a pathSQE input file
# Author: Aiden Sable. Sept 2024.
########################################################################################################

# To run pathSQE, you need the following files in the same directory as this script:
#       pathSQE_input.py
#       pathSQE_functions.py
#       define_data.py
#       POSCAR ( if pathSQE_params['use_seeKpath_path'] = True )
#       reduce_data_to_MDE_07142023.py
#       slice_utils_07142023.py
#
# Reach out to Aiden for the versions of these files that are known to work properly with pathSQE


############### Load parameters defined in the pathSQE input file ###############
# Find path, find spacegroup, load mde data, and create saveDir's if don't already exist
start_time = time.time()

pathSQE_params = pathSQE_input.define_pathSQE_params()
#print(pathSQE_params)

if not os.path.exists(pathSQE_params['saveDir']):
    os.makedirs(pathSQE_params['saveDir'])
if not os.path.exists(pathSQE_params['saveDir']+'BZ_Reports/'):
    os.makedirs(pathSQE_params['saveDir']+'BZ_Reports/')
if not os.path.exists(pathSQE_params['saveDir']+'nxs_files/'):
    os.makedirs(pathSQE_params['saveDir']+'nxs_files/')
if not os.path.exists(pathSQE_params['saveDir']+'nxs_files/DataAndNorms/'):
    os.makedirs(pathSQE_params['saveDir']+'nxs_files/DataAndNorms/')


# To write output to terminal and an output.txt at the same time, use this class
class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # Flush to ensure the output is written immediately

    def flush(self):
        for f in self.files:
            f.flush()

# Define the output file path
output_file_path = os.path.join(pathSQE_params['saveDir'], 'output.txt')

# Open the file in write mode
output_file = open(output_file_path, 'a')

# Redirect stdout to both the terminal and the file
sys.stdout = Tee(sys.stdout, output_file)


if pathSQE_params['use_seeKpath_path']:
    poscar = pathSQE_functions.simple_read_poscar('POSCAR')
    path = seekpath.get_path(structure=poscar)
else:
    path = pathSQE_params['user_defined_path']

mtd_spacegroup = pathSQE_functions.find_matching_spacegroup(pathSQE_params)

mde_data = define_data.define_data_set()
reduce_data_to_MDE(mde_data)

############### INTRA-BZ folding for each BZ specified ###############
if pathSQE_params['fold_one_zone']:
    BZ_list = [pathSQE_params['BZ_offset']]
elif pathSQE_params['fold_all_zones']:    
    BZ_list_init = pathSQE_functions.gen_BZ_coverage_list_threshold(mde_data, pathSQE_params['H_bound'], pathSQE_params['K_bound'], pathSQE_params['L_bound'], pathSQE_params['E_bound'])
    
    # Problematic BZ arrays to remove and exclude from folding
    #rows_to_remove = [np.array([0, 0, 2]), np.array([1, 1, 1])]
    # mask to filter out specified rows
    #mask = np.all(BZ_list[:, None] == rows_to_remove, axis=-1)
    # Remove specified rows
    #BZ_list = BZ_list[~mask.any(axis=1)]

print('\nFolding {} BZs'.format(len(BZ_list_init)))
print('\nBZs: ', BZ_list_init)

print("\nChecking for prior folding progress")
BZ_list = pathSQE_functions.resume_analysis(BZ_list_init, pathSQE_params['saveDir']+'BZ_Reports/')


last_BZ_pathSeg_bounds = []
last_BZ_fold_dsl = []
for B, BZ_offset in enumerate(BZ_list):
    print("\nStarting BZ: ", BZ_offset)
    dsl_fold = []
    all_slice_names_in_BZ = []
    all_slice_evals = []
    for i in range(len(path['path'])):
        path_seg = path['path'][i]
        print('\nFolding BZ {}, path segments from {} to {}'.format(BZ_offset, path_seg[0], path_seg[1]))

        pt1_init = np.array(path['point_coords'][path_seg[0]])
        pt2_init = np.array(path['point_coords'][path_seg[1]])

        pt1_array, pt2_array = pathSQE_functions.generate_unique_paths(mtd_spacegroup, path_seg, pt1_init, pt2_init)

        all_slice_names_for_pathSeg = []
        slice_evals_for_pathSeg = []
        good_slice_index = None
        for j in range(pt1_array.shape[0]):
            q_dims_and_bins = pathSQE_functions.choose_dims_and_bins(pathSQE_params=pathSQE_params, point1=pt1_array[j], point2=pt2_array[j], perp_to_path=True, BZ_offset=BZ_offset)
            slice_desc = pathSQE_functions.make_slice_desc(pathSQE_params, q_dims_and_bins, pt1_array[j], pt2_array[j], path_seg)
            all_slice_names_for_pathSeg.append(slice_desc['Name'])
            
            make_slice(mde_data[0], slice_desc, ASCII_slice_folder='', MD_slice_folder='')
            
            # evaluate quality of slice based on filters
            slice_good = pathSQE_functions.evaluate_slice_quality(slice_name=slice_desc['Name'], filters=pathSQE_params['filters'], path_seg=path_seg)
            
            # only use good slices in the folding
            if slice_good:
                slice_desc['good_slice'] = True
                slice_evals_for_pathSeg.append(True)
                if good_slice_index is None:
                    # First good slice, clone the data
                    data = mtd['_data'].clone()
                    norm = mtd['_norm'].clone()
                    dsl_fold.append(slice_desc.copy())
                    dsl_fold[i]['Name'] = '{}2{}_BZ{}_final'.format(path_seg[0], path_seg[1], BZ_offset)
                    good_slice_index = j  # Save the index of the first good slice
                    
                    # save endpoints of each segment if this is the last BZ
                    if B == len(BZ_list)-1:
                        last_BZ_pathSeg_bounds.append(slice_desc['Dimension0Binning'])
                        last_BZ_fold_dsl.append(slice_desc.copy())
                    
                else:
                    # Subsequent good slices, combine data
                    data, norm = pathSQE_functions.combine_data_within_bz(data, norm, mtd['_data'], mtd['_norm'])
            else:
                slice_desc['good_slice'] = False
                slice_evals_for_pathSeg.append(False)

        # Save the data for the last good slice (if any)
        if good_slice_index is not None:
            saveDir = pathSQE_params['saveDir']+'nxs_files/'
            SaveMD(data, Filename=pathSQE_params['saveDir'] + 'nxs_files/DataAndNorms/{}2{}_BZ{}_data.nxs'.format(path_seg[0], path_seg[1], BZ_offset), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
            SaveMD(norm, Filename=pathSQE_params['saveDir'] + 'nxs_files/DataAndNorms/{}2{}_BZ{}_norm.nxs'.format(path_seg[0], path_seg[1], BZ_offset), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)

            output_ws = data / norm
            name = '{}2{}_BZ{}_final.nxs'.format(path_seg[0], path_seg[1], BZ_offset)
            SaveMD(output_ws, Filename=pathSQE_params['saveDir'] + 'nxs_files/' + name, SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
            LoadMD(Filename=pathSQE_params['saveDir'] + 'nxs_files/' + name, OutputWorkspace=name, LoadHistory=False)
            
        all_slice_names_in_BZ.append(all_slice_names_for_pathSeg)
        all_slice_evals.append(slice_evals_for_pathSeg)

    if pathSQE_params['BZ_report']:
        print('\nGenerating BZ Report for BZ {}'.format(BZ_offset))
        pathSQE_functions.generate_BZ_Report_morePages(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold)




############### INTER-BZ folding over all BZ specified ###############
if pathSQE_params['fold_all_zones']: 
    print('\nCombining all {} BZ'.format(len(BZ_list_init)))
    data_and_norm_dir = pathSQE_params['saveDir']+'nxs_files/DataAndNorms/'
    all_files = os.listdir(data_and_norm_dir)

    last_BZ = BZ_list_init[-1]
    last_BZ_files = [filename for filename in all_files if str(last_BZ) in filename]
    BZ_list_files = [filename for filename in all_files if filename not in last_BZ_files]

    # generate seg names
    seg_names = []
    for i in range(len(path['path'])):
        seg_names.append('{}2{}'.format(path['path'][i][0], path['path'][i][1]))

    #print(last_BZ_pathSeg_bounds)
    # for each seg, combine all BZ
    for path_seg in seg_names:
        # Load one good zone to start (data and norm) 
        pathSeg_files_lastBZ = [file for file in last_BZ_files if path_seg in file]
        first_data_file =  [file for file in pathSeg_files_lastBZ if 'data' in file]
        first_norm_file = [file for file in pathSeg_files_lastBZ if 'norm' in file]
        data=LoadMD(Filename=data_and_norm_dir+first_data_file[0])
        norm=LoadMD(Filename=data_and_norm_dir+first_norm_file[0])
        #print(path_seg, last_BZ, mtd['data'].getNonIntegratedDimensions)

        # Get list of all file names in the folder for data and norm (except for last BZ)
        all_files = [file for file in BZ_list_files if path_seg in file]
        data_files = [file for file in all_files if 'data' in file]
        norm_files = [file for file in all_files if 'norm' in file]
        #print(len(data_files))

        # For loop over all other zones/directions: 
        for i in range(len(data_files)):
            data_temp=LoadMD(Filename=data_and_norm_dir+data_files[i], LoadHistory=False) 
            norm_temp=LoadMD(Filename=data_and_norm_dir+norm_files[i], LoadHistory=False)  
            
            data.setSignalArray(data.getSignalArray() + data_temp.getSignalArray())
            norm.setSignalArray(norm.getSignalArray() + norm_temp.getSignalArray())
            data.setErrorSquaredArray(data.getErrorSquaredArray() + data_temp.getErrorSquaredArray())
            
        # Normalize the summed data and save
        folded_ws = data / norm
        SaveMD(folded_ws, Filename=pathSQE_params['saveDir']+'{}_folded.nxs'.format(path_seg))
        LoadMD(Filename=pathSQE_params['saveDir']+'{}_folded.nxs'.format(path_seg), OutputWorkspace=path_seg, LoadHistory=False)

    # Adaptive colorbars
    step = float(pathSQE_params['E_bins'].split(',')[1])
    ind_aboveElasticLine = int(np.ceil(3/step))
    pathMax = -1
    pathMin = 1e4
    for seg in seg_names:              
        slice = mtd[seg].getSignalArray()[:,0,0,ind_aboveElasticLine:]
        # Mask NaN and zero values
        mask = ~np.isnan(slice) & (slice != 0)
        segMin = np.percentile(slice[mask],10)
        segMax = np.percentile(slice[mask],95)
        if segMax > pathMax:
            pathMax = segMax
        if segMin < pathMin:
            pathMin = segMin

    fig = pathSQE_functions.plot_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, Ei=pathSQE_params['Ei'], vmi=pathMin, vma=pathMax, cma=pathSQE_params['cmap'])
    fig.savefig(pathSQE_params['saveDir']+'folded_path_plot.png')
    fig2 = pathSQE_functions.plot_SED_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, pathSQE_params=pathSQE_params)
    fig2.savefig(pathSQE_params['saveDir']+'folded_SED_path_plot.png')

print('\nFolded {} BZ in {} seconds'.format(len(BZ_list_init), time.time()-start_time))

# Remember to close the file when done
output_file.close()



# Define the current folder and the output folder from pathSQE_params
current_folder = os.getcwd()  # Get the current working directory
outputs_folder = pathSQE_params['saveDir']  # Use the pre-defined output path

# Loop over all files in the current folder
for filename in os.listdir(current_folder):
    file_path = os.path.join(current_folder, filename)

    # Check if it's a file (not a directory)
    if os.path.isfile(file_path):
        # Copy the file to the outputs folder
        shutil.copy(file_path, outputs_folder)

print(f"All files from {current_folder} copied to {outputs_folder}")