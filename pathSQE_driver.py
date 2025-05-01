########################################################################################################
# Driver script for running pathSQE main program in conjuction with a pathSQE input file
# Author: Aiden Sable. April 2025.
########################################################################################################


import numpy as np
import seekpath
import os
import sys
sys.stdout = sys.__stdout__  # Reset stdout before the script exits
from mantid.simpleapi import *
import time
import shutil


import importlib as imp
import pathSQE_input
imp.reload(pathSQE_input)
import define_data
imp.reload(define_data)

# Import refactored code from pathSQE_utils/ folder
import pathSQE_utils.core as pathSQE_core
imp.reload(pathSQE_core)

import pathSQE_utils.plotting_and_reports as pathSQE_plotting_and_reports
imp.reload(pathSQE_plotting_and_reports)

import pathSQE_utils.helper as pathSQE_helper
imp.reload(pathSQE_helper)

import pathSQE_utils.simulations as pathSQE_simulations
imp.reload(pathSQE_simulations)

import pathSQE_utils.filter_functions as pathSQE_filter_functions
imp.reload(pathSQE_filter_functions)

import pathSQE_utils.reduce_data_to_MDE_07142023 as reduce_data_to_MDE_07142023
imp.reload(reduce_data_to_MDE_07142023)

import pathSQE_utils.slice_utils_07142023 as slice_utils_07142023
imp.reload(slice_utils_07142023)

#######################################################################################################
        ################## LOAD PARAMETERS, DATA, AND PREPARE FOR PROCESSING ###################
#######################################################################################################

# Find path, find spacegroup, load mde data, and create saveDir's if don't already exist, etc.
start_time = time.time()

pathSQE_params = pathSQE_input.define_pathSQE_params()
#print(pathSQE_params)

if not os.path.exists(pathSQE_params['output directory']):
    os.makedirs(pathSQE_params['output directory'])
if not os.path.exists(os.path.join(pathSQE_params['output directory'],'all_slices/')):
    os.makedirs(os.path.join(pathSQE_params['output directory'],'all_slices/'))
if not os.path.exists(os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/')):
    os.makedirs(os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/'))
if not os.path.exists(os.path.join(pathSQE_params['output directory'],'allSym_sim_files/')) and pathSQE_params['perform simulations']:
    os.makedirs(os.path.join(pathSQE_params['output directory'],'allSym_sim_files/'))

# Copy scripts being used to output folder
current_folder = os.getcwd()  # Get the current working directory
outputs_folder = pathSQE_params['output directory']  # Use the pre-defined output path
all_slices_dir = os.path.join(pathSQE_params['output directory'],'all_slices/')

print(f"All files from {current_folder} being copied to {outputs_folder}")
# Loop over all files in the current folder
for filename in os.listdir(current_folder):
    file_path = os.path.join(current_folder, filename)

    # Check if it's a file (not a directory)
    if os.path.isfile(file_path):
        # Copy the file to the outputs folder
        shutil.copy(file_path, outputs_folder)


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
output_file_path = os.path.join(pathSQE_params['output directory'], 'output.txt')

# Open the file in write mode and Redirect stdout to both the terminal and the file
output_file = open(output_file_path, 'a')
sys.stdout = Tee(sys.stdout, output_file)


if pathSQE_params['use seeKpath path']:
    poscar = pathSQE_helper.simple_read_poscar('POSCAR')
    user_defined_Qpoints = seekpath.get_path(structure=poscar)
else:
    user_defined_Qpoints = pathSQE_params['user defined Qpoints']

mtd_spacegroup = pathSQE_helper.find_matching_spacegroup(pathSQE_params)

###################################################
mde_data = define_data.define_data_set(pathSQE_params['T and Ei conditions'])
print('Processing ', len(mde_data), ' datasets.')
reduce_data_to_MDE_07142023.reduce_data_to_MDE(mde_data)


if pathSQE_params['all symmetric 2d slices'] or pathSQE_params['all symmetric 1d cuts']:
    if not os.path.exists(os.path.join(pathSQE_params['output directory'],'Reports/')):
        os.makedirs(os.path.join(pathSQE_params['output directory'],'Reports/'))

    if isinstance(pathSQE_params['BZ to process'], list):
        BZ_list_init = [np.array(BZ) for BZ in pathSQE_params['BZ to process']] 
    else:
        soft_time = time.time()
        BZ_list_init, BZ_fractional_coverages, BZ_full_array = pathSQE_core.find_BZ_with_data(mde_data[0], pathSQE_params)
        print('BZ cov and thresh in ', time.time()-soft_time)

    BZ_list_init = np.array(BZ_list_init)
    print("Initial BZ list: ", len(BZ_list_init), BZ_list_init.tolist())
    print("\nChecking for prior progress")
    BZ_list, removed_BZ = pathSQE_helper.resume_analysis(BZ_list_init, os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/'), len(user_defined_Qpoints['path']), all_slices_dir)
    if len(BZ_list) == 0:
        print("All BZ already processed - Terminating execution.")
        sys.exit(1)  # Exit with error status

    print('\nData being processed in {} BZs'.format(len(BZ_list)))
    print('\nBZs:', BZ_list.tolist())

sliceID = pathSQE_helper.get_next_sliceID(all_slices_dir)


#######################################################################################################
        ###################################### all symmetric 2d slices ###################################
#######################################################################################################

if pathSQE_params['all symmetric 2d slices']:
########## INTRA-BZ folding for each BZ specified ##########
    if not os.path.exists(os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/DataAndNorms/')):
        os.makedirs(os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/DataAndNorms/'))
    if not os.path.exists(os.path.join(pathSQE_params['output directory'],'folding_progess/')):
        os.makedirs(os.path.join(pathSQE_params['output directory'],'folding_progess/'))
    if pathSQE_params['perform simulations']:
        fully_folded_sims_data, fully_folded_sims_norm = pathSQE_simulations.load_previous_simulation_progress(os.path.join(pathSQE_params['output directory'],'allSym_sim_files/'), user_defined_Qpoints, removed_BZ)
        print('loaded data and norm for sim lens ', len(fully_folded_sims_data), len(fully_folded_sims_norm))

    fold_timer = time.time()
    num_total_slices = pathSQE_helper.count_total_bzfold_slices(pathSQE_params, BZ_list, mtd_spacegroup)
    print(f"\n\n Total expected slices: {num_total_slices}\n\n")
    num_processed_slices = 0
    last_update_step = 0

    for B, BZ_offset in enumerate(BZ_list):
        last_BZ_fold_dsl = []
        dsl_fold = []
        all_slice_names_in_BZ = []
        all_slice_evals = []
        all_sims_in_BZ = [[] for i in range(len(user_defined_Qpoints['path']))]
        for i in range(len(user_defined_Qpoints['path'])):
            path_seg = user_defined_Qpoints['path'][i]
            print('Processing BZ {}, path segments from {} to {}'.format(BZ_offset, path_seg[0], path_seg[1]))

            pt1_init = np.array(user_defined_Qpoints['point_coords'][path_seg[0]])
            pt2_init = np.array(user_defined_Qpoints['point_coords'][path_seg[1]])

            pt1_array, pt2_array = pathSQE_core.generate_unique_paths(mtd_spacegroup, pt1_init, pt2_init, pathSQE_params['primitive to mantid transformation matrix'])
            print("\n{} symmetrically equivalent path segments".format(pt1_array.shape[0]))
            print(np.hstack((pt1_array,pt2_array)))

            all_slice_names_for_pathSeg = []
            slice_evals_for_pathSeg = []
            good_slice_index = None
            for j in range(pt1_array.shape[0]):
                q_dims_and_bins = pathSQE_core.choose_dims_and_bins(pathSQE_params, pt1_array[j], pt2_array[j],  pathSQE_params['u_vec'], pathSQE_params['v_vec'], BZ_offset=BZ_offset)
                slice_desc = pathSQE_helper.make_slice_desc(pathSQE_params, q_dims_and_bins, pt1_array[j], pt2_array[j], path_seg)
                all_slice_names_for_pathSeg.append((slice_desc['Name'],sliceID))
                
                # make slice and evaluate quality based on filters
                slice_utils_07142023.make_slice(mde_data[0], slice_desc, ASCII_slice_folder='', MD_slice_folder='')
                SaveMD(mtd['_data'], Filename=os.path.join(all_slices_dir,'BZ{}_{}2{}_{}to{}_ID{}_data.nxs'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                SaveMD(mtd['_norm'], Filename=os.path.join(all_slices_dir,'BZ{}_{}2{}_{}to{}_ID{}_norm.nxs'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                slice_good = pathSQE_filter_functions.evaluate_slice_quality(slice_name=slice_desc['Name'], filters=pathSQE_params['slice filter functions'], path_seg=path_seg)
                
                # optionally do simulation workflow for analogous slice
                if pathSQE_params['perform simulations']:
                    Qpoints = pathSQE_simulations.construct_sim_Qpts(pt1_array[j], pt2_array[j], q_dims_and_bins[0], pathSQE_params['qdim0 step size'], pathSQE_params['primitive to mantid transformation matrix'], BZ_offset)
                    #print(len(Qpoints))
                    sim_SQE_output = pathSQE_simulations.sim_SQE(pathSQE_params, Qpoints, Temperature=pathSQE_params['T and Ei conditions'][0][0])
                    #print(sim_SQE_output.shape)
                    BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum(pathSQE_params, sim_SQE_output)
                    #print(BinnedSQE.shape)

                    if pathSQE_params['use experimental coverage mask']:
                        # Get the experimental signal array
                        slice_data = np.squeeze(mtd[slice_desc['Name']].getSignalArray())

                        # Ensure shapes match before masking
                        if slice_data.shape == BinnedSQE.shape:
                            BinnedSQE[np.isnan(slice_data)] = np.nan
                        else:
                            print(f"Warning: Shape mismatch between experimental data {slice_data.shape} and simulated data {BinnedSQE.shape}. Masking not applied.")

                    all_sims_in_BZ[i].append(BinnedSQE)
                    np.save(os.path.join(all_slices_dir,'BZ{}_{}2{}_{}to{}_ID{}_sim.npy'.format(BZ_offset, path_seg[0], path_seg[1], pt1_array[j]+BZ_offset, pt2_array[j]+BZ_offset, sliceID)), BinnedSQE)
                    #print('num nans ', np.isnan(BinnedSQE).sum())  # Count NaN values

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
                        
                        # save slice desc of each segment for plotting if this is the recently processed BZ
                        last_BZ_fold_dsl.append(slice_desc.copy())
                        
                    else:
                        # Subsequent good slices, combine data
                        data, norm = pathSQE_core.combine_data_within_bz(data, norm, mtd['_data'], mtd['_norm'])
                else:
                    slice_desc['good_slice'] = False
                    slice_evals_for_pathSeg.append(False)

                sliceID += 1

                # Check if we should print an update
                num_processed_slices += 1
                should_print, last_update_step = pathSQE_helper.should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=2)
                if should_print:
                    pathSQE_helper.print_progress(num_processed_slices, num_total_slices, fold_timer, "BZ Folding")

            # Save the data for the last good slice (if any)
            if good_slice_index is not None:
                saveDir = os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/')
                SaveMD(data, Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/DataAndNorms/{}2{}_BZ{}_data.nxs'.format(path_seg[0], path_seg[1], BZ_offset)), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                SaveMD(norm, Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/DataAndNorms/{}2{}_BZ{}_norm.nxs'.format(path_seg[0], path_seg[1], BZ_offset)), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)

                output_ws = data / norm
                name = '{}2{}_BZ{}_final.nxs'.format(path_seg[0], path_seg[1], BZ_offset)
                SaveMD(output_ws, Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',name), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                LoadMD(Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',name), OutputWorkspace=name, LoadHistory=False)
                
            all_slice_names_in_BZ.append(all_slice_names_for_pathSeg)
            all_slice_evals.append(slice_evals_for_pathSeg)

        if pathSQE_params['perform simulations']:
            folded_sim_results = pathSQE_simulations.fold_symmetric_simulations(all_sims_in_BZ)

            allSym_sim_files = os.path.join(pathSQE_params['output directory'], 'allSym_sim_files/')
            np.savez(os.path.join(allSym_sim_files, 'BZ{}_foldedSim_data.npz'.format(BZ_offset)), *[result[0] for result in folded_sim_results])
            np.savez(os.path.join(allSym_sim_files, 'BZ{}_foldedSim_norm.npz'.format(BZ_offset)), *[result[1] for result in folded_sim_results])
            np.savez(os.path.join(allSym_sim_files, 'BZ{}_foldedSim_sim.npz'.format(BZ_offset)), *[result[2] for result in folded_sim_results])
        
            if len(fully_folded_sims_data)==0:
                #print('first BZ so cloned?')
                fully_folded_sims_data = [results[0] for results in folded_sim_results]
                fully_folded_sims_norm = [results[1] for results in folded_sim_results]
                #print(len(fully_folded_sims_data))
            else:
                #print('folding more BZ segs')
                for tic in range(len(fully_folded_sims_data)):
                    #print('seg ', tic)
                    fully_folded_sims_data[tic] += folded_sim_results[tic][0]
                    fully_folded_sims_norm[tic] += folded_sim_results[tic][1]
            
        print('\nGenerating Report for BZ {}'.format(BZ_offset))
        if pathSQE_params['perform simulations']:
            pathSQE_plotting_and_reports.generate_BZ_Report_with_Sims(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold, folded_sim_results, all_sims_in_BZ)
        else:
            pathSQE_plotting_and_reports.generate_BZ_Report_morePages(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold)

        if not isinstance(pathSQE_params['BZ to process'], list):
            pathSQE_plotting_and_reports.update_BZ_folding_plot(BZ_full_array, BZ_offset, pathSQE_params)

        ########### INTER-BZ folding exp data over all BZ specified ###########
        if len(BZ_list_init)>1: 
            print('\nCombining all {} BZ'.format(len(removed_BZ)+B))
            data_and_norm_dir = os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/DataAndNorms/')
            all_files = os.listdir(data_and_norm_dir)

            last_BZ = BZ_offset
            last_BZ_files = [filename for filename in all_files if str(last_BZ) in filename]
            BZ_list_files = [filename for filename in all_files if filename not in last_BZ_files]

            # generate seg names
            seg_names = []
            for i in range(len(user_defined_Qpoints['path'])):
                seg_names.append('{}2{}'.format(user_defined_Qpoints['path'][i][0], user_defined_Qpoints['path'][i][1]))

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
                SaveMD(folded_ws, Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files','{}_folded.nxs'.format(path_seg)))
                LoadMD(Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files','{}_folded.nxs'.format(path_seg)), OutputWorkspace=path_seg, LoadHistory=False)

            # Adaptive colorbars
            step = float(pathSQE_params['E bins'].split(',')[1])
            ind_aboveElasticLine = int(np.ceil(3 / step))
            pathMax = -1
            pathMin = 1e4

            folded_array = []
            for seg in seg_names:
                folded_array.append(np.squeeze(mtd[seg].getSignalArray()))           
                slice = mtd[seg].getSignalArray()[:, 0, 0, ind_aboveElasticLine:]
                
                # Mask NaN and zero values
                mask = ~np.isnan(slice) & (slice != 0)
                
                try:
                    # Check if there are any valid values before calculating percentiles
                    if np.any(mask):
                        segMin = np.percentile(slice[mask], 10)
                        segMax = np.percentile(slice[mask], 95)

                        # Update pathMin and pathMax
                        if segMax > pathMax:
                            pathMax = segMax
                        if segMin < pathMin:
                            pathMin = segMin
                    else:
                        print(f"Warning: No valid data in segment {seg}, skipping percentile calculation.")
                
                except Exception as e:
                    print(f"Error processing segment {seg}: {e}")
                    
            np.save(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_path_plot_{}.npy'.format(B+len(removed_BZ))), np.vstack(folded_array))
            fig = pathSQE_plotting_and_reports.plot_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin, vma=pathMax, cma=pathSQE_params['cmap'])
            fig.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_path_plot_{}.png'.format(B+len(removed_BZ))))
            fig2 = pathSQE_plotting_and_reports.plot_SED_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=last_BZ_fold_dsl, pathSQE_params=pathSQE_params)
            fig2.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_SED_path_plot_{}.png'.format(B+len(removed_BZ))))

            if pathSQE_params['perform simulations']: 
                folded_sim_final = []
                for ps in range(len(fully_folded_sims_data)):
                    final_norm = fully_folded_sims_norm[ps].copy()
                    final_norm[final_norm == 0] = np.nan
                    folded_sim_final.append( fully_folded_sims_data[ps] / final_norm )

                simMax = -1
                simMin = 1e4
                # Find min/max percentiles across all slices
                for path_sim in folded_sim_final:
                    mask = ~np.isnan(path_sim) & (path_sim > 0)
                    if np.any(mask):
                        segMin = np.percentile(path_sim[mask],60)
                        segMax = np.max(path_sim[mask])
                        simMax = max(pathMax, segMax)
                        simMin = min(pathMin, segMin)

                fig3 = pathSQE_plotting_and_reports.plot_along_path_foldedBZ_sim(folded_sim_final,last_BZ_fold_dsl,pathSQE_params['T and Ei conditions'][0][1],simMin,simMax,pathSQE_params['cmap'])
                fig3.savefig(os.path.join(pathSQE_params['output directory'],'folding_progess/','folded_sim_path_plot_{}.png'.format(B+len(removed_BZ))))
                np.savez(os.path.join(pathSQE_params['output directory'], 'folding_progess/', 'folded_sim_path_plot_{}.npz'.format(B+len(removed_BZ))), *folded_sim_final)

            print('\nBZ {} done. Total time elapsed {}'.format(BZ_offset, pathSQE_helper.format_time(time.time()-start_time)))

    print('\nFolded {} BZ in {}'.format(len(BZ_list_init), pathSQE_helper.format_time(time.time()-start_time)))


#######################################################################################################
        #################################### all symmetric 1d cuts ##################################
#######################################################################################################

if pathSQE_params['all symmetric 1d cuts']:
    sym_points_timer = time.time()
    num_total_slices = pathSQE_helper.count_total_symPt_slices(pathSQE_params, BZ_list, mtd_spacegroup, len(mde_data))
    print(f"Total expected slices: {num_total_slices}")
    num_processed_slices = 0
    last_update_step = 0

    for i in range(len(user_defined_Qpoints['1d_points'])):
        symPt_slice_names_allTemp = [[] for _ in range(len(mde_data))]
        symPt_slice_arrays = []
        
        # find unique sym points (arbitrary BZ)
        symPoint = user_defined_Qpoints['1d_points'][i]
        symPt_init = np.array(user_defined_Qpoints['point_coords'][symPoint])
        symPts_array = pathSQE_core.generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init)
        print('\nPoint', symPoint)
        print('Equivalent points in a single BZ: ', symPts_array)

        # find all unique sym points in BZ coverage (avoid repeating same point in different BZ)
        allSymPts_array = pathSQE_core.generate_unique_symPts_inAllBZ(BZ_list, symPts_array)
        np.save(os.path.join(pathSQE_params['output directory'], 'probed_{}_pts.npy'.format(symPoint)), allSymPts_array)

        for j in range(allSymPts_array.shape[0]):
            for k in range(len(mde_data)):
                if k == 0:
                    slice_desc = pathSQE_helper.make_slice_desc_1DSymPoints(pathSQE_params, symPoint, allSymPts_array[j])
                    slice_desc_temp2og = slice_desc.copy()
                    slice_desc['Name'] = slice_desc['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                    slice_utils_07142023.make_slice(mde_data[k], slice_desc, ASCII_slice_folder='', MD_slice_folder='') 
                    
                    if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:
                        symPt_slice_arrays.append(allSymPts_array[j])
                        symPt_slice_names_allTemp[k].append(slice_desc['Name'])
                        SaveMD(mtd[slice_desc['Name']], Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',slice_desc['Name']+'.nxs'), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
            
                else:
                    # for other temps
                    slice_desc_temp2 = slice_desc_temp2og.copy()
                    slice_desc_temp2['Name'] = slice_desc_temp2['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                    slice_utils_07142023.make_slice(mde_data[k], slice_desc_temp2, ASCII_slice_folder='', MD_slice_folder='') 
                
                    if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:     
                        symPt_slice_names_allTemp[k].append(slice_desc_temp2['Name'])
                        SaveMD(mtd[slice_desc_temp2['Name']], Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',slice_desc_temp2['Name']+'.nxs'), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
                          
                # Check if we should print an update
                num_processed_slices += 1
                should_print, last_update_step = pathSQE_helper.should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=2)
                if should_print:
                    pathSQE_helper.print_progress(num_processed_slices, num_total_slices, sym_points_timer, "1D Sym Point")


        #print('sym point slice names: ', symPt_slice_names_allTemp)
        print('Generating Report for {} point'.format(symPoint))
        pathSQE_plotting_and_reports.generate_BZ_Report_1DSymPoints(pathSQE_params, symPoint, symPt_slice_arrays, symPt_slice_names_allTemp, pathSQE_params['u_vec'], pathSQE_params['v_vec'])


#######################################################################################################
        ################################ SIMPLE 2D PATH ####################################
#######################################################################################################

if not pathSQE_params['all symmetric 2d slices'] and user_defined_Qpoints['path'] and len(user_defined_Qpoints['path']) > 0:

    dsl=[]
    seg_names = []
    all_sims = []
    for path_seg in user_defined_Qpoints['path']:
        pt1 = np.array(user_defined_Qpoints['point_coords'][path_seg[0]])
        pt2 = np.array(user_defined_Qpoints['point_coords'][path_seg[1]])
        q_dims_and_bins = pathSQE_core.choose_dims_and_bins(pathSQE_params, pt1, pt2, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  
        slice_desc = pathSQE_helper.make_slice_desc(pathSQE_params, q_dims_and_bins, pt1, pt2, path_seg)
        dsl.append(slice_desc)
        seg_names.append(slice_desc['Name'])
        slice_utils_07142023.make_slice(mde_data[0], slice_desc, ASCII_slice_folder='', MD_slice_folder='')
        SaveMD(slice_desc['Name'], Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/{}_to_{}.nxs'.format(path_seg[0], path_seg[1])), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)

        # optionally do simulation workflow for analogous slice
        if pathSQE_params['perform simulations']:
            Qpoints = pathSQE_simulations.construct_sim_Qpts(pt1, pt2, q_dims_and_bins[0], pathSQE_params['qdim0 step size'], pathSQE_params['primitive to mantid transformation matrix'])
            sim_SQE_output = pathSQE_simulations.sim_SQE(pathSQE_params, Qpoints, Temperature=pathSQE_params['T and Ei conditions'][0][0])
            BinnedSQE = pathSQE_simulations.SQE_to_2d_spectrum(pathSQE_params, sim_SQE_output)

            if pathSQE_params['use experimental coverage mask']:
                # Get the experimental signal array
                slice_data = np.squeeze(mtd[slice_desc['Name']].getSignalArray())

                # Ensure shapes match before masking
                if slice_data.shape == BinnedSQE.shape:
                    BinnedSQE[np.isnan(slice_data)] = np.nan
                else:
                    print(f"Warning: Shape mismatch between experimental data {slice_data.shape} and simulated data {BinnedSQE.shape}. Masking not applied.")

            all_sims.append(BinnedSQE)

    # Adaptive colorbars
    step = float(pathSQE_params['E bins'].split(',')[1])
    ind_aboveElasticLine = int(np.ceil(3 / step))
    pathMax = -1
    pathMin = 1e4

    combined_array = []
    for seg in seg_names:
        combined_array.append(np.squeeze(mtd[seg].getSignalArray()))           
        slice = mtd[seg].getSignalArray()[:, 0, 0, ind_aboveElasticLine:]
        
        # Mask NaN and zero values
        mask = ~np.isnan(slice) & (slice != 0)
        
        try:
            # Check if there are any valid values before calculating percentiles
            if np.any(mask):
                segMin = np.percentile(slice[mask], 10)
                segMax = np.percentile(slice[mask], 95)

                # Update pathMin and pathMax
                if segMax > pathMax:
                    pathMax = segMax
                if segMin < pathMin:
                    pathMin = segMin
            else:
                print(f"Warning: No valid data in segment {seg}, skipping percentile calculation.")
        
        except Exception as e:
            print(f"Error processing segment {seg}: {e}")
            
    fig = pathSQE_plotting_and_reports.plot_along_path_foldedBZ(foldedSegNames=seg_names, dsl_fold=dsl, Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin, vma=pathMax, cma=pathSQE_params['cmap'])
    fig.savefig(os.path.join(pathSQE_params['output directory'],'exp_path_plot.png'))
    np.savez(os.path.join(pathSQE_params['output directory'],'exp_path_plot.npz'), *combined_array)

    if pathSQE_params['perform simulations']:
        simMax = -1
        simMin = 1e4
        # Find min/max percentiles across all slices
        for path_sim in all_sims:
            mask = ~np.isnan(path_sim) & (path_sim > 0)
            if np.any(mask):
                segMin = np.percentile(path_sim[mask],60)
                segMax = np.max(path_sim[mask])
                simMax = max(pathMax, segMax)
                simMin = min(pathMin, segMin)

        fig3 = pathSQE_plotting_and_reports.plot_along_path_foldedBZ_sim(all_sims,dsl,pathSQE_params['T and Ei conditions'][0][1],simMin,simMax,pathSQE_params['cmap'])
        fig3.savefig(os.path.join(pathSQE_params['output directory'],'sim_path_plot.png'))
        np.savez(os.path.join(pathSQE_params['output directory'],'sim_path_plot.npz'), *all_sims)


#######################################################################################################
        ################################ SIMPLE 1D POINTS ####################################
#######################################################################################################

if not pathSQE_params['all symmetric 1d cuts'] and user_defined_Qpoints['1d_points'] and len(user_defined_Qpoints['1d_points']) > 0:

    for i in range(len(user_defined_Qpoints['1d_points'])):
        slice_names = []
        
        # find unique sym points (arbitrary BZ)
        symPoint = user_defined_Qpoints['1d_points'][i]
        symPt_init = np.array(user_defined_Qpoints['point_coords'][symPoint])
        print('\nPoint', symPoint, symPt_init)

        for k in range(len(mde_data)):
            if k == 0:
                slice_desc = pathSQE_helper.make_slice_desc_1DSymPoints(pathSQE_params, symPoint, symPt_init)
                slice_desc_temp2og = slice_desc.copy()
                slice_desc['Name'] = slice_desc['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                slice_names.append(slice_desc['Name'])
                slice_utils_07142023.make_slice(mde_data[k], slice_desc, ASCII_slice_folder='', MD_slice_folder='') 
                
                if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:
                    SaveMD(slice_desc['Name'], Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',slice_desc['Name']+'.nxs'), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
        
            else:
                # for other temps
                slice_desc_temp2 = slice_desc_temp2og.copy()
                slice_desc_temp2['Name'] = slice_desc_temp2['Name']+'_{}K'.format(mde_data[k]['SampleLogVariables']['Temperature'])
                slice_names.append(slice_desc_temp2['Name'])
                slice_utils_07142023.make_slice(mde_data[k], slice_desc_temp2, ASCII_slice_folder='', MD_slice_folder='') 
            
                if np.count_nonzero((mtd[slice_desc['Name']].getSignalArray() != 0) & (~np.isnan(mtd[slice_desc['Name']].getSignalArray()))) > 5:     
                    SaveMD(slice_desc_temp2['Name'], Filename=os.path.join(pathSQE_params['output directory'],'allSym_nxs_files/',slice_desc_temp2['Name']+'.nxs'), SaveHistory=False, SaveInstrument=False, SaveSample=False, SaveLogs=False)
        
        pathSQE_plotting_and_reports.plot_simple_1d_pt(pathSQE_params, symPoint, symPt_init, slice_names)




# Remember to close the output file when done
time.sleep(5)
output_file.close()
sys.stdout = sys.__stdout__  # Reset stdout before the script exits
