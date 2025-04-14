import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.cm import ScalarMappable, coolwarm, viridis
from matplotlib.backends.backend_pdf import PdfPages
from mantid.simpleapi import *
from mantid.geometry import SpaceGroupFactory, PointGroupFactory
from slice_utils_07142023 import make_slice
from mantid import plots
import time
import re


########################################################################################################
# All necessary functions for running pathSQE
# Author: Aiden Sable. April 2025.
########################################################################################################

# To run pathSQE_driver.py, you need the following files in addition to this script:
#       pathSQE_input.py
#       pathSQE_driver.py
#       define_data.py
#       reduce_data_to_MDE.py
#       slice_utils.py
#       POSCAR (if 'use seeKpath path'=True or 'perform simulations=True)
#       FORCE_CONSTANTS (if 'perform simulations=True)



##############################################################################################################################
################################################ FUNCTIONS FOR pathSQE DRIVER ################################################
##############################################################################################################################


def simple_read_poscar(fname):
    """Read a POSCAR file."""
    with open(fname) as f:
        lines = [l.partition("!")[0] for l in f.readlines()]

    alat = float(lines[1])
    v1 = [float(_) * alat for _ in lines[2].split()]
    v2 = [float(_) * alat for _ in lines[3].split()]
    v3 = [float(_) * alat for _ in lines[4].split()]
    cell = [v1, v2, v3]

    species = lines[5].split()
    num_atoms = [int(_) for _ in lines[6].split()]

    next_line = lines[7]
    if next_line.strip().lower() != "direct":
        raise ValueError("This simple routine can only deal with 'direct' POSCARs")
    # Note: to support also cartesian, remember to multiply the coordinates by alat

    positions = []
    atomic_numbers = []
    cnt = 8
    for el, num in zip(species, num_atoms):
        atom_num = atoms_num_dict[el.capitalize()]
        for _ in range(num):
            atomic_numbers.append(atom_num)
            positions.append([float(_) for _ in lines[cnt].split()])
            cnt += 1

    return (cell, positions, atomic_numbers)

'''
def find_qdim_1and2(qdim0):
    hkl = np.array([[1,0,0],[0,1,0],[0,0,1]])
    
    # if perp to path method is wanted
    if np.count_nonzero(qdim0) == 1:
        # if diff=dim0 is 1D, then dim1 and dim2 will also just be 1D h,k,or l
        qdim0_dir = np.nonzero(qdim0)[0]
        hkl_avail = np.delete(hkl, qdim0_dir, axis=0)
        qdim1 = hkl_avail[0]
        qdim2 = hkl_avail[1]
    elif (np.count_nonzero(qdim0) == 2) or (np.count_nonzero(qdim0) == 3):
        nzeroInds = np.nonzero(qdim0)[0]
        qdim1 = qdim0.copy()
        # swap h and k and make one negative to give perpen direction within the horizontal plane
        # If qdim0 is 3D, need to set l to 0 so that it is perpen and in the horizontal plane 
        qdim1[nzeroInds[0]], qdim1[nzeroInds[1]] = qdim0[nzeroInds[1]], -qdim0[nzeroInds[0]]
        if np.count_nonzero(qdim0) == 3:
            qdim1[2] = 0 
        # qdim2 is simply the direction perpen to qdim0 and qdim1
        qdim2 = np.cross(qdim0, qdim1)
        qdim2 = qdim2 / np.max(np.absolute(qdim2))
        # MAYBE DON'T (messing w bin offset) - normalize the two qdims so that mag of each is 1 A-1 and can easily do bins +/- 0.15
        #qdim1 = qdim1 / np.linalg.norm(qdim1)
        #qdim2 = qdim2 / np.linalg.norm(qdim2)
    
    return qdim1, qdim2
    '''



def find_qdim_1and2(qdim0, u, v):
    """
    Find qdim1 and qdim2 given qdim0, ensuring qdim1 lies in the plane defined by u and v.

    Parameters:
        qdim0 (array): The difference vector between input Q points, defining the path direction.
        u (array): A vector defining the first basis direction of the horizontal scattering plane.
        v (array): A vector defining the second basis direction of the horizontal scattering plane.

    Returns:
        qdim1 (array): A vector perpendicular to qdim0 and lying in the plane defined by (u, v).
        qdim2 (array): A vector perpendicular to both qdim0 and qdim1.
    """
    # Normalize input vectors
    qdim0 = qdim0.astype(float) / np.linalg.norm(qdim0)
    u = u.astype(float) / np.linalg.norm(u)
    v = v.astype(float) / np.linalg.norm(v)

    # Compute the normal to the scattering plane
    z_scattering = np.cross(u, v)
    z_scattering /= np.linalg.norm(z_scattering)

    # Attempt to construct qdim1 as the cross product of qdim0 and z_scattering
    qdim1 = np.cross(z_scattering, qdim0)

    if np.linalg.norm(qdim1) < 1e-10:  # qdim0 is collinear with z_scattering
        # Find a new qdim1 as a linear combination of u and v, ensuring perpendicularity to qdim0
        qdim1 = u - np.dot(u, qdim0) * qdim0  # Start with u and remove its component along qdim0

        # If still degenerate, use v instead
        if np.linalg.norm(qdim1) < 1e-10:
            qdim1 = v - np.dot(v, qdim0) * qdim0

        # Normalize qdim1
        if np.linalg.norm(qdim1) < 1e-10:
            raise ValueError("Failed to construct qdim1. Check input vectors for collinearity issues.")
    
    qdim1 /= np.linalg.norm(qdim1)

    # Compute qdim2 as the cross product (ensuring it is perpendicular to both qdim0 and qdim1)
    qdim2 = np.cross(qdim0, qdim1)

    # Normalize qdim2 safely
    if np.linalg.norm(qdim2) > 1e-10:
        qdim2 /= np.linalg.norm(qdim2)
    else:
        raise ValueError("Failed to construct qdim2. Check for collinearity issues.")

    return qdim1, qdim2



def find_qbins(pathSQE_params, qdim0, qdim1, qdim2, pt1, pt2, BZ_offset):
    # we start with bins of 0 so hopefully can catch it if code isn't working
    qdim1_range = 0
    qdim2_range = 0 

    # find the new path endpoint coords given change of basis from cart to qdim0,1,2
    CoB_mat = np.linalg.inv(np.array([qdim0,qdim1,qdim2]).T) # back to og, edited from Bi version with np.array([qdim0,qdim1,qdim2])
    path_start_transf = CoB_mat @ (pt1 + BZ_offset)
    path_end_transf = CoB_mat @ (pt2 + BZ_offset)
    print('Slicing ', pt1+BZ_offset, ' to ', pt2+BZ_offset)
    #print('find_qbins start and end: ', path_start_transf, path_end_transf)
    
    # Compute absolute step size based on fractional step // commenting bc not sure this is what we want..
    #qdim0_step_size = np.abs((path_end_transf[0] - path_start_transf[0]) * pathSQE_params['qdim0 step size'])

    # Construct the array with the computed step size
    qdim0_range = np.array([path_start_transf[0], pathSQE_params['qdim0 step size'], path_end_transf[0]])
    
    # might need to check bin/path endpoints...
    # old check_bin0_and_path_endpoints function
    
    # fixed integration range for qdim1 and qdim2
    qdim1_int_range = pathSQE_params['qdim1 integration range']
    qdim1_range = np.array([path_start_transf[1]-qdim1_int_range, path_end_transf[1]+qdim1_int_range])
    
    qdim2_int_range = pathSQE_params['qdim2 integration range']
    qdim2_range = np.array([path_start_transf[2]-qdim2_int_range, path_end_transf[2]+qdim2_int_range])
    
    # if qdim1 or qdim2 bin range gets too big raise error because code above may be wrong...
    #if (np.absolute(qdim1_range[1]-qdim1_range[0]) > 0.2) or (np.absolute(qdim2_range[1]-qdim2_range[0]) > 0.2):
    #    raise ValueError("Might be something weird with bins for qdim1 or dqim2")
    
    return qdim0_range, qdim1_range, qdim2_range



def choose_dims_and_bins(pathSQE_params, point1, point2, u, v, BZ_offset=np.array([0,0,0])):
    diff = point2 - point1

    # calculating default QDimension0 axis direction and scaling
    qdim0 = diff / np.max(np.absolute(diff))

    # determine the directions of qdim1 and qdim2
    qdim1, qdim2 = find_qdim_1and2(qdim0, u, v)
    #print('qdims ', qdim0, qdim1, qdim2)
    qdim0_range, qdim1_range, qdim2_range = find_qbins(pathSQE_params, qdim0, qdim1, qdim2, point1, point2, BZ_offset)

    # return list contains various types (mainly np arrays) that are converted to properly formatted strings
    # as a part of the slice description generation process
    q_dims_and_bins = [qdim0, qdim1, qdim2, qdim0_range, qdim1_range, qdim2_range]
    return q_dims_and_bins



def conv_to_desc_string(input_item):
    # works for np arrays, lists, and tuples
    string = ''
    for item in input_item:
        string = string + str(item) +','
    string = string[:-1]
    
    return string



def make_slice_desc(pathSQE_params, q_dims_and_bins, point1, point2, path_seg):
    diff = point2 - point1

    # slice description with desired dims and bins and unique name
    slice_desc={'QDimension0':conv_to_desc_string(q_dims_and_bins[0]),
                'QDimension1':conv_to_desc_string(q_dims_and_bins[1]),
                'QDimension2':conv_to_desc_string(q_dims_and_bins[2]),
                'Dimension0Name':'QDimension0',
                'Dimension0Binning':conv_to_desc_string(q_dims_and_bins[3]),
                'Dimension1Name':'QDimension1',
                'Dimension1Binning':conv_to_desc_string(q_dims_and_bins[4]),
                'Dimension2Name':'QDimension2',
                'Dimension2Binning':conv_to_desc_string(q_dims_and_bins[5]),
                'Dimension3Name':'DeltaE',
                'Dimension3Binning':pathSQE_params['E bins'],
                'SymmetryOperations':'x,y,z',
                'Name':'pathSQE_'+conv_to_desc_string(point1)+'_to_'+conv_to_desc_string(point2),
                'qdim0_range':q_dims_and_bins[3],
                'seg_start_name':path_seg[0],
                'seg_end_name':path_seg[1],
                'inv_angstrom_ratio':np.linalg.norm(diff)/0.5,
                'good_slice':True}

    #slice_desc['Name'] = slice_desc['Name'].replace(',', '_').replace('.', '').replace('-','m')
    
    return slice_desc



def find_matching_spacegroup(pathSQE_params):
    query_spacegroup = pathSQE_params['space group']
    print(f'Searching for spacegroup {query_spacegroup}')
    spacegroup_list = mantid.geometry.SpaceGroupFactory.getAllSpaceGroupSymbols()
    query_spacegroup = ''.join(query_spacegroup.split())  # Remove spacing from the query_spacegroup

    for spacegroup in spacegroup_list:
        stripped_spacegroup = ''.join(spacegroup.split())  # Remove spacing from the spacegroup in the list
        if query_spacegroup == stripped_spacegroup:
            print(f"Matching spacegroup found: {spacegroup}")
            print(SpaceGroupFactory.createSpaceGroup(spacegroup))
            print(SpaceGroupFactory.createSpaceGroup(spacegroup).getPointGroup())
            return spacegroup

    raise ValueError("No matching spacegroup found.")



def generate_unique_paths(mtd_spacegroup, pt1, pt2, prim2mantid):
    # Get the symmetry operations for the spacegroup
    spacegroup = SpaceGroupFactory.createSpaceGroup(mtd_spacegroup)
    pg = spacegroup.getPointGroup()
    #print(pg)
    
    rotations = [np.array([so.transformHKL((1, 0, 0)), so.transformHKL((0, 1, 0)), so.transformHKL((0, 0, 1))])
                 for so in pg.getSymmetryOperations()]

    new_pt1s = rotations @ pt1 #@ prim2mantid
    new_pt2s = rotations @ pt2 #@ prim2mantid

    # Combine new_pt1s and new_pt2s into a single array for comparison
    combined_array = np.column_stack((new_pt1s, new_pt2s))

    # Find unique rows in the combined array
    unique_combined_array = np.unique(combined_array, axis=0)

    # Split the unique_combined_array back into separate arrays for pt1 and pt2
    unique_pt1_array = unique_combined_array[:, :3]
    unique_pt2_array = unique_combined_array[:, 3:]

    return unique_pt1_array, unique_pt2_array



def combine_data_within_bz(data, norm, mtd_data, mtd_norm):
    data.setSignalArray(data.getSignalArray() + mtd_data.getSignalArray())
    norm.setSignalArray(norm.getSignalArray() + mtd_norm.getSignalArray())
    data.setErrorSquaredArray(data.getErrorSquaredArray() + mtd_data.getErrorSquaredArray())
    return data, norm



def get_next_sliceID(all_slices_dir):
    """
    Searches the given directory for files matching the slice pattern and 
    returns the next available unique slice ID (max existing + 1).
    """
    slice_pattern = re.compile(r'_ID(\d+)_data\.nxs$')
    existing_ids = []

    if not os.path.isdir(all_slices_dir):
        os.makedirs(all_slices_dir)
        return 0

    for fname in os.listdir(all_slices_dir):
        match = slice_pattern.search(fname)
        if match:
            existing_ids.append(int(match.group(1)))

    if existing_ids:
        return max(existing_ids) + 1
    else:
        return 0



def project_to_uv_plane(BZ_list, u, v):
    """
    Projects BZ centers onto the (u, v) plane.
    
    Parameters:
        BZ_list: List of tuples, where each tuple contains (BZ center, coverage)
        u: First basis vector defining the plane
        v: Second basis vector defining the plane
        
    Returns:
        np.array of projected (u, v) coordinates.
    """
    u = np.array(u) / np.linalg.norm(u)  # Normalize u
    v = np.array(v) / np.linalg.norm(v)  # Normalize v

    uv_matrix = np.vstack([u, v]).T  # Form transformation matrix
    uv_inv = np.linalg.pinv(uv_matrix)  # Compute pseudo-inverse for projection

    BZ_uv = np.array([uv_inv @ bz[0] for bz in BZ_list])  # Project each BZ center
    return BZ_uv



def find_BZ_with_data(mde_data, pathSQE_params, hkl_range=(-10, 10, -10, 10, -10, 10)): 
    """
    Identify and rank Brillouin Zones (BZs) by fractional data coverage.

    Parameters:
        mde_data: Mantid workspace containing the data.
        pathSQE_params (dict): Dictionary containing parameters, including the transformation matrix.
        hkl_range (tuple): Defines BZ grid in primitive basis.

    Returns:
        list: Ranked list of BZ centers sorted by fractional coverage (descending).
    """
    E_bins = pathSQE_params['E bins'].split(',')
    min_fraction = pathSQE_params['BZ to process']

    # Transform to current basis
    prim2current = pathSQE_params['primtive to mantid T matrix'].copy().T

    # Generate all integer BZ centers in primitive basis
    Q_primitive = np.array([[h, k, l] for h in range(hkl_range[0], hkl_range[1] + 1)
                                      for k in range(hkl_range[2], hkl_range[3] + 1)
                                      for l in range(hkl_range[4], hkl_range[5] + 1)])

    #print(prim2current[0], prim2current[1], prim2current[2])

    slice_desc = {
        'QDimension0': conv_to_desc_string(prim2current[0]),  # make it current2prim columns instead of simple HKL aligned basis
        'QDimension1': conv_to_desc_string(prim2current[1]),
        'QDimension2': conv_to_desc_string(prim2current[2]),
        'Dimension0Name': 'QDimension0',
        'Dimension0Binning': f'{hkl_range[0] - 0.5},0.5,{hkl_range[1] + 0.5}',  # Step size = 0.5 for octants
        'Dimension1Name': 'QDimension1',
        'Dimension1Binning': f'{hkl_range[2] - 0.5},0.5,{hkl_range[3] + 0.5}',
        'Dimension2Name': 'QDimension2',
        'Dimension2Binning': f'{hkl_range[4] - 0.5},0.5,{hkl_range[5] + 0.5}',
        'Dimension3Name': 'DeltaE',
        'Dimension3Binning': f'{E_bins[0]},{0.2*(float(E_bins[2])-float(E_bins[0]))},{E_bins[2]}',
        'Name': 'pathSQE_BZwithData'
    }

    make_slice(mde_data, slice_desc)
    data_slice = mtd[slice_desc['Name']].getSignalArray()
    #print('cov slice shape', data_slice.shape)
    #np.save(os.path.join(pathSQE_params['output directory'],'BZ_coverage_slice.npy'), data_slice)

    # Step 3: Compute BZ fractional coverage
    BZ_list = []
    bin_edges_h = np.arange(hkl_range[0] - 0.5, hkl_range[1] + 1.0, 0.5)
    bin_edges_k = np.arange(hkl_range[2] - 0.5, hkl_range[3] + 1.0, 0.5)
    bin_edges_l = np.arange(hkl_range[4] - 0.5, hkl_range[5] + 1.0, 0.5)

    BZ_full_list = []  # Store all BZ centers with their fractional coverages

    for q_center in Q_primitive:
        Q_current= q_center @ prim2current.T
        h_idx = np.digitize(q_center[0] + np.array([-0.25, 0.25]), bin_edges_h) - 1
        k_idx = np.digitize(q_center[1] + np.array([-0.25, 0.25]), bin_edges_k) - 1
        l_idx = np.digitize(q_center[2] + np.array([-0.25, 0.25]), bin_edges_l) - 1

        num_bins = 0
        num_valid_bins = 0

        for hi in h_idx:
            for ki in k_idx:
                for li in l_idx:
                    octant_data = data_slice[hi, ki, li, :]
                    num_bins += np.size(octant_data)
                    num_valid_bins += np.sum(np.logical_and(~np.isnan(octant_data), octant_data != 0))

        # Compute the fractional coverage
        fraction_coverage = num_valid_bins / num_bins if num_bins > 0 else 0

        # Only include BZs with coverage > 0
        if fraction_coverage > 0:
            BZ_full_list.append((Q_current, fraction_coverage))

        # Keep only those meeting the threshold for later analysis
        if fraction_coverage >= min_fraction:
            BZ_list.append((Q_current, fraction_coverage))

    # Convert to NumPy array and sort (highest to lowest coverage)
    BZ_full_array = np.array(BZ_full_list, dtype=object)  # Keep object type for tuple storage
    BZ_full_array = BZ_full_array[BZ_full_array[:, 1].argsort()[::-1]]  # Sort by coverage (descending)

    # Save to file
    np.savetxt(os.path.join(pathSQE_params['output directory'],"BZ_coverage_sorted.txt"), 
               np.column_stack([BZ_full_array[:, 0], BZ_full_array[:, 1]]), 
               fmt='%s', header="BZ_center   Fraction_Coverage")

    # --- Scatter Plot: Coverage vs. Rank ---
    plt.figure(figsize=(6, 4))
    plt.scatter(range(len(BZ_full_array)), BZ_full_array[:, 1], color='tab:blue', s=10, label="BZ Coverage")
    plt.axhline(y=min_fraction, color='tab:grey', linestyle='--', label=f"Threshold ({min_fraction:.2f})")
    plt.xlabel("BZ Rank")
    plt.ylabel("Fractional Coverage")
    plt.title("BZ Coverage vs. Rank")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(pathSQE_params['output directory'],'BZ_coverage.png'))
    plt.close()

    # --- 2D Scatter Plot: BZs with Data in u,v plane ---
    # Project BZ Centers to (u, v) Plane
    BZ_uv_all = project_to_uv_plane(BZ_full_list, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  # All BZs with data
    BZ_uv_threshold = project_to_uv_plane(BZ_list, pathSQE_params['u_vec'], pathSQE_params['v_vec'])  # BZs above coverage threshold
    np.save(os.path.join(pathSQE_params['output directory'], 'BZ_uv_all.npy'), BZ_uv_all)
    np.save(os.path.join(pathSQE_params['output directory'], 'BZ_uv_threshold.npy'), BZ_uv_threshold)

    plt.figure(figsize=(6, 6))
    plt.scatter(BZ_uv_all[:, 0], BZ_uv_all[:, 1], color='black', s=10, label="BZ with Data")
    plt.scatter(BZ_uv_threshold[:, 0], BZ_uv_threshold[:, 1], color='red', s=20, label=f"BZ > {min_fraction:.2f} Coverage")

    plt.xlabel("{}".format(pathSQE_params['u_vec']))
    plt.ylabel("{}".format(pathSQE_params['v_vec']))
    plt.title("BZ Centers with Data (Projected to u-v Plane)")
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.savefig(os.path.join(pathSQE_params['output directory'], 'BZ_uv_scatter.png'))
    plt.close()

    # Step 4: Sort and return the ranked BZ list
    BZ_list.sort(key=lambda x: x[1], reverse=True)  # Sort descending by coverage

    # Extract separate lists for BZ centers and their corresponding fractional coverages
    BZ_centers = [bz[0] for bz in BZ_list]
    fractional_coverages = [bz[1] for bz in BZ_list]

    return BZ_centers, fractional_coverages, BZ_full_array 



def update_BZ_folding_plot(BZ_full_array, BZ_offset, pathSQE_params):
    """
    Updates and saves the BZ Coverage scatter plot during iterative folding.
    
    Parameters:
        BZ_full_list (list): List of (BZ_center, Fractional_Coverage) tuples.
        BZ_offset (array): The BZ center that was just processed.
        pathSQE_params (dict): Dictionary containing 'output directory' path.
        min_fraction (float): Minimum fraction threshold line.
    """
    # Find the index of the last processed BZ (BZ_offset) in the sorted array
    min_fraction = pathSQE_params['BZ to process']
    BZ_centers = np.vstack(BZ_full_array[:, 0])  # Extract BZ center coordinates
    match_idx = np.where(np.all(BZ_centers == BZ_offset, axis=1))[0]

    if match_idx.size == 0:
        print("Warning: BZ_offset not found in BZ_full_array. It will not be highlighted.")
        match_idx = None
    else:
        match_idx = match_idx[0]  # Convert to integer index

    # --- Scatter Plot: Coverage vs. Rank ---
    plt.figure(figsize=(6, 4))
    plt.scatter(range(len(BZ_full_array)), BZ_full_array[:, 1], color='tab:blue', s=10, label="BZ Coverage")
    plt.axhline(y=min_fraction, color='tab:grey', linestyle='--', label=f"Threshold ({min_fraction:.2f})")

    # Mark the last processed BZ with a red star
    if match_idx is not None:
        plt.scatter(match_idx, BZ_full_array[match_idx, 1], color='red', marker='*', s=30, label="Last Processed BZ")

    # Labels and title
    plt.xlabel("BZ Rank")
    plt.ylabel("Fractional Coverage")
    plt.title("BZ Coverage vs. Rank")
    plt.legend()
    plt.grid(True)

    # Save the updated plot
    plt.savefig(os.path.join(pathSQE_params['output directory'], 'BZ_coverage_iteration.png'))
    plt.close()


import os
import re
import numpy as np

def get_existing_BZ_counts(output_folder_path):
    """
    Retrieve BZ arrays that already have processed .nxs files and count occurrences.

    Parameters:
        output_folder_path (str): Path to the folder containing processed BZ files.

    Returns:
        dict: Mapping of tuple(BZ center) -> count of occurrences.
    """
    BZ_counts = {}
    pattern = re.compile(r'_BZ\[(.*?)\]_final\.nxs$')  # Regex to capture the BZ array

    for filename in os.listdir(output_folder_path):
        match = pattern.search(filename)
        if match:
            BZ_str = match.group(1)  # Extract array content inside brackets
            try:
                BZ_array = tuple(float(x) for x in BZ_str.split())  # Convert to tuple for dict keys
                BZ_counts[BZ_array] = BZ_counts.get(BZ_array, 0) + 1
            except ValueError:
                print(f"Warning: Could not parse BZ array from filename '{filename}'")

    return BZ_counts


def resume_analysis(BZ_list, output_folder_path, num_paths):
    """
    Resume analysis by removing already processed BZ arrays from BZ_list, 
    ensuring that a BZ is only removed when all paths have been processed.

    Parameters:
        BZ_list (np.ndarray): List of BZ centers to process.
        output_folder_path (str): Path to the folder containing processed BZ files.
        num_paths (int): Number of paths each BZ should have to be considered complete.

    Returns:
        new_BZ_list (np.ndarray): Filtered list of BZs still needing processing.
        removed_BZ (np.ndarray): List of BZs that were fully processed.
    """
    existing_BZ_counts = get_existing_BZ_counts(output_folder_path)

    # Check for each BZ if all paths are processed
    mask = np.array([
        existing_BZ_counts.get(tuple(bz), 0) < num_paths  # Keep only if less than required paths
        for bz in BZ_list
    ])
    
    new_BZ_list = BZ_list[mask]
    removed_BZ = BZ_list[~mask]

    print(f"{len(removed_BZ)} BZ arrays fully processed and removed:")
    print(removed_BZ.tolist())

    return new_BZ_list, removed_BZ



def count_total_bzfold_slices(pathSQE_params, BZ_list, mtd_spacegroup):
    """
    Calculates the total number of slices expected before processing starts.

    Parameters:
    - BZ_list: list of BZ offsets
    - user_defined_Qpoints: dictionary containing 'path' and 'point_coords'
    - mtd_spacegroup: Mantid space group object used for path generation
    - pathSQE_params: dictionary containing parameters including 'primtive to mantid T matrix'

    Returns:
    - total_slices: int, total number of slices to be processed
    """

    num_BZ = len(BZ_list)  # Number of BZs
    unique_slices_per_BZ = 0  # Count slices for one BZ

    # Only iterate through path segments once (since they are the same for all BZs)
    for path_seg in pathSQE_params['user defined Qpoints']['path']:
        pt1_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][path_seg[0]])
        pt2_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][path_seg[1]])

        # Generate unique paths for this segment
        pt1_array, _ = generate_unique_paths(mtd_spacegroup, pt1_init, pt2_init, pathSQE_params['primtive to mantid T matrix'])
        print("\n{} {} symmetrically equivalent path segments".format(pt1_array.shape[0], path_seg))

        # Add the number of unique paths generated (which corresponds to slices)
        unique_slices_per_BZ += pt1_array.shape[0]

    # Multiply by the number of BZs to get the total slices
    num_total_slices = unique_slices_per_BZ * num_BZ

    return num_total_slices



def count_total_symPt_slices(pathSQE_params, BZ_list, mtd_spacegroup, num_datasets):
    """
    Calculates the total number of slices expected before processing starts.

    Parameters:
    - pathSQE_params: dictionary containing parameters, including 'user defined Qpoints'
    - BZ_list: list of BZ offsets
    - mtd_spacegroup: Mantid space group object used for symmetry generation
    - num_datasets: int, number of datasets being processed

    Returns:
    - total_slices: int, estimated total number of slices to be processed
    """

    num_BZ = len(BZ_list)  # Number of BZs
    unique_slices_per_BZ = 0  # Count unique slices in a single BZ

    # Iterate over all defined 1D symmetry points
    for symPoint in pathSQE_params['user defined Qpoints']['1d_points']:
        symPt_init = np.array(pathSQE_params['user defined Qpoints']['point_coords'][symPoint])

        # Find symmetry-equivalent points in a single BZ
        symPts_array = generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init)

        # Find unique symmetry points across all BZs
        allSymPts_array = generate_unique_symPts_inAllBZ(BZ_list, symPts_array)

        # Count total symmetry-equivalent slices
        unique_slices_per_BZ += allSymPts_array.shape[0]

    # Multiply by the number of datasets since each dataset requires slicing
    num_total_slices = unique_slices_per_BZ * num_datasets

    return num_total_slices




def should_print_update(num_processed_slices, num_total_slices, last_update_step, update_percent=5):
    """
    Determines whether to print a progress update based on percentage completion.

    Ensures that the first slice always triggers a print.

    Parameters:
    - num_processed_slices: int, number of slices processed so far.
    - num_total_slices: int, total number of slices.
    - last_update_step: int, last progress checkpoint (percentage completed).
    - update_percent: int, interval at which to print updates (default 5%).

    Returns:
    - should_print: bool, whether to print an update.
    - new_update_step: int, updated progress checkpoint.
    """
    progress = (num_processed_slices / num_total_slices) * 100
    
    # Always print after first slice
    if num_processed_slices == 1:
        return True, last_update_step

    # Print at each update_percent step
    if progress >= last_update_step + update_percent:
        return True, (last_update_step + update_percent)
    
    return False, last_update_step


# Convert to readable format (hh:mm:ss)
def format_time(seconds):
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{int(h):02d}:{int(m):02d}:{int(s):02d}"



def print_progress(num_processed_slices, num_total_slices, fold_timer, processing_type):
    """
    Prints the estimated time to completion based on the rolling average processing speed.

    Parameters:
    - num_processed_slices: int, number of slices processed so far
    - num_total_slices: int, total number of slices to be processed
    - fold_timer: float, timestamp when processing started (time.time())
    """

    elapsed_time = time.time() - fold_timer  # Time since start
    avg_time_per_slice = elapsed_time / max(num_processed_slices, 1)  # Avoid division by zero
    remaining_slices = num_total_slices - num_processed_slices
    estimated_remaining_time = remaining_slices * avg_time_per_slice

    print(f"\n\nProcessed: {num_processed_slices}/{num_total_slices} "
          f"({(num_processed_slices / num_total_slices) * 100:.2f}%) | "
          f"Elapsed: {format_time(elapsed_time)} | "
          f"{processing_type} estimated remaining time: {format_time(estimated_remaining_time)}\n\n")



def load_previous_simulation_progress(output_folder, user_defined_Qpoints, removed_BZ):
    """
    Loads and combines previously saved simulation data and norm files, but only for BZs that are fully processed.

    Parameters:
        output_folder (str): Path to the folder containing saved simulation results.
        user_defined_Qpoints (dict): Dictionary containing the path information.
        removed_BZ (np.ndarray): Nx3 array of BZ centers that were fully processed.

    Returns:
        tuple: (datas, norms) where:
            - datas: List of summed data arrays, one per path segment.
            - norms: List of summed norm arrays, one per path segment.
            - Returns ([], []) if no valid data is found.
    """
    datas = []
    norms = []

    # Convert removed_BZ to string format for easier matching in filenames
    removed_BZ_strs = {str(bz.tolist()) for bz in removed_BZ}
    print('BZ strings for match:', removed_BZ_strs)

    # Initialize empty lists for each segment (but return [] if nothing is found)
    num_segments = len(user_defined_Qpoints['path'])
    found_any_files = False  # Track if any valid files were found
    datas = [None] * num_segments
    norms = [None] * num_segments

    # Search for all BZ-related files in the output folder
    for filename in os.listdir(output_folder):
        match = re.search(r'BZ(\[.*?\])_foldedSim_data.npz', filename)  # Extract BZ array from filename
        if match:
            print('Match found:', filename)
            BZ_str = match.group(1)  # Extract the array portion inside brackets
            BZ_array = np.array([float(x) for x in BZ_str.strip("[]").split()])  # Convert to array

            # Check if this BZ is fully processed
            if str(BZ_array.tolist()) in removed_BZ_strs:
                data_path = os.path.join(output_folder, filename)
                norm_path = data_path.replace("_foldedSim_data.npz", "_foldedSim_norm.npz")

                print('Checking files:', data_path, norm_path)

                if os.path.exists(data_path) and os.path.exists(norm_path):
                    print('Files exist - loading data')
                    found_any_files = True  # At least one file was found

                    with np.load(data_path, allow_pickle=True) as data_npz:
                        data_arrays = [data_npz[key] for key in data_npz]  # Extract all arrays

                    with np.load(norm_path, allow_pickle=True) as norm_npz:
                        norm_arrays = [norm_npz[key] for key in norm_npz]  # Extract all arrays

                    # Add data to corresponding segment index
                    for seg_idx in range(num_segments):
                        if datas[seg_idx] is None:
                            datas[seg_idx] = data_arrays[seg_idx]
                            norms[seg_idx] = norm_arrays[seg_idx]
                        else:
                            datas[seg_idx] += data_arrays[seg_idx]
                            norms[seg_idx] += norm_arrays[seg_idx]

    # If no valid files were found, return empty lists
    if not found_any_files:
        print('No previous simulation data found.')
        return [], []

    print('Final data and norm lists populated')
    return datas, norms




##############################################################################################################################
############################################### FUNCTIONS FOR SLICE FILTERTING ###############################################
##############################################################################################################################



def evaluate_slice_quality(slice_name, filters, path_seg):
    # if not using filters, just say all slices are good
    if len(filters) == 0:
        return True
    
    filter_evals = []
    # filter based on fractional coverage of slice
    if 'coverage' in filters:
        slice_data = mtd[slice_name].getSignalArray()
        fractional_coverage = np.count_nonzero(~np.isnan(slice_data)) / slice_data.size
        print('fc: ', fractional_coverage)
        
        if fractional_coverage > 0.54:
            filter_evals.append(True)
        else:
            filter_evals.append(False)
    
    # UNDER DEV - filter looking at distribution of pixel intensities (contrast/sharpness) above elastic line (assumed bottom 3 rows)
    if 'contrast' in filters:
        slice_data = mtd[slice_name].getSignalArray()
        slice_data_withoutElasticLine = slice_data[:,0,0,3:]
        non_nan_values = slice_data_withoutElasticLine[~np.isnan(slice_data_withoutElasticLine)]
        #print(np.min(non_nan_values), np.max(non_nan_values), np.mean(non_nan_values), non_nan_values.shape)
        # Create histogram with log scale x-axis
        plt.figure()
        # Specify logarithmically spaced bin edges
        bins = np.logspace(-4.5, -1.5, 50)

        # Create histogram with log scale on both axes
        plt.hist(non_nan_values, bins=bins, log=True, edgecolor='black')
        plt.xlabel('Values (log scale)')
        plt.ylabel('Frequency (log scale)')
        plt.title('Histogram of Non-NaN Values (log-log scale)')
        plt.grid(True)
        plt.xscale('log')  # Set log scale on x-axis
        #plt.xlim([1e-4, np.max(non_nan_values)])

        # Save the plot as a PNG file
        plt.savefig(os.path.join('/SNS/ARCS/IPTS-5307/shared/Aiden/comprehensive/pathSQE_outputs/withFilters/',slice_name,'.png'))
        
        if 1 > 0.54:
            filter_evals.append(True)
        else:
            filter_evals.append(False)    


    # filter to get rid of M point artifacts in 40meV FeSi dataset (IPTS-5307)
    if 'M_bragg_tails' in filters: 
        if path_seg==('M','Gamma'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)
            
            # mask where bad signal is
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[:6,10:15] = 1  # E 5-8 and q 1/3 way from M to Gamma
            
            # mask data based on percentile
            E_ind_high = int(np.ceil(3/0.5))  
            slice = slice_data[:, E_ind_high:]
            mask = ~np.isnan(slice) & (slice != 0)
            segMax = np.percentile(slice[mask],98)
            perc_mask = slice_data >= segMax
                    
            if np.sum(feat_mask*perc_mask)>2:
                filter_evals.append(False)
            else:
                filter_evals.append(True)
        else:
            filter_evals.append(True)


    # filter to get rid of M point artifacts in 70meV FeSi dataset (IPTS-21211)
    if 'M_bragg_tails_70meV' in filters: 
        if path_seg==('M','Gamma'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)
            
            # mask where bad signal is
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[:10,10:17] = 1  # E 5.5-9 ish and q 0-0.2 on M to Gamma
            
            # mask data based on percentile
            E_ind_high = int(np.ceil(4/0.5))  
            slice = slice_data[:, E_ind_high:]
            mask = ~np.isnan(slice) & (slice != 0)
            segMax = np.percentile(slice[mask],95)
            perc_mask = slice_data >= segMax
                    
            if np.sum(feat_mask*perc_mask)>5:
                filter_evals.append(False)
            else:
                filter_evals.append(True)
        elif path_seg==('X','M'):
            slice_data = mtd[slice_name].getSignalArray()
            slice_data = np.squeeze(slice_data)
            
            # mask where bad signal is
            shape = slice_data.shape
            feat_mask = np.zeros(shape)
            feat_mask[12:,9:19] = 1  # E 5-10 ish and q 0.3-0.5 on X to M
            
            # mask data based on percentile
            E_ind_high = int(np.ceil(4/0.5))  
            slice = slice_data[:, E_ind_high:]
            mask = ~np.isnan(slice) & (slice != 0)
            segMax = np.percentile(slice[mask],95)
            perc_mask = slice_data >= segMax
                    
            if np.sum(feat_mask*perc_mask)>5:
                filter_evals.append(False)
            else:
                filter_evals.append(True)
        else:
            filter_evals.append(True)           

    
    # more filters...
    
    # final check to see if slice passes all filters
    slice_good = np.all(filter_evals)
    return slice_good



##############################################################################################################################
############################################# FUNCTIONS FOR PLOTTING FOLDED DATA #############################################
##############################################################################################################################



def set_seg_xlabels(seg_axis, i, dsl_fold):
    if i == len(dsl_fold)-1:
        prev_end = dsl_fold[i-1]['seg_end_name']
        curr_start = dsl_fold[i]['seg_start_name']
        curr_end = dsl_fold[i]['seg_end_name']
        
        pt_names = [prev_end, curr_start, curr_end]
        labels = [pt if len(pt) == 1 else '\\' + pt for pt in pt_names]
        
        if curr_start == prev_end:
            seg_axis.set_xticklabels([r'${}$'.format(labels[1]), r'${}$'.format(labels[2])])
        else:
            seg_axis.set_xticklabels([r'${},{}$'.format(labels[0], labels[1]), r'${}$'.format(labels[2])])
    
    elif i == 0:
        curr_start = dsl_fold[i]['seg_start_name']
        
        pt_names = [curr_start]
        labels = [pt if len(pt) == 1 else '\\' + pt for pt in pt_names]
        seg_axis.set_xticklabels([r'${}$'.format(labels[0]),''])
        
    else:
        prev_end = dsl_fold[i-1]['seg_end_name']
        curr_start = dsl_fold[i]['seg_start_name']
        
        pt_names = [prev_end, curr_start]
        labels = [pt if len(pt) == 1 else '\\' + pt for pt in pt_names]
        
        if curr_start == prev_end:
            seg_axis.set_xticklabels([r'${}$'.format(labels[1]),''])
        else:
            seg_axis.set_xticklabels([r'${},{}$'.format(labels[0], labels[1]),''])


def plot_along_path_foldedBZ(foldedSegNames, dsl_fold, Ei, vmi, vma, cma='jet'):
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams.update({'font.size': 28})
    
    num_segments = len(foldedSegNames)
    ratios = [dsl_fold[i]['inv_angstrom_ratio'] for i in range(num_segments)]

    
    # Add space for the colorbar
    fig, axes = plt.subplots(1, num_segments, gridspec_kw={'width_ratios': ratios}, figsize=(12, 5), subplot_kw={'projection':'mantid'}) #constrained_layout=True
    
    colormesh_pars = {
        'norm': SymLogNorm(linthresh=vmi, vmin=vmi, vmax=vma),
        'cmap': cma
    }
    
    for i in range(num_segments):
        x_start = dsl_fold[i]['qdim0_range'][0]
        x_end = dsl_fold[i]['qdim0_range'][2]
        
        if i == 0:
            ax1 = axes[i]
            ax1.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            #ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction='in')
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel('')
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction='in')
        else:
            ax = axes[i]
            ax.pcolormesh(mtd[foldedSegNames[i]], rasterized=True, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel('')
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels(ax, i, dsl_fold)
            ax.tick_params(direction='in')
    
    plt.subplots_adjust(wspace=0, hspace=0)

    # Create the colorbar
    cb_ax = fig.add_axes([.91,.124,.02,.754])
    cbar = fig.colorbar(mappable,orientation='vertical',cax=cb_ax)
    cbar.ax.tick_params(labelsize=15)

    plt.rcParams.update({'font.size': 10})

    return fig
    


def generate_BZ_Report_morePages(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold):
    pdf_filename = os.path.join(pathSQE_params['output directory'],'Reports/BZ_Report_{}.pdf'.format(BZ_offset))
    with PdfPages(pdf_filename) as pdf:
        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ, ax_path_foldedBZ = plt.subplots(subplot_kw={'projection': 'mantid'}, figsize=(8, 6))

        # plot the final folded path
        pathNames = [ds['Name'] for ds in dsl_fold]
        sortedWorkspaceNames = pathNames.copy()
        sortedWorkspaceNames = [name + '.nxs' for name in sortedWorkspaceNames]
    
        # Adaptive colorbars
        step = float(pathSQE_params['E bins'].split(',')[1])
        ind_aboveElasticLine = int(np.ceil(3 / step))
        pathMax = -1
        pathMin = 1e4

        for seg in sortedWorkspaceNames:              
            slice = mtd[seg].getSignalArray()[:, 0, 0, ind_aboveElasticLine:]
            
            # Mask NaN and zero values
            mask = ~np.isnan(slice) & (slice != 0)
            
            try:
                # Check if there are any valid values before calculating percentiles
                if np.any(mask):
                    segMin = np.percentile(slice[mask], 10)
                    segMax = np.percentile(slice[mask], 97)

                    # Update pathMin and pathMax
                    if segMax > pathMax:
                        pathMax = segMax
                    if segMin < pathMin:
                        pathMin = segMin
                else:
                    print(f"Warning: No valid data in segment {seg}, skipping percentile calculation.")
            
            except Exception as e:
                print(f"Error processing segment {seg}: {e}")

        
        fig_path_foldedBZ = plot_along_path_foldedBZ(foldedSegNames=sortedWorkspaceNames, dsl_fold=dsl_fold,
                                                    Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin,
                                                    vma=pathMax, cma=pathSQE_params['cmap'])

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ.suptitle(f'Folded data along path in BZ {BZ_offset}', fontsize=18)

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ)
        plt.close(fig_path_foldedBZ)

        # now make page for each path seg showing all of the individual slices that were combined
        for i in range(len(all_slice_names_in_BZ)):  # range(len(path['path'])):
            path_seg = [dsl_fold[i]['seg_start_name'], dsl_fold[i]['seg_end_name']]

            # Calculate the number of pages needed for this path segment
            num_pages = int(np.ceil(len(all_slice_names_in_BZ[i])/12))

            for page_num in range(num_pages):
                # Calculate the start and end index for slices on this page
                start_idx = page_num * 12
                end_idx = min((page_num + 1) * 12, len(all_slice_names_in_BZ[i]))

                # Create a new figure for each path_seg
                fig, axes = plt.subplots(nrows=3, ncols=4, subplot_kw={'projection': 'mantid'}, figsize=(12, 9))
                plt.rcParams['figure.dpi'] = 300
                plt.rcParams['savefig.dpi'] = 300
                plt.set_cmap(pathSQE_params['cmap'])

                for j, slice_info in enumerate(all_slice_names_in_BZ[i][start_idx:end_idx]):
                    slice_name = slice_info[0]
                    slice_ID = slice_info[1]
                    
                    # Adjust the subplot position based on the number of rows and columns
                    row_index = j // 4
                    col_index = j % 4
                        
                    step = float(pathSQE_params['E bins'].split(',')[1])
                    ind_aboveElasticLine = int(np.ceil(3/step))
                    slice = mtd[slice_name].getSignalArray()[:,0,0,ind_aboveElasticLine:]
                    # Mask NaN and zero values
                    mask = ~np.isnan(slice) & (slice != 0)
                    
                    if np.any(mask != 0):
                        segMin = np.percentile(slice[mask],10)
                        segMax = np.percentile(slice[mask],97)
                        axes[row_index, col_index].pcolormesh(mtd[slice_name], vmin=segMin, vmax=segMax, rasterized=True)
                        axes[row_index, col_index].set_title('ID {}'.format(slice_ID))
                        if not all_slice_evals[i][j]:
                            axes[row_index, col_index].set_title('Flagged Slice (ID {})'.format(slice_ID), color='red')

                        # x ticks and labels
                        point1_str, point2_str = slice_name.replace("pathSQE_", "").split("_to_")
                        point1 = np.array([float(x) for x in point1_str.split(",")])
                        point2 = np.array([float(x) for x in point2_str.split(",")])
                        #print(slice_name, point1_str, point2_str, point1, point2)
                        x_min_exp, x_max_exp = axes[row_index, col_index].get_xlim()
                        #print(x_min_exp, x_max_exp)
                        axes[row_index, col_index].set_xticks([x_min_exp, x_max_exp])
                        axes[row_index, col_index].set_xticklabels([f'{tick}' for tick in [point1, point2]])

                    else:
                        axes[row_index, col_index].text(0.5, 0.5, 'No data along this segment', ha='center',
                                                       va='center', fontsize=10, color='red')

                # Add BZ_offset and path_seg information at the top of each page
                fig.suptitle(f'BZ: {BZ_offset}, Path Segment: {path_seg[0]} to {path_seg[1]} - Page {page_num + 1}',
                             y=0.98, fontsize=14)

                # Adjust layout to prevent overlapping labels
                fig.tight_layout(h_pad=0.3, w_pad=0.3)

                # Save the current figure to the PDF
                pdf.savefig(fig)
                plt.close('all')



def IqE_to_SED(SED_ws,T,Ebins):
    SED_ws_shape = SED_ws.getSignalArray().shape
    k = 8.617e-2
    # Split the string to extract min, step, and max values
    min_val, step_val, max_val = map(float, Ebins.split(','))

    # Generate the vector using numpy's arange function
    E_vals = np.arange(min_val + step_val/2, max_val + step_val/2, step_val)
    E_terms = np.exp(E_vals/(k*T))
    factors = E_vals / (1 + (1/(E_terms-1)))
    reshaped_factors = np.tile(factors.reshape(1, 1, 1, len(factors)), (SED_ws_shape[0], SED_ws_shape[1], SED_ws_shape[2], 1))
    SED_ws.setSignalArray(SED_ws.getSignalArray() * reshaped_factors)
        
    return SED_ws


def plot_SED_along_path_foldedBZ(foldedSegNames,dsl_fold,pathSQE_params):
    T= pathSQE_params['T and Ei conditions'][0][0]
    Ebins=pathSQE_params['E bins']
    Ei= pathSQE_params['T and Ei conditions'][0][1]
    cma=pathSQE_params['cmap']
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams.update({'font.size': 28})

    fig=plt.figure(figsize=(12,5))
    colormesh_pars={}
    colormesh_pars['cmap']=cma
    vmi=0
    vma=0
        
    num_segments = len(foldedSegNames)
    ratios = []
    for i in range(num_segments):
        ratios.append(dsl_fold[i]['inv_angstrom_ratio'])
    
    gs = GridSpec(1, num_segments,width_ratios=ratios,wspace=0)

    for i in range(num_segments):
        x_start = dsl_fold[i]['qdim0_range'][0]
        x_end = dsl_fold[i]['qdim0_range'][2]
        
        SED_ws = mtd[foldedSegNames[i]].clone()
        SED_ws = IqE_to_SED(SED_ws=SED_ws,T=T,Ebins=Ebins)
        SaveMD(SED_ws, Filename=os.path.join(pathSQE_params['output directory'],'bzfold_nxs_files','{}2{}_folded_SED.nxs'.format(dsl_fold[i]['seg_start_name'],dsl_fold[i]['seg_end_name'])))
        
        if i == 0:
            slice = SED_ws.getSignalArray()
            mask = ~np.isnan(slice) & (slice != 0)
            vmi = np.percentile(slice[mask],10)
            vma = np.percentile(slice[mask],97)
            colormesh_pars['norm']=SymLogNorm(linthresh=vmi,vmin=vmi,vmax=vma)#(linthresh=9e-6,vmin=9e-6,vmax=2e-3)
                        
            ax1 = plt.subplot(gs[i],projection='mantid')
            ax1.pcolormesh(SED_ws, **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            #ax1.set_ylim(0.,Ei)
            ax1.set_xticks([x_start, x_end])
            ax1.set_xticklabels(['{}'.format(dsl_fold[i]['seg_start_name']),''])
            ax1.tick_params(direction='in')
        elif i == num_segments-1:
            axL = plt.subplot(gs[i],sharey=ax1,projection='mantid')
            cb = axL.pcolormesh(SED_ws, **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel('')
            axL.set_xticks([x_start, x_end])
            axL.set_xticklabels(['${}$'.format(dsl_fold[i]['seg_start_name']), '${}$'.format(dsl_fold[i]['seg_end_name'])])
            axL.tick_params(direction='in')
        else:
            ax = plt.subplot(gs[i],sharey=ax1,projection='mantid')
            ax.pcolormesh(SED_ws, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel('')
            ax.set_xticks([x_start, x_end])
            ax.set_xticklabels(['{}'.format(dsl_fold[i]['seg_start_name']),''])
            ax.tick_params(direction='in')

    fig.colorbar(cb)
    plt.rcParams.update({'font.size': 10})

    return fig



########################################## ADDED FOR 1D SYM POINTS ##########################################



def generate_unique_symPts(mtd_spacegroup, symPoint, symPt_init):
    spacegroup = SpaceGroupFactory.createSpaceGroup(mtd_spacegroup)
    pg = spacegroup.getPointGroup()
    rotations = [np.array([so.transformHKL((1, 0, 0)), so.transformHKL((0, 1, 0)), so.transformHKL((0, 0, 1))])
                 for so in pg.getSymmetryOperations()]

    new_symPts = np.dot(rotations, symPt_init)

    # Find unique rows
    unique_symPts_array = np.unique(new_symPts, axis=0)

    return unique_symPts_array



def generate_unique_symPts_inAllBZ(BZ_list, symPts_array):
    unique_points_set = set()
    
    # Iterate through each BZ center, which is a NumPy array
    for bz_center in BZ_list:
        # Translate symPts_array for the current BZ center
        translated_points = symPts_array + bz_center
        
        # Add each translated point to the set (as a tuple to ensure hashability)
        for point in translated_points:
            unique_points_set.add(tuple(point))
    
    # Convert the set back to an Nx3 numpy array
    allSymPts_array = np.array(list(unique_points_set))
    
    return allSymPts_array



def make_slice_desc_1DSymPoints(pathSQE_params, symPoint, point):
    # slice description with desired dims and bins and unique name
    slice_desc={'QDimension0':'1,0,0',
                'QDimension1':'0,1,0',
                'QDimension2':'0,0,1',
                'Dimension0Name':'QDimension0',
                'Dimension0Binning':conv_to_desc_string([point[0]-pathSQE_params['qdim0 step size'], point[0]+pathSQE_params['qdim0 step size']]),
                'Dimension1Name':'QDimension1',
                'Dimension1Binning':conv_to_desc_string([point[1]-pathSQE_params['qdim1 integration range'], point[1]+pathSQE_params['qdim1 integration range']]),
                'Dimension2Name':'QDimension2',
                'Dimension2Binning':conv_to_desc_string([point[2]-pathSQE_params['qdim2 integration range'], point[2]+pathSQE_params['qdim2 integration range']]),
                'Dimension3Name':'DeltaE',
                'Dimension3Binning':pathSQE_params['E bins'],
                'Name':'pathSQE_{}point_{}'.format(symPoint, point)}

    #slice_desc['Name'] = slice_desc['Name'].replace(',', '_').replace('.', '').replace('-','m')
    
    return slice_desc



def generate_BZ_Report_1DSymPoints(pathSQE_params, symPoint, symPt_slice_arrays, symPt_slice_names_allTemp, u, v):
    pdf_filename = os.path.join(pathSQE_params['output directory'],'Reports/{}point_Report.pdf'.format(symPoint))
    with PdfPages(pdf_filename) as pdf:
        fig_scatter, ax_scatter = plt.subplots(figsize=(8, 6))

        # Convert list to NumPy array for sorting
        symPt_slice_array = np.array(symPt_slice_arrays)
        
        # Compute out-of-plane unit normal
        u = u.astype(float) / np.linalg.norm(u)
        v = v.astype(float) / np.linalg.norm(v)
        normal_vec = np.cross(u, v)
        normal_vec /= np.linalg.norm(normal_vec)  # Normalize

        # Transform (H, K, L) points into (u, v) plane coordinates
        x = np.dot(symPt_slice_array[:, :3], u)  # Project onto u
        y = np.dot(symPt_slice_array[:, :3], v)  # Project onto v
        z = np.dot(symPt_slice_array[:, :3], normal_vec)  # Out-of-plane coordinate

        # Sort by out-of-plane coordinate
        sorted_indices = np.argsort(z)
        x, y, z = x[sorted_indices], y[sorted_indices], z[sorted_indices]

        # Normalize z values between 0 and 1
        norm = Normalize(vmin=np.min(z), vmax=np.max(z))
        norm_z = norm(z)

        # Define marker size based on z values (larger for closer points)
        marker_size = 500 * (1 - norm_z) ** 2 + 30  

        # Define colors based on z values
        colors = viridis(norm_z)

        # Plot scatter plot
        scatter = ax_scatter.scatter(x, y, c=colors, s=marker_size, alpha=1)

        # Add color bar
        sm = ScalarMappable(norm=norm, cmap=viridis)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax_scatter, label='Out-of-Plane Coord')

        # Set grid lines aligned with u and v
        x_half_integer_locs = np.arange(np.floor(np.min(x))-0.5, np.ceil(np.max(x)) + 1.5, 1)
        y_half_integer_locs = np.arange(np.floor(np.min(y))-0.5, np.ceil(np.max(y)) + 1.5, 1)
        ax_scatter.xaxis.set_major_locator(plt.FixedLocator(x_half_integer_locs))
        ax_scatter.yaxis.set_major_locator(plt.FixedLocator(y_half_integer_locs))
        ax_scatter.xaxis.set_minor_locator(plt.MultipleLocator(1))
        ax_scatter.yaxis.set_minor_locator(plt.MultipleLocator(1))
        ax_scatter.xaxis.set_minor_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax_scatter.yaxis.set_minor_formatter(plt.FuncFormatter(lambda y, _: f"{int(y)}"))
        ax_scatter.grid(which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.3)

        # Set axis labels based on u and v vectors
        ax_scatter.set_xlabel(f"u = {u}")
        ax_scatter.set_ylabel(f"v = {v}")
        fig_scatter.suptitle('{} {} points with data'.format(len(symPt_slice_arrays), symPoint), y=0.95, fontsize=18)

        # Save figure to PDF
        pdf.savefig(fig_scatter)
        plt.close(fig_scatter)
        

        # Get e range vector for plotting
        spectrum_labels = pathSQE_params['T and Ei conditions']
        start, step, end = map(float, pathSQE_params['E bins'].split(','))
        evec_tot = np.arange(start + 0.5 * step, end - 0.5 * step + step, step)
        
        # Calculate number of pages needed
        num_pages = int(np.ceil(len(symPt_slice_arrays) / 12))

        for page_num in range(num_pages):
            start_idx = page_num * 12
            end_idx = min((page_num + 1) * 12, len(symPt_slice_arrays))
            symPts_subset = symPt_slice_arrays[start_idx:end_idx]

            # Create figure for each page
            fig, axes = plt.subplots(nrows=3, ncols=4, subplot_kw={'projection': 'mantid'}, figsize=(12, 9))
            plt.rcParams['figure.dpi'] = 300
            plt.rcParams['savefig.dpi'] = 300
            
            for j in range(len(symPts_subset)):
                row_index = j // 4
                col_index = j % 4

                # Get colormap
                cmap = plt.get_cmap('coolwarm')

                for k in range(len(symPt_slice_names_allTemp)):
                    slice_names_subset = symPt_slice_names_allTemp[k][start_idx:end_idx]
                    slice_name = slice_names_subset[j]

                    color = None
                    if len(symPt_slice_names_allTemp)>1:
                        norm_k = k / (len(symPt_slice_names_allTemp) - 1)
                        color = cmap(norm_k)
                    else:
                        color='tab:green'

                    if k == 0:
                        # First temperature: Normalize using percentile
                        signal_tot = mtd[slice_name].getSignalArray()[0, 0, 0, :]
                        nonNaN_indices = np.where(~np.isnan(signal_tot))[0]
                        signal = signal_tot[nonNaN_indices]
                        
                        vert_offset = 1.1 * np.percentile(signal, 98)  # Use 95th percentile instead of max
                        peak_scaling = np.percentile(signal, 98)  # Scale sim by 95th percentile

                        signal += vert_offset
                        error = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
                        error = error[nonNaN_indices]
                        e_vec = evec_tot[nonNaN_indices]

                        axes[row_index, col_index].errorbar(x=e_vec, y=signal, yerr=error, color=color)
                        axes[row_index, col_index].axhline(y=vert_offset, color='black', alpha=0.15)

                        # **Only perform simulation if pathSQE_params['perform simulations'] is True**
                        if pathSQE_params['perform simulations']:
                            #print('T? ', pathSQE_params['T and Ei conditions'][k][0])
                            #print(symPts_subset[j])
                            sim_SQE_output = sim_SQE(pathSQE_params, [np.linalg.inv(pathSQE_params['primtive to mantid T matrix']) @ np.array(symPts_subset[j])], Temperature=pathSQE_params['T and Ei conditions'][k][0])
                            unique_freq, summed_intensities, freq_range, spectrum = SQE_to_1d_spectrum(pathSQE_params, sim_SQE_output)

                            sim_mask = (unique_freq > np.min(e_vec)) & (unique_freq < np.max(e_vec))
                            unique_freq = unique_freq[sim_mask]
                            summed_intensities = summed_intensities[sim_mask]

                            sim_mask = (freq_range > np.min(e_vec)) & (freq_range < np.max(e_vec))
                            freq_range = freq_range[sim_mask]
                            spectrum = spectrum[sim_mask]

                            summed_intensities *= peak_scaling / np.max(spectrum)
                            spectrum *= peak_scaling / np.max(spectrum)

                            axes[row_index, col_index].scatter(unique_freq, summed_intensities, color='black', s=5)
                            axes[row_index, col_index].plot(freq_range, spectrum, color='black')

                            if symPoint == 'R':
                                np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name,'_simSpectrum.npy'), np.vstack((freq_range, spectrum)))
                                np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name,'_simScatter.npy'), np.vstack((unique_freq, summed_intensities)))

                    else:
                        # Other temperatures
                        signal_tot2 = mtd[slice_name].getSignalArray()[0, 0, 0, :]
                        nonNaN_indices2 = np.where(~np.isnan(signal_tot2))[0]
                        signal2 = signal_tot2[nonNaN_indices2]
                        signal2 += (k + 1) * vert_offset

                        error2 = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
                        error2 = error2[nonNaN_indices2]
                        e_vec2 = evec_tot[nonNaN_indices2]

                        axes[row_index, col_index].errorbar(x=e_vec2, y=signal2, yerr=error2, color=color)
                        axes[row_index, col_index].axhline(y=(k + 1) * vert_offset, color='black', alpha=0.15)

                        if k == len(symPt_slice_names_allTemp) - 1:
                            axes[row_index, col_index].set_ylim([0, (k+3)*vert_offset])

                    axes[row_index, col_index].set_xlim([start, end])
                    axes[row_index, col_index].set_xlabel('Energy (meV)')
                    axes[row_index, col_index].set_title('{}'.format(symPts_subset[j]))

                # **Add a legend to the top-left subplot**
                if j == 0 and len(symPt_slice_names_allTemp)>1:  # Top-left plot on each page
                    legend_handles = [plt.Line2D([0], [0], color=cmap(i / (len(spectrum_labels) - 1)), lw=3) 
                                      for i in range(len(spectrum_labels))]
                    axes[row_index, col_index].legend(legend_handles, spectrum_labels, loc='upper right', fontsize=8)

            # Adjust layout and save the PDF page
            fig.tight_layout(h_pad=0.3, w_pad=0.3)
            pdf.savefig(fig)
            plt.close('all')



def plot_simple_1d_pt(pathSQE_params, symPoint, symPt_init, slice_names):
    # Get e range vector for plotting
    spectrum_labels = pathSQE_params['T and Ei conditions']
    start, step, end = map(float, pathSQE_params['E bins'].split(','))
    evec_tot = np.arange(start + 0.5 * step, end - 0.5 * step + step, step)

    # Create figure for each page
    fig, axes = plt.subplots(subplot_kw={'projection': 'mantid'})
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

    # Get colormap
    cmap = plt.get_cmap('coolwarm')

    for k in range(len(pathSQE_params['T and Ei conditions'])):
        slice_name = slice_names[k]

        color = None
        if len(pathSQE_params['T and Ei conditions'])>1:
            norm_k = k / (len(pathSQE_params['T and Ei conditions']) - 1)
            color = cmap(norm_k)
        else:
            color='tab:green'

        if k == 0:
            # First temperature: Normalize using percentile
            signal_tot = mtd[slice_name].getSignalArray()[0, 0, 0, :]
            nonNaN_indices = np.where(~np.isnan(signal_tot))[0]
            signal = signal_tot[nonNaN_indices]
            
            vert_offset = 1.1 * np.percentile(signal, 98)  # Use 95th percentile instead of max
            peak_scaling = np.percentile(signal, 98)  # Scale sim by 95th percentile

            signal += vert_offset
            error = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
            error = error[nonNaN_indices]
            e_vec = evec_tot[nonNaN_indices]

            axes.errorbar(x=e_vec, y=signal, yerr=error, color=color)
            axes.axhline(y=vert_offset, color='black', alpha=0.15)

            # **Only perform simulation if pathSQE_params['perform simulations'] is True**
            if pathSQE_params['perform simulations']:
                #print('T? ', pathSQE_params['T and Ei conditions'][k][0])
                #print(symPts_subset[j])
                sim_SQE_output = sim_SQE(pathSQE_params, [np.linalg.inv(pathSQE_params['primtive to mantid T matrix']) @ symPt_init], Temperature=pathSQE_params['T and Ei conditions'][k][0])
                unique_freq, summed_intensities, freq_range, spectrum = SQE_to_1d_spectrum(pathSQE_params, sim_SQE_output)

                sim_mask = (unique_freq > np.min(e_vec)) & (unique_freq < np.max(e_vec))
                unique_freq = unique_freq[sim_mask]
                summed_intensities = summed_intensities[sim_mask]

                sim_mask = (freq_range > np.min(e_vec)) & (freq_range < np.max(e_vec))
                freq_range = freq_range[sim_mask]
                spectrum = spectrum[sim_mask]

                summed_intensities *= peak_scaling / np.max(spectrum)
                spectrum *= peak_scaling / np.max(spectrum)

                axes.scatter(unique_freq, summed_intensities, color='black', s=5)
                axes.plot(freq_range, spectrum, color='black')

                np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name+'_simSpectrum.npy'), np.vstack((freq_range, spectrum)))
                np.save(os.path.join(pathSQE_params['output directory'],'nxs_files/',slice_name+'_simScatter.npy'), np.vstack((unique_freq, summed_intensities)))

        else:
            # Other temperatures
            signal_tot2 = mtd[slice_name].getSignalArray()[0, 0, 0, :]
            nonNaN_indices2 = np.where(~np.isnan(signal_tot2))[0]
            signal2 = signal_tot2[nonNaN_indices2]
            signal2 += (k + 1) * vert_offset

            error2 = np.sqrt(mtd[slice_name].getErrorSquaredArray()[0, 0, 0, :])
            error2 = error2[nonNaN_indices2]
            e_vec2 = evec_tot[nonNaN_indices2]

            axes.errorbar(x=e_vec2, y=signal2, yerr=error2, color=color)
            axes.axhline(y=(k + 1) * vert_offset, color='black', alpha=0.15)

            if k == len(pathSQE_params['T and Ei conditions']) - 1:
                axes.set_ylim([0, (k+3)*vert_offset])

        axes.set_xlim([start, end])
        axes.set_xlabel('Energy (meV)')
        axes.set_title('{}'.format(symPt_init))

    # **Add a legend if applicable**
    if len(pathSQE_params['T and Ei conditions'])>1: 
        legend_handles = [plt.Line2D([0], [0], color=cmap(i / (len(spectrum_labels) - 1)), lw=3) 
                            for i in range(len(spectrum_labels))]
        axes.legend(legend_handles, spectrum_labels, loc='upper right', fontsize=8)

    # Adjust layout and save the PDF page
    fig.savefig(os.path.join(pathSQE_params['output directory'],'{}point_{}.png'.format(symPoint, symPt_init)))
    plt.close('all')







#################################### SIMULATION STUFF ####################################




def determine_forces_file():
    """Determine whether FORCE_CONSTANTS or FORCE_SETS is present in the directory."""
    if os.path.exists("FORCE_CONSTANTS"):
        forceconstants_file = "FORCE_CONSTANTS"
        forcesets_file = None
    elif os.path.exists("FORCE_SETS"):
        forceconstants_file = None
        forcesets_file = "FORCE_SETS"
    else:
        raise FileNotFoundError("Error: Neither FORCE_CONSTANTS nor FORCE_SETS found in the current directory.")

    return forceconstants_file, forcesets_file



def get_coherent_scattering_lengths(sample_formula, filename="ScattLengths.txt"):
    """
    Given a sample's chemical formula, find the corresponding coherent scattering lengths.

    Parameters:
        sample_formula (str): Sample composition as a space-separated string (e.g., "Fe Si").
        filename (str): Name of the scattering length data file.

    Returns:
        dict: Coherent scattering lengths for each element in the sample.
    """
    # Check if the scattering lengths file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Error: {filename} not found.")

    # Read the scattering lengths file and store values in a dictionary
    coh_scatter_length = {}

    with open(filename, "r") as file:
        lines = file.readlines()

        for line in lines:
            # Skip headers or malformed lines
            if line.startswith("atom_name") or line.strip() == "":
                continue
            
            parts = line.split()
            if len(parts) < 3:
                continue  # Ensure enough columns exist
            
            element = parts[0]  # First column is the element name
            coh_b = parts[1]  # Second column is the coherent scattering length

            # Some entries may have complex values; extract real part if needed
            try:
                coh_b = float(coh_b.split("-")[0])  # Extract real part (ignoring imaginary)
            except ValueError:
                continue  # Skip if conversion fails (e.g., for weird formats)

            coh_scatter_length[element] = coh_b

    # Extract scattering lengths for the sample elements
    sample_elements = sample_formula.split()
    sample_scattering = {elem: coh_scatter_length[elem] for elem in sample_elements if elem in coh_scatter_length}

    return sample_scattering


''' WORK FOR FeSi BUT NOT Ge
def construct_sim_Qpts(pt1, BZ_offset, q_diff, step_size):
    """
    Constructs simulated Q points at finer resolution before binning.

    Parameters:
    - pt1: ndarray, starting point of the Q trajectory
    - BZ_offset: ndarray, offset to apply to the start point
    - q_diff: ndarray, direction vector along which Q points are generated
    - step_size: float, experimental step size along the q_diff direction

    Returns:
    - Qpoints: list of ndarrays, generated Q points with finer step size
    """
    from scipy.linalg import norm

    q_start = pt1 + BZ_offset
    finer_step = step_size / 2  # Finer grid for simulation

    # Generate finer Q-points
    Qpoints = [q_start + i * q_diff for i in np.arange(finer_step / 2, 0.5, finer_step)]

    # Avoid DW singularity at Q = [0,0,0]
    Qpoints_abs = norm(Qpoints, axis=1)
    zero_index = np.where(Qpoints_abs == 0)[0]
    if zero_index.shape[0] != 0:
        Qpoints[zero_index[0]] = np.array([1e-6, 1e-6, 1e-6])

    return Qpoints
'''


def construct_sim_Qpts(pt1, pt2, q_diff, step_size, prim2mantid, BZ_offset=np.array([0,0,0])):
    """
    Constructs simulated Q points along the direction of q_diff.

    Parameters:
    - pt1: ndarray, starting point of the Q trajectory
    - pt2: ndarray, ending point of the Q trajectory
    - BZ_offset: ndarray, offset to apply to both start and end points
    - q_diff: ndarray, direction vector along which Q points are generated
    - step_size: float, step size along the q_diff direction

    Returns:
    - Qpoints: list of ndarrays, generated Q points
    """
    from scipy.linalg import norm

    q_start = pt1 + BZ_offset
    q_end = pt2 + BZ_offset  # Explicitly set the correct endpoint

    # Compute the projection of q_end onto q_start + q_diff * t
    total_distance = np.max(np.abs(q_end - q_start))
    step_size = step_size / 4 # for finer sim                            # EDIT USED TO BE /2
    num_steps = int(np.round(total_distance / step_size))

    # Generate Q points
    Qpoints = [q_start + i * step_size * q_diff for i in np.arange(0.5, num_steps)]
    #print('len Q pts ', len(Qpoints), Qpoints[0], Qpoints[-1])

    # Avoid DW singularity at Q = [0,0,0]
    Qpoints = np.array(Qpoints)
    #print('Q pt array shape ', Qpoints.shape)
    zero_index = np.where(norm(Qpoints, axis=1) == 0)[0]
    if zero_index.size > 0:
        Qpoints[zero_index[0]] = np.array([1e-6, 1e-6, 1e-6])

    Qpoints = Qpoints @ np.linalg.inv(prim2mantid)

    return Qpoints




def run(phonon, Qpoints, temperature, atomic_form_factor_func=None, scattering_lengths=None):
    from phonopy import load

    # Transformation to the Q-points in reciprocal primitive basis vectors
    Q_prim = np.dot(Qpoints, phonon.primitive_matrix)
    # Q_prim must be passed to the phonopy dynamical structure factor code.
    phonon.run_dynamic_structure_factor(
        Q_prim,
        temperature,
        atomic_form_factor_func=atomic_form_factor_func,
        scattering_lengths=scattering_lengths,
        freq_min=8e-2)
    dsf = phonon.dynamic_structure_factor
    q_cartesian = np.dot(dsf.qpoints,
                         np.linalg.inv(phonon.primitive.get_cell()).T)
    distances = np.sqrt((q_cartesian ** 2).sum(axis=1))

    SandE = np.array([dsf.frequencies,dsf.dynamic_structure_factors])
    return SandE



def sim_SQE(pathSQE_params, Qpoint, Temperature):
    from phonopy import load

    #### START OF USER INPUTS #####
    #print('Q pt ', Qpoint)   # needs to be a list of arrays so like [array([-1.5, -1.5,  0.5])] when printed
    ## Q inputs ##
    primitive_cell = [[1,0,0],[0,1,0],[0,0,1]]
    supercell = pathSQE_params['supercell dimensions']
    forceconstants_file, forcesets_file = determine_forces_file()

    # Neutron coherent scattering length can be found at https://www.ncnr.nist.gov/resources/n-lengths/
    coh_scatter_length = get_coherent_scattering_lengths(pathSQE_params['sample'])
    THz2meV = 4.1357

    #### END OF USER INPUTS ####
    phonon = load(supercell_matrix=supercell,
                  primitive_matrix=primitive_cell,
                  unitcell_filename="POSCAR",
                  force_sets_filename=forcesets_file,
                  force_constants_filename=forceconstants_file
                  )

    #print(Qpoints)
    # Mesh sampling phonon calculation is needed for Debye-Waller factor.
    # This must be done with is_mesh_symmetry=False and with_eigenvectors=True.
    mesh = [11, 11, 11]
    phonon.run_mesh(mesh,
                    is_mesh_symmetry=False, # symmetry must be off
                    with_eigenvectors=True) # eigenvectors must be true
    temperature = Temperature

    # For INS, scattering length has to be given.
    # The following values is obtained at (Coh b)
    # https://www.nist.gov/ncnr/neutron-scattering-lengths-list
    output = run(phonon,
                 Qpoint,
                 temperature,
                 scattering_lengths=coh_scatter_length)
    
    ## output has shape as (2,len(Qpoints),branches), the [0,:,:] is for frequency and the [1,:,:] is for SQE; The frequency is in THz unit
    for i in range(len(Qpoint)):
        output[0,i,:] *= THz2meV

    return output



def SQE_to_1d_spectrum(pathSQE_params, output):
    frequencies = output[0, 0, :]
    intensities = output[1, 0, :]

    # define spectral binning info
    E_min = float(pathSQE_params['E bins'].split(',')[0])
    E_max = float(pathSQE_params['E bins'].split(',')[2])
    
    e_resolution = 1.177 * pathSQE_params['energy blurring sigma'] * float(pathSQE_params['E bins'].split(',')[1]) # gauss sigma to lorentz hwhm

    # Tolerance for considering frequencies equal (0.1%)
    tolerance = 0.001

    # Initialize arrays to store unique frequencies and their corresponding summed intensities
    unique_freq = []
    summed_intensities = []

    # Iterate over the frequency array
    for freq in np.unique(frequencies):
        # Exclude frequencies less than 0.1 to avoid singularity at gamma
        if freq < 0.1:
            continue
        
        # Check if the current frequency is sufficiently separated from existing unique frequencies
        if all(abs(freq - existing_freq) > tolerance * existing_freq for existing_freq in unique_freq):
            # Find indices where the current frequency occurs in the original array
            indices = np.where(np.abs(frequencies - freq) / freq <= tolerance)[0]
            # Sum corresponding intensities for the current frequency
            total_intensity = np.sum(intensities[indices])
            # Append the unique frequency and its corresponding summed intensity
            unique_freq.append(freq)
            summed_intensities.append(total_intensity)


    # Convert lists to numpy arrays
    unique_freq = np.array(unique_freq)
    summed_intensities = np.array(summed_intensities)
    
    # Generate frequency range
    freq_range = np.arange(E_min, E_max, 0.02)
    
    # Initialize spectrum array
    spectrum = np.zeros(freq_range.shape)
    ind_spec = []

    # Iterate over each peak
    for freq, intensity in zip(unique_freq, summed_intensities):
        # Add the Gaussian function to the spectrum
        #gaussian_peak = intensity * np.exp(-((freq_range - freq) / (2 * e_resolution)) ** 2)
        lorentzian_peak = intensity / (1 + ((freq_range - freq) / e_resolution) ** 2)
        ind_spec.append(lorentzian_peak)
        spectrum += lorentzian_peak


    return unique_freq, summed_intensities, freq_range, spectrum



def SQE_to_2d_spectrum(pathSQE_params, output): 
    """
    Bins simulated SQE data to a finer grid, applies smoothing, 
    and then rebins to match experimental resolution.

    Parameters:
    - pathSQE_params: dict, contains experimental binning information
    - output: ndarray, simulated data (fine resolution)

    Returns:
    - BinnedSQE: ndarray, rebinned SQE to match experimental binning
    """
    from scipy.ndimage import gaussian_filter

    E_min = float(pathSQE_params['E bins'].split(',')[0])
    E_max = float(pathSQE_params['E bins'].split(',')[2])
    E_step = float(pathSQE_params['E bins'].split(',')[1])
    
    nql_fine = output.shape[1]  # Fine Q resolution from simulation
    finer_E_step = E_step / 4   # Simulated data has finer E step                           # THIS AND BELOW 2
    nql = nql_fine // 4         # Experimental Q resolution (factor of 2 binning)
    
    # Generate energy bin edges
    E_bin_edges = np.arange(E_min, E_max + E_step, E_step)
    ne_exp = len(E_bin_edges) - 1

    # Fine binning storage
    FineBinnedSQE = np.zeros((nql_fine, len(np.arange(E_min, E_max, finer_E_step))))
    
    # Step 1: Bin energy values finely
    for ih in range(nql_fine):  # qpoints
        for j in range(len(output[0, 0, :])):  # branches
            energy_val = output[0][ih][j]
            hist, _ = np.histogram([energy_val], bins=np.arange(E_min, E_max + finer_E_step, finer_E_step), weights=[output[1][ih][j]])
            FineBinnedSQE[ih, :] += hist

    # Step 2: Apply Gaussian smoothing
    FineBinnedSQE = gaussian_filter(FineBinnedSQE, sigma=(0.5, pathSQE_params['energy blurring sigma']*4))                               # THIS ONE

    # Step 3: Rebin Q dimension (fine → coarse) using NumPy reshape and sum
    CoarseBinnedSQE = FineBinnedSQE.reshape(nql, 4, -1).sum(axis=1)  # Sum pairs of adjacent fine Q bins             # THIS ALSO 2

    # Step 4: Rebin E dimension (fine → coarse) using NumPy histogram
    BinnedSQE = np.zeros((nql, ne_exp))
    for ih in range(nql):
        hist, _ = np.histogram(np.linspace(E_min, E_max, FineBinnedSQE.shape[1]), bins=E_bin_edges, weights=CoarseBinnedSQE[ih, :])
        BinnedSQE[ih, :] = hist

    return BinnedSQE



def fold_symmetric_simulations(all_sims_in_BZ):
    """
    Folds all symmetrically equivalent simulated data by summing and normalizing.

    Parameters:
    - all_sims_in_BZ: List of lists of binned SQE arrays.
                      Outer list indexes path segments, inner lists contain equivalent simulations.

    Returns:
    - folded_results: List of tuples, each containing:
        (summed_data, summed_norm, final_folded_data) for each segment.
    """
    folded_results = []

    for segment_sims in all_sims_in_BZ:
        # Ensure the segment contains at least one simulation
        if not segment_sims:
            continue  
        
        # Initialize summed arrays with zeros
        summed_data = np.zeros_like(segment_sims[0])
        summed_norm = np.zeros_like(segment_sims[0])

        for sim in segment_sims:
            norm = np.where(~np.isnan(sim), 1, 0)  # Norm is 1 where data is valid, 0 otherwise
            sim_filled = np.nan_to_num(sim)  # Replace NaNs with zero for summation
            
            summed_data += sim_filled
            summed_norm += norm
        
        # Avoid division by zero by setting invalid norms to NaN
        summed_norm[summed_norm == 0] = np.nan
        
        # Compute final folded data
        final_folded_data = summed_data / summed_norm

        folded_results.append((np.nan_to_num(summed_data, nan=0.0), np.nan_to_num(summed_norm, nan=0.0), final_folded_data))

    return folded_results




def plot_along_path_foldedBZ_sim(sim_folded_data, dsl_fold, Ei, vmi, vma, cma='jet'):
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams.update({'font.size': 28})
    
    num_segments = len(dsl_fold)
    ratios = [dsl_fold[i]['inv_angstrom_ratio'] for i in range(num_segments)]

    # Extract binning info and generate the correct energy values
    E_min, E_step, E_max  = map(float, dsl_fold[0]['Dimension3Binning'].split(','))
    num_y_bins = sim_folded_data[0].shape[1]  # Assuming (x_bins, y_bins)
    E_values = np.linspace(E_min, E_max, num_y_bins + 1)  # Use num_y_bins + 1 for bin edges

    # Add space for the colorbar
    fig, axes = plt.subplots(1, num_segments, gridspec_kw={'width_ratios': ratios}, figsize=(12, 5), subplot_kw={'projection':'mantid'}) #constrained_layout=True
    
    colormesh_pars = {
        'norm': SymLogNorm(linthresh=vmi, vmin=vmi, vmax=vma),
        'cmap': cma
    }
    
    for i in range(num_segments):
        x_start = 0
        x_end = sim_folded_data[i].shape[0]
        
        if i == 0:
            ax1 = axes[i]
            ax1.pcolormesh(np.arange(sim_folded_data[i].shape[0] + 1), E_values, sim_folded_data[i].T, rasterized=True, **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            #ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction='in')
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(np.arange(sim_folded_data[i].shape[0] + 1), E_values, sim_folded_data[i].T, rasterized=True, **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel('')
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction='in')
        else:
            ax = axes[i]
            ax.pcolormesh(np.arange(sim_folded_data[i].shape[0] + 1), E_values, sim_folded_data[i].T, rasterized=True, **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel('')
            ax.set_xticks([x_start, x_end])
            set_seg_xlabels(ax, i, dsl_fold)
            ax.tick_params(direction='in')
    
    plt.subplots_adjust(wspace=0, hspace=0)

    # Create the colorbar
    cb_ax = fig.add_axes([.91,.124,.02,.754])
    cbar = fig.colorbar(mappable,orientation='vertical',cax=cb_ax)
    cbar.ax.tick_params(labelsize=15)

    plt.rcParams.update({'font.size': 10})

    return fig



def generate_BZ_Report_with_Sims(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold, folded_results, all_sims_in_BZ):
    pdf_filename = os.path.join(pathSQE_params['output directory'],f'Reports/BZ_Report_{BZ_offset}.pdf')
    
    with PdfPages(pdf_filename) as pdf:
        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ, ax_path_foldedBZ = plt.subplots(subplot_kw={'projection': 'mantid'}, figsize=(8, 6))

        pathNames = [ds['Name'] for ds in dsl_fold]
        sortedWorkspaceNames = [name + '.nxs' for name in pathNames]

        # Adaptive colorbars
        step = float(pathSQE_params['E bins'].split(',')[1])
        ind_aboveElasticLine = int(np.ceil(3 / step))
        pathMax = -1
        pathMin = 1e4

        # Find min/max percentiles across all slices
        for seg in sortedWorkspaceNames:
            slice = mtd[seg].getSignalArray()[:, 0, 0, ind_aboveElasticLine:]
            mask = ~np.isnan(slice) & (slice != 0)
            if np.any(mask):
                segMin = np.percentile(slice[mask], 10)
                segMax = np.percentile(slice[mask], 97)
                pathMax = max(pathMax, segMax)
                pathMin = min(pathMin, segMin)

        # Plot experimental folded data        
        fig_path_foldedBZ = plot_along_path_foldedBZ(foldedSegNames=sortedWorkspaceNames, dsl_fold=dsl_fold,
                                                    Ei=pathSQE_params['T and Ei conditions'][0][1], vmi=pathMin,
                                                    vma=pathMax, cma=pathSQE_params['cmap'])

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ.suptitle(f'Folded experimental data along path in BZ {BZ_offset}', fontsize=18)

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ)
        plt.close(fig_path_foldedBZ)


        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ2, ax_path_foldedBZ2 = plt.subplots(subplot_kw={'projection': 'mantid'}, figsize=(8, 6))

        # Plot simulated folded data
        sim_folded_data = [result[2] for result in folded_results]  # Extract final folded sim data
        simMax = -1
        simMin = 1e4
        # Find min/max percentiles across all slices
        for path_sim in sim_folded_data:
            mask = ~np.isnan(path_sim) & (path_sim > 0)
            if np.any(mask):
                segMin = np.percentile(path_sim[mask],60)
                segMax = np.max(path_sim[mask])
                simMax = max(pathMax, segMax)
                simMin = min(pathMin, segMin)

        fig_path_foldedBZ2 = plot_along_path_foldedBZ_sim(sim_folded_data,dsl_fold,pathSQE_params['T and Ei conditions'][0][1],simMin,simMax,pathSQE_params['cmap'])

        # Add BZ_offset information at the top of the final folded path figure
        fig_path_foldedBZ2.suptitle(f'Folded simulated data along path in BZ {BZ_offset}', fontsize=18)

        # Save the final folded path figure to the PDF
        pdf.savefig(fig_path_foldedBZ2)
        plt.close(fig_path_foldedBZ2)


        # Loop over path segments and create slice comparison pages
        for i in range(len(all_slice_names_in_BZ)):
            path_seg = [dsl_fold[i]['seg_start_name'], dsl_fold[i]['seg_end_name']]
            num_pages = int(np.ceil(len(all_slice_names_in_BZ[i]) / 6))  # 6 pairs (12 subplots) per page

            for page_num in range(num_pages):
                start_idx = page_num * 6
                end_idx = min((page_num + 1) * 6, len(all_slice_names_in_BZ[i]))

                fig, axes = plt.subplots(nrows=3, ncols=4, subplot_kw={'projection': 'mantid'}, figsize=(12, 9))
                plt.rcParams['figure.dpi'] = 300
                plt.rcParams['savefig.dpi'] = 300
                plt.set_cmap(pathSQE_params['cmap'])

                for j, slice_idx in enumerate(range(start_idx, end_idx)):
                    slice_info = all_slice_names_in_BZ[i][slice_idx]
                    slice_name = slice_info[0]
                    slice_ID = slice_info[1]
                    sim_data = all_sims_in_BZ[i][slice_idx]  # Extract final folded simulation data
                    row_index = j // 2  # 3 rows
                    col_index = (j % 2) * 2  # Exp goes in 0, 2; Sim goes in 1, 3

                    # --- Experimental Slice ---
                    step = float(pathSQE_params['E bins'].split(',')[1])
                    ind_aboveElasticLine = int(np.ceil(3 / step))
                    slice = mtd[slice_name].getSignalArray()[:, 0, 0, ind_aboveElasticLine:]
                    mask = ~np.isnan(slice) & (slice != 0)

                    if np.any(mask):
                        segMin = np.percentile(slice[mask], 10)
                        segMax = np.percentile(slice[mask], 97)
                        exp_plot = axes[row_index, col_index].pcolormesh(mtd[slice_name], vmin=segMin, vmax=segMax, rasterized=True)
                        axes[row_index, col_index].set_title('ID {}'.format(slice_ID))
                        if not all_slice_evals[i][slice_idx]:
                            axes[row_index, col_index].set_title('Flagged Slice (ID {})'.format(slice_ID), color='red')

                        # x ticks and labels
                        point1_str, point2_str = slice_name.replace("pathSQE_", "").split("_to_")
                        point1 = np.array([float(x) for x in point1_str.split(",")])
                        point2 = np.array([float(x) for x in point2_str.split(",")])
                        #print(slice_name, point1_str, point2_str, point1, point2)
                        x_min_exp, x_max_exp = axes[row_index, col_index].get_xlim()
                        #print(x_min_exp, x_max_exp)
                        axes[row_index, col_index].set_xticks([x_min_exp, x_max_exp])
                        axes[row_index, col_index].set_xticklabels([f'{tick}' for tick in [point1, point2]])

                        # y ticks and labels
                        exp_yticks = axes[row_index, col_index].get_yticks()
                        exp_yticklabels = [tick.get_text() for tick in axes[row_index, col_index].get_yticklabels()]
                        y_min_real, y_max_real = axes[row_index, col_index].get_ylim()
                        _, num_y_bins = sim_data.shape
                        sim_yticks = np.interp(exp_yticks, [y_min_real, y_max_real], [0, num_y_bins - 1])
                        sim_yticks = np.round(sim_yticks).astype(int)
                        sim_yticks = np.clip(sim_yticks, 0, num_y_bins - 1)
                        axes[row_index, col_index + 1].set_yticks(sim_yticks)
                        axes[row_index, col_index + 1].set_yticklabels([f'{float(tick.strip().replace("−", "-")):.2g}' for tick in exp_yticklabels])

                        # Copy axis labels
                        axes[row_index, col_index + 1].set_xlabel(axes[row_index, col_index].get_xlabel())
                        axes[row_index, col_index + 1].set_ylabel(axes[row_index, col_index].get_ylabel())

                        # --- Simulated Slice (Right of Exp) ---
                        mask = ~np.isnan(sim_data) & (sim_data > 0)

                        if np.any(mask):
                            segMin = np.percentile(sim_data[mask], 60)
                            segMax = np.percentile(sim_data[mask], 99)

                        if sim_data is not None:
                            sim_plot = axes[row_index, col_index + 1].pcolormesh(sim_data.T, vmin=segMin, vmax=segMax, rasterized=True)
                            x_min_sim, x_max_sim = axes[row_index, col_index + 1].get_xlim()
                            #print(x_min_sim, x_max_sim)
                            axes[row_index, col_index + 1].set_xticks([x_min_sim, x_max_sim])
                            axes[row_index, col_index + 1].set_xticklabels([f'{tick}' for tick in [point1, point2]])
                    else:
                        axes[row_index, col_index].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=10, color='red')
                        axes[row_index, col_index+1].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=10, color='red')


                fig.suptitle(f'BZ: {BZ_offset}, Path Segment: {path_seg[0]} to {path_seg[1]} - Page {page_num + 1}', fontsize=14)
                fig.tight_layout(h_pad=0.3, w_pad=0.3)

                pdf.savefig(fig)
                plt.close(fig)
