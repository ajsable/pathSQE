import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfReader, PdfWriter
import seekpath
from seekpath.util import atoms_num_dict
from mantid.simpleapi import *
from mantid.geometry import SpaceGroupFactory

########################################################################################################
# All necessary functions for running pathSQE
# Author: Aiden Sable. Sept 2024.
########################################################################################################

# To run pathSQE_driver.py, you need the following files in addition to this script:
#       pathSQE_driver.py
#       pathSQE_input.py
#       define_data.py
#       reduce_data_to_MDE.py
#       slice_utils.py
#       POSCAR (if pathSQE_params['use_seeKpath_path']=True)
#
# Reach out to Aiden for copies of these files that are known to work properly with pathSQE



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



def calc_para_scores(qdim0):
    hkl = np.array([[1,0,0],[0,1,0],[0,0,1]])
    unit_qdim0 = qdim0 / np.linalg.norm(qdim0)
    # calculate degree of parallelism btw qdim0 and hkl vecs
    para_scores = np.zeros((3,1))
    for i in range(3):
        para_scores[i] = np.dot(unit_qdim0,hkl[i])
    
    return para_scores



def find_qdim_1and2(qdim0, perp_to_path=True):
    hkl = np.array([[1,0,0],[0,1,0],[0,0,1]])
    # if perp to path method is wanted
    if perp_to_path:
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
        else:
            raise ValueError("The selected symmetry path contains a segment which is greater than 3 dimensions") 
    else:
        # if perp to path isn't wanted, do it based on which 2 hkl are most perp (i.e. least para) to dim0
        para_scores = calc_para_scores(qdim0)
        # want two directions most perpen to qdim0 to be qdim1 and qdim2
        qdim1_ind = np.argmin(para_scores)
        qdim1 = hkl[qdim1_ind]
        # set qdim1 as highest so it won't be chosen for qdim2
        para_scores[qdim1_ind] = 2
        qdim2_ind = np.argmin(para_scores)
        qdim2 = hkl[qdim2_ind]
    
    return qdim1, qdim2



def find_qbins(pathSQE_params, qdim0, qdim1, qdim2, pt1, pt2, BZ_offset):
    # we start with bins of 0 so hopefully can catch it if code isn't working
    qdim1_range = 0
    qdim2_range = 0 

    # find the new path endpoint coords given change of basis from cart to qdim0,1,2
    CoB_mat = np.linalg.inv(np.array([qdim0,qdim1,qdim2]).T)
    path_start_transf = CoB_mat @ (pt1 + BZ_offset)
    path_end_transf = CoB_mat @ (pt2 + BZ_offset)
    
    # fixed step size of 0.025 for qdim0 from path start to end
    qdim0_range = np.array([path_start_transf[0], pathSQE_params['qdim0_step_size'], path_end_transf[0]])
    
    # might need to check bin/path endpoints...
    # old check_bin0_and_path_endpoints function
    
    # fixed integration range for qdim1 and qdim2
    qdim1_int_range = pathSQE_params['qdim1_int_range']
    qdim1_range = np.array([path_start_transf[1]-qdim1_int_range, path_end_transf[1]+qdim1_int_range])
    
    qdim2_int_range = pathSQE_params['qdim2_int_range']
    qdim2_range = np.array([path_start_transf[2]-qdim2_int_range, path_end_transf[2]+qdim2_int_range])
    
    # if qdim1 or qdim2 bin range gets too big raise error because code above may be wrong...
    #if (np.absolute(qdim1_range[1]-qdim1_range[0]) > 0.2) or (np.absolute(qdim2_range[1]-qdim2_range[0]) > 0.2):
    #    raise ValueError("Might be something weird with bins for qdim1 or dqim2")
    
    return qdim0_range, qdim1_range, qdim2_range



def choose_dims_and_bins(pathSQE_params, point1, point2, perp_to_path=True, BZ_offset=np.array([0,0,0])):
    diff = point2 - point1

    # calculating default QDimension0 axis direction and scaling
    qdim0 = diff / np.max(np.absolute(diff))

    # determine the directions of qdim1 and qdim2
    qdim1, qdim2 = find_qdim_1and2(qdim0, perp_to_path)
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
                'Dimension3Binning':pathSQE_params['E_bins'],
                'SymmetryOperations':'x,y,z',
                'Name':'pathSQE_'+conv_to_desc_string(point1)+'_to_'+conv_to_desc_string(point2),
                'qdim0_range':q_dims_and_bins[3],
                'seg_start_name':path_seg[0],
                'seg_end_name':path_seg[1],
                'inv_angstrom_ratio':np.linalg.norm(diff)/0.5,
                'good_slice':True}

    #slice_desc['Name'] = slice_desc['Name'].replace(',', '_').replace('.', '').replace('-','m')
    
    return slice_desc



def make_all_slice_descs(path_to_poscar, BZ_offset=np.array([0,0,0])):
    poscar = simple_read_poscar(path_to_poscar)
    path = seekpath.get_path(structure=poscar)

    dsl=[]
    # APPEARS IT MAY NOT BE ABLE TO HANDLE DISCONTINUITIES IN PATH YET (e.g. R>X then R>M (path[5]))
    for i in range(5): #range(len(path['path'])):
        path_seg = path['path'][i]

        pt1 = np.array(path['point_coords'][path_seg[0]])
        pt2 = np.array(path['point_coords'][path_seg[1]])
        q_dims_and_bins = choose_dims_and_bins(pt1, pt2, perp_to_path=True, BZ_offset=np.array([0,0,0]))   
        slice_desc = make_slice_desc(q_dims_and_bins, pt1, pt2, path_seg)
        dsl.append(slice_desc)
        
    return dsl



def find_matching_spacegroup(pathSQE_params):
    query_spacegroup = pathSQE_params['spacegroup']
    print(f'Searching for spacegroup {query_spacegroup}')
    spacegroup_list = mantid.geometry.SpaceGroupFactory.getAllSpaceGroupSymbols()
    query_spacegroup = ''.join(query_spacegroup.split())  # Remove spacing from the query_spacegroup

    for spacegroup in spacegroup_list:
        stripped_spacegroup = ''.join(spacegroup.split())  # Remove spacing from the spacegroup in the list
        if query_spacegroup == stripped_spacegroup:
            print(f"Matching spacegroup found: {spacegroup}")
            return spacegroup

    raise ValueError("No matching spacegroup found.")



def generate_unique_paths(mtd_spacegroup, path_seg, pt1, pt2):
    spacegroup = SpaceGroupFactory.createSpaceGroup(mtd_spacegroup)
    pg = spacegroup.getPointGroup()
    rotations = [np.array([so.transformHKL((1, 0, 0)), so.transformHKL((0, 1, 0)), so.transformHKL((0, 0, 1))])
                 for so in pg.getSymmetryOperations()]

    new_pt1s = np.dot(rotations, pt1)
    new_pt2s = np.dot(rotations, pt2)

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



##############################################################################################################################
############################################# FUNCTIONS FOR GENERATING BZ LISTS ##############################################
##############################################################################################################################



def gen_BZ_coverage_list(mde_data, H_bound, K_bound, L_bound, E_bound):
    from slice_utils_07142023 import make_slice
    
    # specify bounds of data you want to include in folding (BZ list will be found from these)

    # slice description with desired bounds to find data in each BZ
    slice_desc_totalIntE = {'QDimension0':'1,0,0',
                            'QDimension1':'0,1,0',
                            'QDimension2':'0,0,1',
                            'Dimension0Name':'QDimension0',
                            'Dimension0Binning':'{},1,{}'.format(H_bound[0]-0.5, H_bound[1]+0.5),
                            'Dimension1Name':'QDimension1',
                            'Dimension1Binning':'{},1,{}'.format(K_bound[0]-0.5, K_bound[1]+0.5),
                            'Dimension2Name':'QDimension2',
                            'Dimension2Binning':'{},1,{}'.format(L_bound[0]-0.5, L_bound[1]+0.5),
                            'Dimension3Name':'DeltaE',
                            'Dimension3Binning':'{},{}'.format(E_bound[0], E_bound[1]),
                            'Name':'pathSQE_totalIntE'}

    make_slice(mde_data[0], slice_desc_totalIntE)
    data_totalIntE = mtd[slice_desc_totalIntE['Name']].getSignalArray()

    # look for BZ with data and keep them in list
    BZ_H_vec = np.arange(H_bound[0],H_bound[1]+1,1)
    BZ_K_vec = np.arange(K_bound[0],K_bound[1]+1,1)
    BZ_L_vec = np.arange(L_bound[0],L_bound[1]+1,1)

    BZ_list = []
    for h in range(len(BZ_H_vec)):
        for k in range(len(BZ_K_vec)):
            for l in range(len(BZ_L_vec)):
                if ~np.isnan(data_totalIntE[h,k,l]) and data_totalIntE[h,k,l] > 0:
                    BZ_list.append(np.array([BZ_H_vec[h], BZ_K_vec[k], BZ_L_vec[l]]))

    return BZ_list



def gen_BZ_coverage_list_threshold(mde_data, H_bound, K_bound, L_bound, E_bound):
    from slice_utils_07142023 import make_slice
    
    # specify bounds of data you want to include in folding (BZ list will be found from these)

    # slice description with desired bounds to find data in each BZ
    slice_desc_BZcovThreshold = {'QDimension0':'1,0,0',
                            'QDimension1':'0,1,0',
                            'QDimension2':'0,0,1',
                            'Dimension0Name':'QDimension0',
                            'Dimension0Binning':'{},0.5,{}'.format(H_bound[0]-0.5, H_bound[1]+0.5),
                            'Dimension1Name':'QDimension1',
                            'Dimension1Binning':'{},0.5,{}'.format(K_bound[0]-0.5, K_bound[1]+0.5),
                            'Dimension2Name':'QDimension2',
                            'Dimension2Binning':'{},0.5,{}'.format(L_bound[0]-0.5, L_bound[1]+0.5),
                            'Dimension3Name':'DeltaE',
                            'Dimension3Binning':'{},5,{}'.format(E_bound[0], E_bound[1]),
                            'Name':'pathSQE_BZcovThreshold'}

    make_slice(mde_data[0], slice_desc_BZcovThreshold)
    data_BZcovThreshold = mtd[slice_desc_BZcovThreshold['Name']].getSignalArray()

    # look for BZ with data and keep them in list
    BZ_H_vec = np.arange(H_bound[0],H_bound[1]+1,1)
    BZ_K_vec = np.arange(K_bound[0],K_bound[1]+1,1)
    BZ_L_vec = np.arange(L_bound[0],L_bound[1]+1,1)

    BZ_list = []
    for h in range(len(BZ_H_vec)):
        for k in range(len(BZ_K_vec)):
            for l in range(len(BZ_L_vec)):
                BZ_data = data_BZcovThreshold[2*h:2*h+2, 2*k:2*k+2, 2*l:2*l+2, :]
                frac_nonNaN_nonZero = np.sum(np.logical_and(~np.isnan(BZ_data), BZ_data != 0)) / np.size(BZ_data)
                if frac_nonNaN_nonZero > 0.5:
                    BZ_list.append(np.array([BZ_H_vec[h], BZ_K_vec[k], BZ_L_vec[l]]))
    BZ_list = np.array(BZ_list) 

    return BZ_list



def gen_BZ_list_from_FCC_BraggPts():
    
    h_range = np.arange(-7,8,1)
    k_range = np.arange(-7,8,1)
    l_range = np.arange(-7,8,1)

    # bragg peaks in conv basis w fcc selection rules and constraints from elastic maps
    bragg_pt_conv_fcc = []
    for h in h_range:
        for k in k_range:
            for l in l_range:
                # selection rules
                if (h%2==0 and k%2==0 and l%2==0) or (h%2!=0 and k%2!=0 and l%2!=0):
                    # additional experimental constraints
                    if np.linalg.norm([h,k,l])<5.2 and np.abs(np.array([h,k,l])@np.array([-0.5,0.5,0]))<1 and h>-1 and h<4 and k>-1 and k<4 and l>-1 and l<6:
                        bragg_pt_conv_fcc.append(np.array([h,k,l]))
    bragg_pt_conv_fcc = np.array(bragg_pt_conv_fcc)    

    return bragg_pt_conv_fcc



def gen_BZ_list_from_BCC_BraggPts():
    
    h_range = np.arange(-4,2,1)
    k_range = np.arange(-4,2,1)
    l_range = np.arange(-5,5,1)

    # bragg peaks in conv basis w bcc selection rules and constraints from elastic maps
    bragg_pt_conv_bcc = []
    for h in h_range:
        for k in k_range:
            for l in l_range:
                # selection rules
                if (h+k+l)%2==0:
                    # additional experimental constraints
                    if np.linalg.norm([h,k,l]) < 5 and np.abs(np.array([h,k,l])@np.array([-0.5,0.5,0]))<=1:
                        bragg_pt_conv_bcc.append(np.array([h,k,l]))
    bragg_pt_conv_bcc = np.array(bragg_pt_conv_bcc)    


    return bragg_pt_conv_bcc



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
        plt.savefig('/SNS/ARCS/IPTS-5307/shared/Aiden/comprehensive/pathSQE_outputs/withFilters/'+slice_name+'.png')
        
        if 1 > 0.54:
            filter_evals.append(True)
        else:
            filter_evals.append(False)    


    # filter to get rid of M point artifacts
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
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from mantid import plots
    
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
            ax1.pcolormesh(mtd[foldedSegNames[i]], **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            ax1.set_ylim(0., Ei)
            ax1.set_xticks([x_start, x_end])
            set_seg_xlabels(ax1, i, dsl_fold)
            ax1.tick_params(direction='in')
        elif i == num_segments - 1:
            axL = axes[i]
            mappable = axL.pcolormesh(mtd[foldedSegNames[i]], **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel('')
            axL.set_xticks([x_start, x_end])
            set_seg_xlabels(axL, i, dsl_fold)
            axL.tick_params(direction='in')
        else:
            ax = axes[i]
            ax.pcolormesh(mtd[foldedSegNames[i]], **colormesh_pars)
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
    


def reorder_pdf(pdf_filename, pdf_filename_reordered):
    # Merge the last page to the front using PyPDF2
    with open(pdf_filename, 'rb') as original, open(pdf_filename_reordered, 'wb') as reordered:
        reader = PdfReader(original)
        writer = PdfWriter()

        # Add the last page to the front
        writer.add_page(reader.pages[-1])

        # Add the rest of the pages
        for page in reader.pages[:-1]:
            writer.add_page(page)

        # Save the reordered PDF
        writer.write(reordered)

    # Delete the original PDF
    os.remove(pdf_filename)



def generate_BZ_Report_morePages(pathSQE_params, BZ_offset, all_slice_names_in_BZ, all_slice_evals, dsl_fold):
    pdf_filename = pathSQE_params['saveDir']+'BZ_Reports/BZ_Report_{}.pdf'.format(BZ_offset)
    with PdfPages(pdf_filename) as pdf:
        # Create a new figure for the final folded path in given BZ
        fig_path_foldedBZ, ax_path_foldedBZ = plt.subplots(subplot_kw={'projection': 'mantid'}, figsize=(8, 6))

        # plot the final folded path
        pathNames = [ds['Name'] for ds in dsl_fold]
        sortedWorkspaceNames = pathNames.copy()
        sortedWorkspaceNames = [name + '.nxs' for name in sortedWorkspaceNames]
    
        # Adaptive colorbars
        step = float(pathSQE_params['E_bins'].split(',')[1])
        ind_aboveElasticLine = int(np.ceil(3/step))
        pathMax = -1
        pathMin = 1e4
        for seg in sortedWorkspaceNames:              
            slice = mtd[seg].getSignalArray()[:,0,0,ind_aboveElasticLine:]
            # Mask NaN and zero values
            mask = ~np.isnan(slice) & (slice != 0)
            segMin = np.percentile(slice[mask],10)
            segMax = np.percentile(slice[mask],95)
            if segMax > pathMax:
                pathMax = segMax
            if segMin < pathMin:
                pathMin = segMin
        
        fig_path_foldedBZ = plot_along_path_foldedBZ(foldedSegNames=sortedWorkspaceNames, dsl_fold=dsl_fold,
                                                    Ei=pathSQE_params['Ei'], vmi=pathMin,
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

                for j, slice_name in enumerate(all_slice_names_in_BZ[i][start_idx:end_idx]):
                    # Adjust the subplot position based on the number of rows and columns
                    row_index = j // 4
                    col_index = j % 4
                        
                    step = float(pathSQE_params['E_bins'].split(',')[1])
                    ind_aboveElasticLine = int(np.ceil(3/step))
                    slice = mtd[slice_name].getSignalArray()[:,0,0,ind_aboveElasticLine:]
                    # Mask NaN and zero values
                    mask = ~np.isnan(slice) & (slice != 0)
                    
                    if np.any(mask != 0):
                        segMin = np.percentile(slice[mask],10)
                        segMax = np.percentile(slice[mask],95)
                        axes[row_index, col_index].pcolormesh(mtd[slice_name], vmin=segMin,
                                                               vmax=segMax)
                        if not all_slice_evals[i][j]:
                            axes[row_index, col_index].set_title('Bad slice', color='red')
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
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.gridspec import GridSpec
    from mantid import plots
        
    T=pathSQE_params['Temp']
    Ebins=pathSQE_params['E_bins']
    Ei=pathSQE_params['Ei']
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
        SaveMD(SED_ws, Filename=pathSQE_params['saveDir']+'{}_folded_SED.nxs'.format(i))
        
        if i == 0:
            slice = SED_ws.getSignalArray()
            mask = ~np.isnan(slice) & (slice != 0)
            vmi = np.percentile(slice[mask],10)
            vma = np.percentile(slice[mask],95)
            colormesh_pars['norm']=SymLogNorm(linthresh=vmi,vmin=vmi,vmax=vma)#(linthresh=9e-6,vmin=9e-6,vmax=2e-3)
                        
            ax1 = plt.subplot(gs[i],projection='mantid')
            ax1.pcolormesh(SED_ws, **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            ax1.set_ylim(0.,Ei)
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
