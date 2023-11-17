import numpy as np
import seekpath
from seekpath.util import atoms_num_dict
from mantid.simpleapi import *


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



def find_qdim_1and2(qdim0, perp_to_path=False):
    hkl = np.array([[1,0,0],[0,1,0],[0,0,1]])
    # if non-default perp to path method is wanted
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



def find_qbins(qdim0, qdim1, qdim2, pt1, pt2, BZ_offset):
    # we start with bins of 0 so hopefully can catch it if code isn't working
    qdim1_range = 0
    qdim2_range = 0 

    # find the new path endpoint coords given change of basis from cart to qdim0,1,2
    CoB_mat = np.linalg.inv(np.array([qdim0,qdim1,qdim2]).T)
    path_start_transf = CoB_mat @ (pt1 + BZ_offset)
    path_end_transf = CoB_mat @ (pt2 + BZ_offset)
    
    # fixed step size of 0.025 for qdim0 from path start to end
    qdim0_range = np.array([path_start_transf[0], 0.025, path_end_transf[0]])
    
    # might need to check bin/path endpoints...
    # old check_bin0_and_path_endpoints function
    
    # fixed integration range of +/-0.1 for qdim1 and qdim2
    qdim1_range = np.array([path_start_transf[1]-0.1, path_end_transf[1]+0.1])
    qdim2_range = np.array([path_start_transf[2]-0.1, path_end_transf[2]+0.1])
    
    # if qdim1 or qdim2 bin range gets too big raise error because code above may be wrong...
    #if (np.absolute(qdim1_range[1]-qdim1_range[0]) > 0.2) or (np.absolute(qdim2_range[1]-qdim2_range[0]) > 0.2):
    #    raise ValueError("Might be something weird with bins for qdim1 or dqim2")
    
    return qdim0_range, qdim1_range, qdim2_range



def choose_dims_and_bins(path_segment, point_coords, perp_to_path=False, BZ_offset=np.array([0,0,0])):
    pt1 = np.array(point_coords[path_segment[0]])
    pt2 = np.array(point_coords[path_segment[1]])
    diff = pt2 - pt1

    # calculating default QDimension0 axis direction and scaling
    qdim0 = diff / np.max(np.absolute(diff))

    # determine the directions of qdim1 and qdim2
    qdim1, qdim2 = find_qdim_1and2(qdim0,perp_to_path)
    qdim0_range, qdim1_range, qdim2_range = find_qbins(qdim0, qdim1, qdim2, pt1, pt2, BZ_offset)

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



def make_slice_desc(q_dims_and_bins, path_seg, point_coords):
    pt1 = np.array(point_coords[path_seg[0]])
    pt2 = np.array(point_coords[path_seg[1]])
    diff = pt2 - pt1
    
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
                'Dimension3Binning':'-35.5,1,66.5',
                'SymmetryOperations':'x,y,z',
                'Name':'pathSQE_'+path_seg[0]+'_to_'+path_seg[1],
                'qdim0_range':q_dims_and_bins[3],
                'seg_start_name':path_seg[0],
                'seg_end_name':path_seg[1],
                'inv_angstrom_ratio':np.linalg.norm(diff)/0.5}

    return slice_desc




def make_all_slice_descs(path_to_poscar, BZ_offset=np.array([0,0,0])):
    poscar = simple_read_poscar(path_to_poscar)
    path = seekpath.get_path(structure=poscar)

    dsl=[]
    # APPEARS IT MAY NOT BE ABLE TO HANDLE DISCONTINUITIES IN PATH YET (e.g. R>X then R>M (path[5]))
    for i in range(5): #range(len(path['path'])):
        path_seg = path['path'][i]
        q_dims_and_bins = choose_dims_and_bins(path_segment=path_seg, point_coords=path['point_coords'], perp_to_path=True, BZ_offset=BZ_offset)    
        slice_desc = make_slice_desc(q_dims_and_bins, path_seg, path['point_coords'])
        dsl.append(slice_desc)
        
    return dsl



def plot_along_path(dsl,vmi,vma,cma='jet'):
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.gridspec import GridSpec
    from mantid import plots
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams.update({'font.size': 28})

    fig=plt.figure(figsize=(12,5))

    colormesh_pars={}
    colormesh_pars['norm']=SymLogNorm(linthresh=vmi,vmin=vmi,vmax=vma)#(linthresh=9e-6,vmin=9e-6,vmax=2e-3)
    colormesh_pars['cmap']=cma
    
    num_segments = len(dsl)
    
    ratios = []
    for i in range(num_segments):
        ratios.append(dsl[i]['inv_angstrom_ratio'])
    
    gs = GridSpec(1, num_segments,width_ratios=ratios,wspace=0)

    for i in range(num_segments):
        x_start = dsl[i]['qdim0_range'][0]
        x_end = dsl[i]['qdim0_range'][2]
        
        if i == 0:
            ax1 = plt.subplot(gs[i],projection='mantid')
            ax1.pcolormesh(mtd[dsl[i]['Name']], **colormesh_pars)
            ax1.set_ylabel('E (meV)')
            ax1.set_xlabel('')
            ax1.set_ylim(0.,70.)
            ax1.set_xticks([x_start, x_end])
            ax1.set_xticklabels(['{}'.format(dsl[i]['seg_start_name']),''])
            ax1.tick_params(direction='in')
        elif i == num_segments-1:
            axL = plt.subplot(gs[i],sharey=ax1,projection='mantid')
            cb = axL.pcolormesh(mtd[dsl[i]['Name']], **colormesh_pars)
            axL.get_yaxis().set_visible(False)
            axL.set_xlabel('')
            axL.set_xticks([x_start, x_end])
            axL.set_xticklabels(['${}$'.format(dsl[i]['seg_start_name']), '${}$'.format(dsl[i]['seg_end_name'])])
            axL.tick_params(direction='in')
        else:
            ax = plt.subplot(gs[i],sharey=ax1,projection='mantid')
            ax.pcolormesh(mtd[dsl[i]['Name']], **colormesh_pars)
            ax.get_yaxis().set_visible(False)
            ax.set_xlabel('')
            ax.set_xticks([x_start, x_end])
            ax.set_xticklabels(['{}'.format(dsl[i]['seg_start_name']),''])
            ax.tick_params(direction='in')

    fig.colorbar(cb)
    plt.rcParams.update({'font.size': 10})

    return fig