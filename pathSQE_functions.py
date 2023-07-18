import numpy as np
import seekpath
from seekpath.util import atoms_num_dict



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



def find_qbins_1and2(qdim0, qdim1, qdim2, pt1, pt2):
    # we start with bins of 0 so hopefully can catch it if code isn't working
    qdim1_range = 0
    qdim2_range = 0 

    # if any h,k,l of both points on the segment is not 0 but that same h,k,l of qdim0
    # is 0, then we need some offset in bins for qdim 1/2. E.g. for X>M qdim0=[1,0,0] but
    # both X ([0,0.5,0])and M ([0.5,0.5,0]) have nonzero k=0.5. Thus we need offset in bins 
    # (in this case qdim1=[0,1,0] so nonzero offset in qbin1) 
    qdim0_zeroInds = np.where(qdim0 == 0)[0]
    if len(qdim0_zeroInds) == 0:
        # the default qdim1 and 2 bin range about 0
        qdim1_range = np.array([-0.15,0.15])
        qdim2_range = np.array([-0.15,0.15]) 
        return qdim1_range, qdim2_range
    
    needed_offset = np.zeros((3,1))
    nonzero_hkl = []
    for ind in qdim0_zeroInds:
        if (pt1[ind] != 0) and (pt2[ind] != 0):
            needed_offset[ind]=pt1[ind]
            nonzero_hkl.append(ind)
    
    if not np.any(needed_offset):
        #print("No offset is needed")
        qdim1_range = np.array([-0.15,0.15])
        qdim2_range = np.array([-0.15,0.15])
        return qdim1_range, qdim2_range
        
    # MIGHT NOT WORK FOR KEEPING 2 CONST... NEED TO CHECK e.g. M>R
    for ind in nonzero_hkl: # can prob get rid of this and just put nonzero_hkl wherever ind is
        if qdim1[ind] != 0:
            qdim_coeff_for_bin1 = needed_offset[ind] / qdim1[ind]
            qdim1_range = np.array([0.7,1.3]) * qdim_coeff_for_bin1
            qdim2_range = np.array([-0.15,0.15])
        elif qdim2[ind] != 0:
            qdim_coeff_for_bin2 = needed_offset[ind] / qdim2[ind]
            qdim2_range = np.array([0.7,1.3]) * qdim_coeff_for_bin2
            qdim1_range = np.array([-0.15,0.15])
  
    return qdim1_range, qdim2_range



def find_qdim0_range(qdim0, qdim1, qdim2, qdim1_range, qdim2_range, pt1, pt2):
    # calculating default QDimension0 range and bin size
    # bounds chosen to cover first BZ and bin width is from most mantid examples
    def_range = np.array([0,0.025,0.5])

    if check_bin0_and_path_endpoints(qdim0, qdim1, qdim2, qdim1_range, qdim2_range, pt1, pt2, bin_range=def_range):
        return def_range
    else:
        rev_range = def_range[::-1] 
        rev_range[0] = -rev_range[0]
        if check_bin0_and_path_endpoints(qdim0, qdim1, qdim2, qdim1_range, qdim2_range, pt1, pt2, bin_range=rev_range):
            return rev_range
        else:
            return 0
        
    # scaling may not be necessary now since normalizing qdim0... later make width in A-1 and equal to instrument resolution
    # current method of fixed fraction step of BZ will break for larger unit cells (smaller BZ in A-1)
    #minDiffNonZero = np.min(np.absolute(qdim0[np.nonzero(qdim0)[0]]))
    # scaler makes sure resolution stays normal with weird qdim0
    #scaler = 1/minDiffNonZero
    #qdim0_range = def_range*scaler
    #return def_range



def check_bin0_and_path_endpoints(qdim0, qdim1, qdim2, qdim1_range, qdim2_range, pt1, pt2, bin_range):
    pt1_check = np.all(((bin_range[0]*qdim0 + np.mean(qdim1_range)*qdim1 + np.mean(qdim2_range)*qdim2) == pt1))
    pt2_check = np.all(((bin_range[2]*qdim0 + np.mean(qdim1_range)*qdim1 + np.mean(qdim2_range)*qdim2) == pt2))
    #print(pt1_check, pt2_check)
    return pt1_check and pt2_check



def choose_dims_and_bins(path_segment, point_coords, perp_to_path=False):
    
    pt1 = np.array(point_coords[path_segment[0]])
    pt2 = np.array(point_coords[path_segment[1]])
    diff = pt2 - pt1

    # calculating default QDimension0 axis direction and scaling
    qdim0 = diff / np.max(np.absolute(diff))

    # determine the directions of qdim1 and qdim2
    qdim1, qdim2 = find_qdim_1and2(qdim0,perp_to_path)
    qdim1_range, qdim2_range = find_qbins_1and2(qdim0, qdim1, qdim2, pt1, pt2)

    # calculating default QDimension0  bin range and bin size
    qdim0_range = find_qdim0_range(qdim0, qdim1, qdim2, qdim1_range, qdim2_range, pt1, pt2)

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



def make_slice_desc(q_dims_and_bins, path_seg):
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
                'Name':'slice_'+path_seg[0]+'_to_'+path_seg[1]}

    return slice_desc



def make_all_slice_descs(path_to_poscar):
    poscar = simple_read_poscar('files_for_nb/POSCAR')
    path = seekpath.get_path(structure=poscar)

    dsl=[]
    # APPEARS IT MAY NOT BE ABLE TO HANDLE DISCONTINUITIES IN PATH YET (e.g. R>X then R>M (path[5]))
    for i in range(5): #range(len(path['path'])):
        #print("on path seg {}".format(i))
        path_seg = path['path'][i]
        #print(path_seg)
        q_dims_and_bins = choose_dims_and_bins(path_seg,path['point_coords'],perp_to_path=True)
        slice_desc = make_slice_desc(q_dims_and_bins, path_seg)
        dsl.append(slice_desc)
    
    return dsl
    