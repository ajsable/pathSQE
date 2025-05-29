import time
import numpy as np
from phonopy import load
from phonopy.spectrum.dynamic_structure_factor import atomic_form_factor_WK1995
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
from matplotlib.colors import ListedColormap, SymLogNorm
import scipy
from scipy import ndimage
from matplotlib import ticker
THz2meV = 4.1357  #, neut_dict
import scipy.linalg
from Resolution import Instrument, Resolution

start = time.time()


######################## START OF INPUTS #########################

## Q and E inputs ##
Q_start = np.array([0,8,0]) 
Q_end = np.array([4,8,0])

conv_to_prim = np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
Q_start = conv_to_prim@Q_start 
Q_end = conv_to_prim@Q_end

save_as = 'H80_sqe'

Q_steps = 150

E_min = 0
E_max = 40
E_step = 0.25


## resolution related ##
ElasticFWHM = 1.2  # in meV, elastic line energy FWHM if instrument, else this is the fwhm used
QResolution = 0.02 # in rlu, gauss sigma (i.e. ~FWHM/2.35), constant values only supported at this time

ResolutionType='instrument' #Can be either 
                            #'constant' - a constant resolution in energy transfer
                            #'instrument' - calculates the estimated energy resolution using a windsor approximation and instrument parameters
                            #'polynomial' - uses fourth order polynomial to calculate resolution (provided polynomial coefficients by InstrumentScientist)

InstrumentName='ARCS'  #Can be either ARCS, CNCS, HERIX, Generic; Used for estimated energy resolution convolution
                   #Generic is Gaussian Width
IncidentEnergy = 40  #Needed only when ResolutionType is set to Instrument. In this case energy max should not exceed incident energy, lest the values are meaningless.

PolynomialCoef = [ 1.02188553396572,-0.0594951491635698,0.000791281838645456,2.58258569481579e-05,-1.64952077361938e-07  ]  #coeffiecient required for polynomial fitting of resolution function


## simulation details ##
primitive_cell = [[1,0,0],[0,1,0],[0,0,1]]
supercell = [5,5,5]
forcesets_file = None
forceconstants_file = 'FORCE_CONSTANTS' 

Temperature = 10
# Neutron coherent scattering length can be found at https://www.ncnr.nist.gov/resources/n-lengths/
coh_scatter_length ={'Ge': 8.185}



######################## END OF INPUTS ########################


## The script below are mainly copied from Phonopy website
def run(phonon,
        Qpoints,
        temperature,
        atomic_form_factor_func=None,
        scattering_lengths=None):
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

    Qpoints = np.array(Qpoints)

    SandE = np.array([dsf.frequencies,dsf.dynamic_structure_factors])
    return SandE

if __name__ == '__main__':
    phonon = load(supercell_matrix=supercell,
                primitive_matrix=primitive_cell,
                unitcell_filename="POSCAR",
                force_sets_filename=forcesets_file,
                force_constants_filename=forceconstants_file
                )

    q_start = Q_start
    q_end = Q_end
    band = []
    Qpoints = []
    for i in range(Q_steps+1):
        band.append(list(q_start+(q_end-q_start)/Q_steps*i))
    Qpoints = band
    ### To avoid DW singularity at Q=[0,0,0]
    Qpoints_abs = scipy.linalg.norm(Qpoints,axis=1)
    zero_index = np.where(Qpoints_abs==0)[0]
    if zero_index.shape[0] !=0:
        Qpoints[zero_index[0]] = np.array([1e-6,1e-6,1e-6])
    
    #print(len(Qpoints))

    mesh = [11, 11, 11]
    phonon.run_mesh(mesh,
                    is_mesh_symmetry=False, # symmetry must be off
                    with_eigenvectors=True) # eigenvectors must be true
    temperature = Temperature

    # For INS, scattering length has to be given.
    # The following values is obtained at (Coh b)
    # https://www.nist.gov/ncnr/neutron-scattering-lengths-list
    output = run(phonon,
                Qpoints,
                temperature,
                scattering_lengths=coh_scatter_length) # At this point, we need to input coh scattering length manually. But it's good for isotopes.
    
    ## output has shape as (2,len(Qpoints),branches), the [0,:,:] is for frequency and the [1,:,:] is for SQE; The frequency is in THz unit
    for i in range(len(Qpoints)):
        output[0,i,:] *= THz2meV

    ## Bin the SQE ##
    ne = int(1 + (E_max - E_min)/E_step)
    Evec = [E_min,E_step,ne, E_max]
    evalues = np.arange(E_min+0.5*E_step,ne*E_step,E_step)    

    BinnedSQE=np.zeros((int(Q_steps+1),int(ne)))

    #create the resolution function object
    inst = Instrument(InstrumentName,ElasticFWHM,IncidentEnergy)
    resolution = Resolution(type=ResolutionType,inst=inst,sigmaq=QResolution,poly=PolynomialCoef)

    # Energy binning with convolution
    for ih in range(Q_steps+1):  # q-points
        for j in range(len(output[0, 0, :])):  # phonon branches
            E_phonon = output[0][ih][j]  # in meV (already converted from THz)
            intensity = output[1][ih][j]
            #print(ih, j, intensity)
            
            if ResolutionType == 'instrument':
                conv = resolution.Gauss(evalues, 0, E_phonon, 0, 1, 1)
            else:
                conv = resolution.Gauss(evalues, 0, E_phonon, 0)
            
            BinnedSQE[ih, :] += intensity * conv

    if ResolutionType == 'instrument':
        qvecs = np.array(Qpoints)  # Ensure qvecs are in array format
        for iE in range(BinnedSQE.shape[1]):
            profile = BinnedSQE[:, iE]  # shape (nq,)
            blurred = np.zeros_like(profile)

            for i in range(len(qvecs)):
                qcenter = qvecs[i]
                conv_q = resolution.GaussQ(qvecs, qcenter, qscale=1)  # adjust qscale as needed
                blurred[i] = np.sum(profile * conv_q)

            BinnedSQE[:, iE] = blurred

    np.save(save_as+'.npy',BinnedSQE.T)
    
    plt.figure()
    plt.pcolormesh(BinnedSQE.T, norm=SymLogNorm(linthresh=1, vmin=1, vmax=5e2))
    plt.savefig(save_as+'.png')
    plt.close()

    plt.figure()
    plt.plot(BinnedSQE[40])
    plt.savefig(save_as+'_1d.png')
    plt.close()

    print('done in ', time.time()-start)
