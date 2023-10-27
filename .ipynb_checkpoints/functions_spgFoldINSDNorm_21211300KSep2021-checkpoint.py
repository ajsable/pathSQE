# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import SpaceGroupFactory
import math
from scipy.spatial import cKDTree

######################################################################################################################
############################################## START FUNCTIONS ####################################################
######################################################################################################################

def makeSlice(MDSpaceName,Vanadium,ProjMatrix,BinParams,Folding='x,y,z',Smoothing=0,OutputWorkspace='result'):
    """
    Makes a beautiful plot from a slide or cut
    Returns result
    """
# ## Plot a Quick slice
    MDNorm(InputWorkspace=MDSpaceName,
           SolidAngleWorkspace=Vanadium, # This is where vanadium and p.c. normalization is performed
           QDimension0=ProjMatrix[0], QDimension1=ProjMatrix[1], QDimension2=ProjMatrix[2],
           Dimension0Name='QDimension0', Dimension0Binning=BinParams[0],
           Dimension1Name='QDimension1', Dimension1Binning=BinParams[1],
           Dimension2Name='QDimension2', Dimension2Binning=BinParams[2],
           Dimension3Name='DeltaE', Dimension3Binning=BinParams[3],
           #SymmetryOperations='P -3 m 1', # Symmetry by a space group
           SymmetryOperations=Folding, # Real space symmetry operations
           OutputWorkspace=OutputWorkspace,
           OutputDataWorkspace='dataMD', 
           OutputNormalizationWorkspace='normMD')
           
    if Smoothing!=0:
        print('Smoothing data...')
        SmoothMD(InputWorkspace='dataMD', WidthVector=Smoothing, Function='Gaussian', InputNormalizationWorkspace='normMD', OutputWorkspace='dataMD')
        SmoothMD(InputWorkspace='normMD', WidthVector=Smoothing, Function='Gaussian', InputNormalizationWorkspace='normMD', OutputWorkspace='normMD')
        DivideMD(LHSWorkspace='dataMD', RHSWorkspace='normMD', OutputWorkspace=OutputWorkspace)
    return BinParams

######################################################################################################################

def dim2array(d,center=True):
    """
    Create a numpy array containing bin centers along the dimension d
    input: d - IMDDimension
    return: numpy array, from min+st/2 to max-st/2 with step st  
    """
    dmin=d.getMinimum()
    dmax=d.getMaximum()
    dstep=d.getX(1)-d.getX(0)
    if center:
        return np.arange(dmin+dstep/2,dmax,dstep)
    else:
        return np.linspace(dmin,dmax,d.getNBins()+1)

######################################################################################################################

def bzfoldold(sqw,sqwerror,sqwnorm,ProjectionMatrix,plattice,kx,ky,kz,e,Opdir):
    import time
    start_time = time.time()

    ## Define Projection Matrix
    P = np.transpose(np.array([proju,projv,projw]))
    Pinv = np.linalg.inv(P)

    #Create Normalized Primitive Reciprocal Basis and its transpose:
    rplatticet = np.linalg.inv(np.transpose(plattice))
    rplattice = np.transpose(rplatticet)
    # FIND NON-EMPTY BZ'S.

    taus = []
    bx = np.arange(-20,20,1)
    by = np.arange(-20,20,1)
    bz = np.arange(-20,20,1)

    #GAMMA SEARCH
    taus.append([np.dot(rplattice,[x,y,z]) for x in bx for y in by for z in bz if kx[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]<kx[-1] and ky[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]<ky[-1]\
                     and kz[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]<kz[-1]])
    #tauspcart.append([np.dot(rplattice,[x,y,z]) for x in bx for y in by for z in bz if np.greater(np.sum(sqw[np.abs(kx-np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]).argmin(),\
    #                                                                                                   np.abs(ky-np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]).argmin(),\
    #                                                                                                   np.abs(kz-np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]).argmin(),:],axis=0),0)])
 
    #tauspcart.append([np.dot(rplattice,[x,y,z]) for x in bx for y in by for z in bz if kx[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]<kx[-1]\
    #and ky[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]<ky[-1] and kz[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]<kz[-1] and\
    #np.greater(np.sum(sqw[np.abs(kx-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0])).argmin(),\
    #np.abs(ky-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1])).argmin(),\
    #np.abs(kz-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2])).argmin(),:]),0)])

    #SUPERCELL SEARCH
    #tauspcart.append([np.dot(rplattice,[x,y,z]) for x in bx for y in by for z in bz if kx[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]<kx[-1]\
    #and ky[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]<ky[-1] and kz[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]<kz[-1] and\
    #np.greater(np.sum(sqw[np.abs(kx-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]+np.min(pbzhorace[:,0]))).argmin():\
    #np.abs(kx-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]+np.max(pbzhorace[:,0]))).argmin(),\
    #np.abs(ky-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]+np.min(pbzhorace[:,1]))).argmin():\
    #np.abs(ky-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]+np.max(pbzhorace[:,1]))).argmin(),\
    #np.abs(kz-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]+np.min(pbzhorace[:,2]))).argmin():\
    #np.abs(kz-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]+np.max(pbzhorace[:,2]))).argmin(),:]),0)])


    Tau = np.vstack(taus)
#    print(Taupcart)
    Tau = np.delete(Tau,np.where(np.isclose(Tau,[0,0,0]).all(axis=1)),0)
 
    # Find all nearest neighbors to the origin by computing the distances and finding the nearest neighbors:
    distTau = np.array([math.sqrt(t[0]**2+t[1]**2+t[2]**2) for t in Tau])
    nearestcart = Tau[np.where(distTau==np.min(distTau))]
    for i in range(len(nearestcart)):
        nearestcart[i][np.where(nearestcart[i]<0)]-=0.071
        nearestcart[i][np.where(nearestcart[i]>0)]+=0.071

    print(nearestcart)

    #Add the origin again:
    Tau = np.append(Tau,[[0,0,0]],axis=0)
    nearestcart = np.append(nearestcart,[[0,0,0]],axis=0)
    # Create a list of q-points inside the PBZ w.r.t co-ordinate P:
    nearesthorace = [np.dot(Pinv,i) for i in nearestcart]
    nearesthorace = np.vstack(nearesthorace)
    

    ########## HORACE PBZ ######################    
    # Define a grid spanning nearest neighbors.
    kxnearesthorace = kx[np.abs(kx-np.min(nearesthorace[:,0])).argmin()+1:np.abs(kx-np.max(nearesthorace[:,0])).argmin()]
    kynearesthorace = ky[np.abs(ky-np.min(nearesthorace[:,1])).argmin()+1:np.abs(ky-np.max(nearesthorace[:,1])).argmin()]
    kznearesthorace = kz[np.abs(kz-np.min(nearesthorace[:,2])).argmin()+1:np.abs(kz-np.max(nearesthorace[:,2])).argmin()]
    print(kxnearesthorace.shape)
    print(kynearesthorace.shape)
    print(kznearesthorace.shape)


    #Define PBZ:
    pbzcart = [np.dot(P,[i,j,k]) for i in kxnearesthorace for j in kynearesthorace for k in kznearesthorace\
               if np.equal(nearestcart[np.argmin(np.linalg.norm(np.dot(P,[i,j,k])-nearestcart,axis=1))],[0,0,0]).all()]


    pbzhorace = [np.dot(Pinv,i) for i in pbzcart]


    pbzhorace = np.vstack(pbzhorace)
    print(np.min(pbzhorace[:,0]),np.max(pbzhorace[:,0]))
    print(np.min(pbzhorace[:,1]),np.max(pbzhorace[:,1]))
    print(np.min(pbzhorace[:,2]),np.max(pbzhorace[:,2]))

    print("PBZ defined")

    #### Find non-empty Taus
    tauspcart = []
    #PBZ SEARCH
    tauspcart.append([np.dot(rplattice,[x,y,z]) for x in bx for y in by for z in bz if kx[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]<kx[-1]\
    and ky[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]<ky[-1] and kz[0]<np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]<kz[-1] and\
    np.greater(np.sum([sqw[np.abs(kx-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[0]+q[0])).argmin(),\
    np.abs(ky-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[1]+q[1])).argmin(),\
    np.abs(kz-(np.dot(Pinv,np.dot(rplattice,[x,y,z]))[2]+q[2])).argmin(),:] for q in pbzhorace]),0)])

    Taupcart = np.vstack(tauspcart)
    Taupcart = Taupcart
    np.save(Opdir+"Taupcart",Taupcart)
    Tauphorace = [np.dot(Pinv,i) for i in Taupcart]
    Tauphorace = np.vstack(Tauphorace)
    print("Dataset Tiled")
    print("#BZ's: %s" %Taupcart.shape[0])
 
    # Fold
    spbz = np.zeros(shape=(kxnearesthorace.shape[0],kynearesthorace.shape[0],kznearesthorace.shape[0],e.shape[0]))
    snumevents = np.zeros(shape=(kxnearesthorace.shape[0],kynearesthorace.shape[0],kznearesthorace.shape[0],e.shape[0]))
    serror = np.zeros(shape=(kxnearesthorace.shape[0],kynearesthorace.shape[0],kznearesthorace.shape[0],e.shape[0]))
    
    for q in pbzhorace:
        spbz[np.abs(kxnearesthorace-q[0]).argmin(),np.abs(kynearesthorace-q[1]).argmin(),np.abs(kznearesthorace-q[2]).argmin(),:] +=\
                    sum(sqw[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin(),:]/np.linalg.norm(np.dot(P,q+K))**2 for K in Tauphorace)
        snumevents[np.abs(kxnearesthorace-q[0]).argmin(),np.abs(kynearesthorace-q[1]).argmin(),np.abs(kznearesthorace-q[2]).argmin(),:] =\
                    sum(sqwnorm[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin()] for K in Tauphorace)
        serror[np.abs(kxnearesthorace-q[0]).argmin(),np.abs(kynearesthorace-q[1]).argmin(),np.abs(kznearesthorace-q[2]).argmin(),:] =\
                    sum(sqwerror[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin()] for K in Tauphorace)

    print("Folding Done")
    print(pbzhorace.shape)



    # Bound and Normalize
    kbx = kxnearesthorace[np.abs(kxnearesthorace-np.min(pbzhorace[:,0])).argmin():np.abs(kxnearesthorace-np.max(pbzhorace[:,0])).argmin()+1]
    kby = kynearesthorace[np.abs(kynearesthorace-np.min(pbzhorace[:,1])).argmin():np.abs(kynearesthorace-np.max(pbzhorace[:,1])).argmin()+1]
    kbz = kznearesthorace[np.abs(kznearesthorace-np.min(pbzhorace[:,2])).argmin():np.abs(kznearesthorace-np.max(pbzhorace[:,2])).argmin()+1]
    
    sbzebound = spbz[np.abs(kxnearesthorace-np.min(pbzhorace[:,0])).argmin():np.abs(kxnearesthorace-np.max(pbzhorace[:,0])).argmin()+1,\
                     np.abs(kynearesthorace-np.min(pbzhorace[:,1])).argmin():np.abs(kynearesthorace-np.max(pbzhorace[:,1])).argmin()+1,\
                     np.abs(kznearesthorace-np.min(pbzhorace[:,2])).argmin():np.abs(kznearesthorace-np.max(pbzhorace[:,2])).argmin()+1,]
    normbound = snumevents[np.abs(kxnearesthorace-np.min(pbzhorace[:,0])).argmin():np.abs(kxnearesthorace-np.max(pbzhorace[:,0])).argmin()+1,\
                     np.abs(kynearesthorace-np.min(pbzhorace[:,1])).argmin():np.abs(kynearesthorace-np.max(pbzhorace[:,1])).argmin()+1,\
                     np.abs(kznearesthorace-np.min(pbzhorace[:,2])).argmin():np.abs(kznearesthorace-np.max(pbzhorace[:,2])).argmin()+1,]
    errorbound = serror[np.abs(kxnearesthorace-np.min(pbzhorace[:,0])).argmin():np.abs(kxnearesthorace-np.max(pbzhorace[:,0])).argmin()+1,\
                     np.abs(kynearesthorace-np.min(pbzhorace[:,1])).argmin():np.abs(kynearesthorace-np.max(pbzhorace[:,1])).argmin()+1,\
                     np.abs(kznearesthorace-np.min(pbzhorace[:,2])).argmin():np.abs(kznearesthorace-np.max(pbzhorace[:,2])).argmin()+1,]

    print(kbx.shape)
    print(kby.shape)
    print(kbz.shape)
    print(sbzebound.shape)


    #spbznorm = np.divide(sbzebound, countbound, out=sbzebound, where=countbound!=0)

    print("Domain Weighing Done. Squareroot of the error done.")
    ##################################
    sbze = (pbzcart, sbzebound,normbound,errorbound,kbx,kby,kbz)
    print("Folding Done in %s seconds" % (time.time() - start_time))                

    return sbze

######################################################################################################################

def symmetrizebz(spacegroup,pbzcart,sbzebound,normbound,errorbound,kbx,kby,kbz,e,ProjectionMatrix):
    import time
    start_time = time.time()

    #3: Symmetrize by Fm-3m:
    print("pg")
    pg=spacegroup.getPointGroup()
    print(pg.getName())
    rotations = []
    for so in pg.getSymmetryOperations():
        rotations.append(np.array([so.transformHKL((1,0,0)),so.transformHKL((0,1,0)),so.transformHKL((0,0,1))]))
    P = ProjectionMatrix
    Pinv = np.linalg.inv(P)

    print(len(rotations))

    #Calculate Projection Norms:

    symsbze = np.zeros(shape=(sbzebound.shape[0],sbzebound.shape[1],sbzebound.shape[2],sbzebound.shape[3]))
    symnorm = np.zeros(shape=(normbound.shape[0],normbound.shape[1],normbound.shape[2],normbound.shape[3]))
    symerror = np.zeros(shape=(errorbound.shape[0],errorbound.shape[1],errorbound.shape[2],errorbound.shape[3]))

    for q in pbzcart:
        symsbze[np.abs(kbx-np.dot(Pinv,q)[0]).argmin(),np.abs(kby-np.dot(Pinv,q)[1]).argmin(),np.abs(kbz-np.dot(Pinv,q)[2]).argmin(),:]+=sum(sbzebound[np.abs(kbx-np.dot(Pinv,np.dot(R,q))[0]).argmin(),np.abs(kby-np.dot(Pinv,np.dot(R,q))[1]).argmin(),np.abs(kbz-np.dot(Pinv,np.dot(R,q))[2]).argmin(),:]\
                                                                                                                                            for R in rotations)
        symnorm[np.abs(kbx-np.dot(Pinv,q)[0]).argmin(),np.abs(kby-np.dot(Pinv,q)[1]).argmin(),np.abs(kbz-np.dot(Pinv,q)[2]).argmin(),]+=sum(normbound[np.abs(kbx-np.dot(Pinv,np.dot(R,q))[0]).argmin(),np.abs(kby-np.dot(Pinv,np.dot(R,q))[1]).argmin(),np.abs(kbz-np.dot(Pinv,np.dot(R,q))[2]).argmin(),:]\
                                                                                                                                             for R in rotations)
        symerror[np.abs(kbx-np.dot(Pinv,q)[0]).argmin(),np.abs(kby-np.dot(Pinv,q)[1]).argmin(),np.abs(kbz-np.dot(Pinv,q)[2]).argmin(),]+=sum(errorbound[np.abs(kbx-np.dot(Pinv,np.dot(R,q))[0]).argmin(),np.abs(kby-np.dot(Pinv,np.dot(R,q))[1]).argmin(),np.abs(kbz-np.dot(Pinv,np.dot(R,q))[2]).argmin(),:]\
                                                                                                                                             for R in rotations)

    symmetricbz=symsbze/len(rotations)
    symnorm/=len(rotations)
    symerror/=len(rotations)

    #symmetricbznorm = symmetricbz#np.multiply(symmetricbz, symnorm)
    symmetricbznorm = np.divide(symmetricbz, symnorm, out=symmetricbz, where=symnorm!=0)
    symmetricerrornorm = np.sqrt(np.divide(symerror, symnorm*symnorm, out=symerror, where=symnorm*symnorm!=0))    
    #symmetricsignaltonoise = np.divide(np.copy(symmetricbznorm), np.copy(symmetricerrornorm), out=np.copy(symmetricbznorm), where=np.copy(symmetricerrornorm)!=0)
    
    #print("Highest Fully Symmetrized signal to noise %s" % symmetricsignaltonoise.max())
    #print("Mean Fully Symmetrized signal to noise %s" % symmetricsignaltonoise.mean())
    print("Symmetrization Done in %s seconds" % (time.time() - start_time))

    return (symmetricbznorm,symmetricerrornorm)

######################################################################################################################

def bose_norm(sbze,e,T):

    kb = 0.086173303 #meV/K
    
    sbzbose = sbze
    for en in e:
        sbzbose[:,:,:,np.abs(e-en).argmin()] *= (1-math.exp(-en/(kb*T))) #*en
    return sbzbose
######## END RUN THESE FILES FIRST #################################



######################################################################################################################
################################################### END FUNCTIONS ###################################################
######################################################################################################################




def bzfold(sqw,sqwerror,sqwnorm,ProjectionMatrix,plattice,kx,ky,kz,e,Opdir):
    import time
    start_time = time.time()

    ## Define Projection Matrix
    P = ProjectionMatrix
    Pinv = np.linalg.inv(P)

    #Create Normalized Primitive Reciprocal Basis and its transpose:
    rplatticet = np.linalg.inv(np.transpose(plattice))
    rplattice = np.transpose(rplatticet)
    rplatticeprojected = np.dot(Pinv,rplattice)

    #Define PBZ 
    # Use CKtree to find nearest gamma point and define PBZ:
    
    ## Create Reciprocal Supercell with size +- 1 rplattice
    px, py, pz = np.tensordot(rplattice,np.mgrid[-1:2, -1:2, -1:2],axes=[0, 0])
    nearest = np.c_[px.ravel(), py.ravel(), pz.ravel()]
    # Generate tree object 
    tree = cKDTree(nearest)
    gamma_region_id = tree.query([0, 0, 0])[1]
    # Grid the projected Supercell and convert to identity co-ordinate
    nearestprojected = np.vstack([np.dot(Pinv,i) for i in nearest])
    kxnearestproj = np.arange(np.min(nearestprojected[:,0]),np.max(nearestprojected[:,0])+0.025,0.025)
    #kx[np.abs(kx-np.min(nearestprojected[:,0])).argmin()-1:np.abs(kx-np.max(nearestprojected[:,0])).argmin()+1]
    kynearestproj = np.arange(np.min(nearestprojected[:,1]),np.max(nearestprojected[:,1])+0.025,0.025)
    #ky[np.abs(ky-np.min(nearestprojected[:,1])).argmin()-1:np.abs(ky-np.max(nearestprojected[:,1])).argmin()+1]
    kznearestproj = np.arange(np.min(nearestprojected[:,2]),np.max(nearestprojected[:,2])+0.025,0.025)
    #kz[np.abs(kz-np.min(nearestprojected[:,2])).argmin()-1:np.abs(kz-np.max(nearestprojected[:,2])).argmin()+1]
    nearestgrid = np.vstack([np.dot(P,[i,j,k]) for i in kxnearestproj for j in kynearestproj for k in kznearestproj])
    # Find q-points inside the PBZ
    pbz_id = tree.query(nearestgrid)[1]
    id_in_pbz = (pbz_id==gamma_region_id)
    pbzcart = nearestgrid[id_in_pbz]
    pbzproj = np.vstack([np.dot(Pinv,i) for i in pbzcart])
    print(pbzproj.shape)

    print(np.min(pbzproj[:,0]),np.max(pbzproj[:,0]))
    print(np.min(pbzproj[:,1]),np.max(pbzproj[:,1]))
    print(np.min(pbzproj[:,2]),np.max(pbzproj[:,2]))

    print("PBZ defined")
    
    # Generate Primitive lattice.
    Tausprojected = []
    bx = np.arange(-20,20,1)
    by = np.arange(-20,20,1)
    bz = np.arange(-20,20,1)
    
    Tausprojected.append([np.dot(rplatticeprojected,[x,y,z]) for x in bx for y in by for z in bz if kx[0]<np.dot(rplatticeprojected,[x,y,z])[0]<kx[-1]\
    and ky[0]<np.dot(rplatticeprojected,[x,y,z])[1]<ky[-1] and kz[0]<np.dot(rplatticeprojected,[x,y,z])[2]<kz[-1] and\
    np.greater(sum(sqw[np.abs(kx-(np.dot(rplatticeprojected,[x,y,z])[0]+q[0])).argmin(),\
    np.abs(ky-(np.dot(rplatticeprojected,[x,y,z])[1]+q[1])).argmin(),\
    np.abs(kz-(np.dot(rplatticeprojected,[x,y,z])[2]+q[2])).argmin(),:] for q in pbzproj),0).any()])
 
    Tausprojected = np.vstack(Tausprojected)
    Taupcart = np.vstack([np.dot(P,i) for i in Tausprojected])
    np.save(Opdir+"Taupcart",Taupcart)
    print("Dataset Tiled")
    print("#BZ's: %s" %Taupcart.shape[0])



    spbz = np.zeros(shape=(kxnearestproj.shape[0],kynearestproj.shape[0],kznearestproj.shape[0],e.shape[0]))
    serror = np.zeros(shape=(kxnearestproj.shape[0],kynearestproj.shape[0],kznearestproj.shape[0],e.shape[0]))
    snorm = np.zeros(shape=(kxnearestproj.shape[0],kynearestproj.shape[0],kznearestproj.shape[0],e.shape[0]))
    
    for q in pbzproj:
        spbz[np.abs(kxnearestproj-q[0]).argmin(),np.abs(kynearestproj-q[1]).argmin(),np.abs(kznearestproj-q[2]).argmin(),:] +=\
                    sum(sqw[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin(),:]/np.linalg.norm(np.dot(P,q+K))**2 for K in Tausprojected)
        serror[np.abs(kxnearestproj-q[0]).argmin(),np.abs(kynearestproj-q[1]).argmin(),np.abs(kznearestproj-q[2]).argmin(),:] =\
                    sum(sqwerror[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin()] for K in Tausprojected)
        snorm[np.abs(kxnearestproj-q[0]).argmin(),np.abs(kynearestproj-q[1]).argmin(),np.abs(kznearestproj-q[2]).argmin(),:] =\
                    sum(sqwnorm[np.abs(kx-(q+K)[0]).argmin(),np.abs(ky-(q+K)[1]).argmin(),np.abs(kz-(q+K)[2]).argmin()] for K in Tausprojected)

    print("Folding Done")



    kbx = kxnearestproj[np.abs(kxnearestproj-np.min(pbzproj[:,0])).argmin():np.abs(kxnearestproj-np.max(pbzproj[:,0])).argmin()]
    kby = kynearestproj[np.abs(kynearestproj-np.min(pbzproj[:,1])).argmin():np.abs(kynearestproj-np.max(pbzproj[:,1])).argmin()+1]
    kbz = kznearestproj[np.abs(kznearestproj-np.min(pbzproj[:,2])).argmin():np.abs(kznearestproj-np.max(pbzproj[:,2])).argmin()+1]
    
    sbzebound = spbz[np.abs(kxnearestproj-np.min(pbzproj[:,0])).argmin():np.abs(kxnearestproj-np.max(pbzproj[:,0])).argmin(),\
                     np.abs(kynearestproj-np.min(pbzproj[:,1])).argmin():np.abs(kynearestproj-np.max(pbzproj[:,1])).argmin()+1,\
                     np.abs(kznearestproj-np.min(pbzproj[:,2])).argmin():np.abs(kznearestproj-np.max(pbzproj[:,2])).argmin()+1,]
    errorbound = serror[np.abs(kxnearestproj-np.min(pbzproj[:,0])).argmin():np.abs(kxnearestproj-np.max(pbzproj[:,0])).argmin(),\
                     np.abs(kynearestproj-np.min(pbzproj[:,1])).argmin():np.abs(kynearestproj-np.max(pbzproj[:,1])).argmin()+1,\
                     np.abs(kznearestproj-np.min(pbzproj[:,2])).argmin():np.abs(kznearestproj-np.max(pbzproj[:,2])).argmin()+1,]
    normbound = snorm[np.abs(kxnearestproj-np.min(pbzproj[:,0])).argmin():np.abs(kxnearestproj-np.max(pbzproj[:,0])).argmin(),\
                     np.abs(kynearestproj-np.min(pbzproj[:,1])).argmin():np.abs(kynearestproj-np.max(pbzproj[:,1])).argmin()+1,\
                     np.abs(kznearestproj-np.min(pbzproj[:,2])).argmin():np.abs(kznearestproj-np.max(pbzproj[:,2])).argmin()+1,]

    #sbzenorm = np.divide(sbzebound, normbound, out=sbzebound, where=normbound!=0)
    print(sbzebound.shape)
    ##################################

    sbze = (pbzcart,sbzebound,normbound,errorbound,kbx,kby,kbz)
    print("Folding Done in %s seconds" % (time.time() - start_time))                

    return sbze