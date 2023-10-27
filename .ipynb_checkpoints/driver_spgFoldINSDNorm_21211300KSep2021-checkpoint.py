# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from mantid.geometry import SpaceGroupFactory
import math
from scipy.spatial import cKDTree
from functions_spgFoldINSDNorm_21211300KSep2021 import *

#### LOAD MD File (MultiDimensional) Takes a while######
import time 
start_time = time.time()

###### Define the sample parameters, paths, directory ###################
IPTS=21211
Instrument = 'ARCS'

data_dir='/SNS/{0}/IPTS-{1}/shared/Aiden/merged_mde/'.format(Instrument,IPTS)
mde_name='FeSi_10K_70meV_UBcorrected.nxs'
processed_vanadium='/SNS/{0}/IPTS-{1}/shared/Dipanshu/van108467_new.nxs'.format(Instrument,IPTS)
van=Load(processed_vanadium)

Opdir='/SNS/{0}/IPTS-{1}/shared/Aiden/comprehensive/BZfolding/'.format(Instrument,IPTS)
slicedir = Opdir+'ConstantESlices/'

########Let Things Load#####################################
merged3=LoadMD(data_dir+mde_name,LoadHistory=False) # Uncomment as necessary
M = mtd['merged3']

# below used to be uncommented, but I did the UB correction manually and do not want to override that work
#a,b,c,alpha,beta,gamma=4.449,4.449,4.449,90,90,90
#u,v= [-2.20298,-2.45527,-0.0922226],[-0.150688,0.258645,-3.2864]
#u = [-0.0562,-0.0797,1]
#v = [0.9809,1,0.1348]
#SetUB(merged3,a,b,c,alpha,beta,gamma,u,v)
#############End Load######################################
print("MD Loading Done in %s seconds" % (time.time() - start_time))                
#############End Load#######################################


######## END RUN THESE FILES SECOND #####################################
###### SOME INPUT FILES ########
### Input Primitve Cell and Spacegroup #####
Mat_name = "FeSi10K"
Temperature = 10
plattice = np.array([[1,0,0],[0,1,0],[0,0,1]])
# How to get the Symmetry Operations in the space group, example on 'R 3 c' (161)
print("Space group no. 198: {}".format(SpaceGroupFactory.subscribedSpaceGroupSymbols(198)))
spacegroup=SpaceGroupFactory.createSpaceGroup('P 21 3')
spacegroupstring = 'P 21 3'
#################################
#
################### BINNING DATA ##########################################
###### Define Data Ranges and ProjectionMatrx
kx2 = [-8,0.05,8]
kx = np.arange(kx2[0],kx2[2],kx2[1])
ky2 = [-8,0.05,8]
ky = np.arange(ky2[0],ky2[2],ky2[1])
kz2 = [-1.6,0.05,1.6]
kz = np.arange(kz2[0],kz2[2],kz2[1])
e2 = [0,0.5,67]
e = np.arange(e2[0],e2[2],e2[1])
np.save(Opdir+"Energy",e)
proju = np.array([0,0,1])
projv = np.array([1,1,0])
projw = np.cross(proju,projv)
ProjectionMatrix = np.transpose(np.array([proju,projv,projw]))

######### Loop through the energy axis and generate 3D Qvolumes
sqw = np.zeros(shape=(kx.shape[0],ky.shape[0],kz.shape[0],e.shape[0]))
sqwerror = np.zeros(shape=(kx.shape[0],ky.shape[0],kz.shape[0],e.shape[0]))
sqwnorm = np.zeros(shape=(kx.shape[0],ky.shape[0],kz.shape[0],e.shape[0]))


start_time2 = time.time()
for en in range(len(e)):
    makeSlice(M,van, ProjMatrix=(proju,projv,projw),BinParams=(kx2,ky2,kz2,[e[en]-e2[1]/2,e[en]+e2[1]/2]),Folding=None,Smoothing=0,OutputWorkspace='%smeV'%e[en]) # make the Slice
    
    wssignal = mtd['%smeV'%e[en]]  
    wsnorm = mtd['normMD']
    data=wssignal.getSignalArray()
    err2=wssignal.getErrorSquaredArray()
    norm=wsnorm.getSignalArray()
    
    #write file
    data2 = data.reshape(kx.shape[0],ky.shape[0],kz.shape[0])
    error2 = err2.reshape(kx.shape[0],ky.shape[0],kz.shape[0])
    norm2 = norm.reshape(kx.shape[0],ky.shape[0],kz.shape[0])
    test = np.nan_to_num(data2)
    test[np.where(test>1e300)] = 0 
    test2 = np.nan_to_num(error2)
    test2[np.where(test2>1e300)] = 0
    test3 = np.nan_to_num(norm2)
    test3[np.where(norm2>1e300)] = 0
    print(test.max())
    print(test2.max())
    print(test3.max())
    sqw[:,:,:,en] = test
    sqwerror[:,:,:,en] = test2
    sqwnorm[:,:,:,en] = test3
    plt.pcolor(kx,ky,np.log10(sqw[:,:,np.abs(kz-0).argmin(),en]).T,cmap="viridis")
    plt.savefig(slicedir+"{0} {1} S(Q,{2}).png".format(Mat_name,(proju+projv),e[en]))
    plt.close()
    plt.pcolor(kx,ky,np.log10(sqwerror[:,:,np.abs(kz-0).argmin(),en]).T,cmap="viridis")
    plt.savefig(slicedir+"{0} {1} squareerror(Q,{2}).png".format(Mat_name,(proju+projv),e[en]))
    plt.close()
    plt.pcolor(kx,ky,np.log10(sqwnorm[:,:,np.abs(kz-0).argmin(),en]).T,cmap="viridis")
    plt.savefig(slicedir+"{0} {1} norm(Q,{2}).png".format(Mat_name,(proju+projv),e[en]))
    plt.close()
    print("%s meV" %e[en])

###### Check Dimensions and range of MD File ########
sqwnormed = np.divide(sqw,sqwnorm,out=sqw, where=sqwnorm!=0)
sqwerrornorm2 = np.divide(sqwerror, sqwnorm*sqwnorm, out=sqwerror, where=sqwnorm*sqwnorm!=0)
sqwerrornorm = np.sqrt(sqwerrornorm2)
sqwsignaltonoise = np.divide(sqwnormed, sqwerrornorm, out=sqwnormed, where=sqwerrornorm!=0)
print("Original Highest Signal to noise %s" % sqwsignaltonoise.max())
print("Original Mean Signal to noise %s" % sqwsignaltonoise.mean())
sqwnorm[np.where(sqwnorm>0)] = 1
print("4D Binning Done in %s seconds" % (time.time() - start_time2))  
################### END BINNING DATA ################################
#np.save(Opdir+"Full Map",sqw)

#np.save(Opdir+"Full Mask",sqwnorm)
print(sqwnorm.shape)
################## FOLDING BEGINS #####################################
pbzcart,sbzebound,normbound,errorbound,kbx,kby,kbz = bzfoldold(sqw,sqwerror,sqwnorm,ProjectionMatrix,plattice,kx,ky,kz,e,Opdir)
np.save(Opdir+'%s_S(PBZ,E)'% Mat_name,sbzebound)
np.save(Opdir+'%s_normmap(PBZ,E)'% Mat_name,normbound)
np.save(Opdir+'%s_squareerror(PBZ,E)'% Mat_name,errorbound)
np.save(Opdir+'kbx',kbx)
np.save(Opdir+'kby',kby)
np.save(Opdir+'kbz',kbz)
np.save(Opdir+'qpbz',pbzcart)


#symmetrizedbz,symmetrizedbznorm,symmetricnorm,symerror,symerrornorm = symmetrizebz(spacegroup,pbzcart,sbzebound,spbznorm,normbound,errorbound,errornorm,kbx,kby,kbz,e,ProjectionMatrix)
symmetrizedbznorm,symerrornorm = symmetrizebz(spacegroup,pbzcart,sbzebound,normbound,errorbound,kbx,kby,kbz,e,ProjectionMatrix)


#np.save(Opdir+'%s_%s_S(PBZ,E)'%(Mat_name, spacegroupstring),symmetrizedbz)
np.save(Opdir+'%s_%s_Snorm(PBZ,E)'%(Mat_name, spacegroupstring),symmetrizedbznorm)
#np.save(Opdir+'%s_%s_squareerror(PBZ,E)'%(Mat_name, spacegroupstring),symerror)
np.save(Opdir+'%s_%s_Errornorm(PBZ,E)'%(Mat_name, spacegroupstring),symerrornorm)
#np.save(Opdir+'%s_%s_signaltonoise(PBZ,E)'% (Mat_name,spacegroupstring),symsignaltonoise)

symmetrizedbzbose = bose_norm(symmetrizedbznorm,e,Temperature)
np.save(Opdir+'%s_%s_Snorm{PBZ}E_BoseNorm'%(Mat_name, spacegroupstring),symmetrizedbzbose)
############ FOLDING ENDS ##########################################
#print(symmetrizedbzbose[:,:,np.abs(kbz-0).argmin(),np.abs(e-1).argmin()].max())
#print(symmetrizedbzbose.max())
#for i in kbz:
    #plt.pcolor(kbx,kby,(symmetrizedbzbose[:,:,np.abs(kbz-i).argmin(),np.abs(e-15).argmin()]).T,cmap="Blues")#,vmax=symmetrizedbzbose.max())
    #plt.show()
    #plt.close()
#### Load Files #####
#kbx = np.load(Opdir+'kbx.npy')
#kby = np.load(Opdir+'kby.npy')
#kbz = np.load(Opdir+'kbz.npy')
#e = np.load(Opdir+'Energy.npy')
#print(e)
#symmetrizedbzbose = np.load(Opdir+'%s_%s_Snorm{PBZ}E_BoseNorm.npy'%(Mat_name, spacegroupstring))

