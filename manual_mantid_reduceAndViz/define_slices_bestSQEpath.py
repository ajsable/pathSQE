import numpy as np
from matplotlib.colors import LogNorm

def make_slices_aiden(extra=''):
    dsl=[]
   # finding path with "best" BZ's manually to display on the ORNL GRSI poster about PathSQE. 
   # Looking for BZ w/ most q coverage and phonon dispersions present.
   
    ## GAMMA>X segment for plotting along a path [full range]
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-0.15,0.15',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-8,0.025,8',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Name':'FeSi_300K_Ei70meV_FR_GAMMA2X'+extra}#;-x,y,z;x,y,-z;-x,y,-z',
    dsl.append(description)
   
   ## GAMMA>X segment for plotting along a path [best BZ, around 1.5-2]
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-0.15,0.15',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'1.5,0.025,2',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Name':'FeSi_300K_Ei70meV_Best_GAMMA2X'+extra}#;-x,y,z;x,y,-z;-x,y,-z',
    dsl.append(description)
   
    ## X>M segment for plotting along a path [full range]
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-8,0.025,8',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'0.35,0.65',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_FR_X2M'+extra}
    dsl.append(description)
   
    ## X>M segment for plotting along a path [best BZ, around 1.5-2]
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'1.5,0.025,2',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'0.35,0.65',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_Best_X2M'+extra}
    dsl.append(description)
   
    ## M>GAMMA segment for plotting along a path [full range]
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-8,0.025,8',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_FR_M2GAMMA'+extra}
    dsl.append(description)
   
    ## M>GAMMA segment for plotting along a path [best BZ, around 4-4.5]
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'4,0.025,4.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_Best_M2GAMMA'+extra}
    dsl.append(description)
   
    ## GAMMA>R segment for plotting along a path [full range]
    description={'QDimension0':'1,1,1',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'-1,-1,2',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-8,0.025,8',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_FR_GAMMA2R'+extra}
    dsl.append(description)
    
    ## GAMMA>R segment for plotting along a path [best BZ, around 3-3.5]
    description={'QDimension0':'1,1,1',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'-1,-1,2',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'3,0.025,3.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_Best_GAMMA2R'+extra}
    dsl.append(description)
    
    ## R>X segment for plotting along a path [full range]
    description={'QDimension0':'1,0,1',
                 'QDimension1':'-1,0,1',
                 'QDimension2':'0,1,0',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-8,0.025,8',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'0.35,0.65',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_FR_R2X'+extra}
    dsl.append(description)
   
    ## R>X segment for plotting along a path [best BZ, around 1.5-2]
    description={'QDimension0':'1,0,1',
                 'QDimension1':'-1,0,1',
                 'QDimension2':'0,1,0',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'1.5,0.025,2',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'0.35,0.65',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_Best_R2X'+extra}
    dsl.append(description)
   
    return dsl

