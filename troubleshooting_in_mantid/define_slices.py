import numpy as np
from matplotlib.colors import LogNorm

def make_slices_aiden(extra=''):
    dsl=[]
   # manually doing the qdims and bins from PathSQE output to see what is going wrong with slices
   
   ## GAMMA>X segment
    description={'QDimension0':'0,0.5,0',
                 'QDimension1':'1,0,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.05,1',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Name':'FeSi_300K_Ei70meV_ts_GAMMA2X'+extra}#;-x,y,z;x,y,-z;-x,y,-z',
    dsl.append(description)
   
    ## X>M segment
    description={'QDimension0':'0.5,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.05,1',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_ts_X2M'+extra}
    dsl.append(description)
   
    ## M>GAMMA segment4-4.5]
    description={'QDimension0':'-0.5,-0.5,0',
                 'QDimension1':'0.707,-0.707,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.05,1',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_ts_M2GAMMA'+extra}
    dsl.append(description)
    
    ## GAMMA>R segment
    description={'QDimension0':'0.5,0.5,0.5',
                 'QDimension1':'-0.707,0.707,0',
                 'QDimension2':'-0.408,-0.408,0.816',
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
                 'Name':'FeSi_300K_Ei70meV_ts_GAMMA2R'+extra}
    dsl.append(description)
    
    ## R>X segment1.5-2]
    description={'QDimension0':'-0.5,0,-0.5',
                 'QDimension1':'0.707,0,-0.707',
                 'QDimension2':'0,-1,-0',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.05,1',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_ts_R2X'+extra}
    dsl.append(description)
   
    return dsl

