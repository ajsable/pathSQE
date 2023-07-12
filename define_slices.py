import numpy as np
from matplotlib.colors import LogNorm

def make_slices_aiden(extra=''):
    dsl=[]
    
    ## trying dimensions that aren't 1 (e.g. h 2h 0)
    description={'QDimension0':'1,2,0',
                 'QDimension1':'2,-1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-4,0.025,4',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 
                 #'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm()},#'norm':LogNorm(vmin=1e-2,vmax=1e0)},
                 'Axes_parameters':{'xrange':(0,7),
                                    'yrange':(0,70),
                                    'xtitle':r'[H,2H,0] \AA^{-1}$',
                                    'ytitle':'Energy (meV)',
                                    'title':'ARCS, Ei = 70 meV, T = 300 K,\n [H,2H,0]',
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_300K_Ei70meV_H2H0'+extra}
    #dsl.append(description)
    
    ## trying dimensions that aren't integers (0.5H 0.5H 0)
    description={'QDimension0':'0.5,0.5,0',
                 'QDimension1':'0.5,-0.5,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-8,0.05,8.11',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 
                 #'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm()},#'norm':LogNorm(vmin=1e-2,vmax=1e0)},
                 'Axes_parameters':{'xrange':(0,7),
                                    'yrange':(0,70),
                                    'xtitle':r'[0.5H,0.5H,0] \AA^{-1}$',
                                    'ytitle':'Energy (meV)',
                                    'title':'ARCS, Ei = 70 meV, T = 300 K,\n [0.5,0.5H,0]',
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_300K_Ei70meV_0.5H0.5H0temp'+extra}
    #dsl.append(description)
   
    ## H H 0 for comparison with above
    description={'QDimension0':'1,1,0',
                 'QDimension1':'1,-1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-4,0.025,4',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 
                 #'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm()},#'norm':LogNorm(vmin=1e-2,vmax=1e0)},
                 'Axes_parameters':{'xrange':(0,7),
                                    'yrange':(0,70),
                                    'xtitle':r'[0.5H,0.5H,0] \AA^{-1}$',
                                    'ytitle':'Energy (meV)',
                                    'title':'ARCS, Ei = 70 meV, T = 300 K,\n [0.5,0.5H,0]',
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_300K_Ei70meV_HH0'+extra}
    #dsl.append(description)
    
    ## GAMMA>X segment for plotting along a path
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-0.15,0.15',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'0,0.025,0.5',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_GAMMA2X'+extra}
    dsl.append(description)
   
    ## X>M segment for plotting along a path
    description={'QDimension0':'1,0,0',
                 'QDimension1':'0,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.025,0.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'0.35,0.65',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_X2M'+extra}
    dsl.append(description)
   
    ## M>GAMMA segment for plotting along a path
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.025,0.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_M2GAMMA'+extra}
    dsl.append(description)
   
    ## GAMMA>R segment for plotting along a path
    description={'QDimension0':'1,1,1',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'-1,-1,2',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.025,0.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-0.15,0.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_GAMMA2R'+extra}
    dsl.append(description)
    
    ## R>X segment for plotting along a path
    description={'QDimension0':'1,0,1',
                 'QDimension1':'-1,0,1',
                 'QDimension2':'0,1,0',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0,0.025,0.5',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.15,0.15',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'0.35,0.65',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-35.5,1,66.5',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',#;-x,y,z;x,y,-z;-x,y,-z',
                 'Name':'FeSi_300K_Ei70meV_R2X'+extra}
    dsl.append(description)
   
    return dsl

