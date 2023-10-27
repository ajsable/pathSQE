import numpy as np
from matplotlib.colors import LogNorm

def make_slices_aiden(T,E,extra=''):
    dsl=[]

    ## Elatic map in HH0 and 00L at various K=-HH0
    for K in np.arange(-1,1+1,step=0.5): # -1,1 and half ints -0.5,0.5
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-6,0.05,6', # good step based on 1D cut of bragg peak width in HH0
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'{},{}'.format(K-0.1,K+0.1), # good range based on 1D cut of bragg peak width in HH0
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'-8,0.075,8', # good step based on 1D cut of bragg peak width in L
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'-0.5,0.5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=1e-2)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'[00L]',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Elastic Map w/ -HH0=K={}'.format(T,E,K),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_elasticwK{}'.format(T,E,K),#;-x,y,z;x,y,-z;-x,y,-z',
                     'Int_ticks':'True'}
        #dsl.append(description)
    
    ## Vertical / out of plane elatic map in HH0 and -HH0 at various L
    for L in np.arange(-7,7+1): # -7,7
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-6,0.05,6', # good step based on 1D cut of bragg peak width in HH0
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-2,0.05,2', # good range based on 1D cut of bragg peak width in HH0
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'{},{}'.format(L-0.075,L+0.075), # good step based on 1D cut of bragg peak width in L
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'-0.5,0.5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=1e-2)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'[-HH0]',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Elastic Map w/ L={}'.format(T,E,L),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_elasticVert_HH0wL{}'.format(T,E,L),
                     'Int_ticks':'True'}
        #dsl.append(description)
    
    ## Vertical elatic map in 00L and -HH0 at various H=HH0
    for H in np.arange(-5,5+1,step=0.5): # -5,5 and half ints -4.5,4.5
        description={'QDimension0':'0,0,1',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'1,1,0',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-8,0.075,8', # good step based on 1D cut of bragg peak width in HH0
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-2,0.05,2', # good range based on 1D cut of bragg peak width in HH0
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'{},{}'.format(H-0.1,H+0.1), # good step based on 1D cut of bragg peak width in L
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'-0.5,0.5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=1e-2)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[00L]',
                                        'ytitle':'[-HH0]',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Elastic Map w HH0=H={}'.format(T,E,H),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_elasticVert_00LwH{}'.format(T,E,H),
                     'Int_ticks':'True'}
        #dsl.append(description)
    
    # 1d cut at 1,1,6 bragg peak along HH0 for bin sizes
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0.8,0.025,1.2',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.05,0.05',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'5.9,6.1',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-0.1,0.1',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 #'Smoothing':1,
                 #'Plot_parameters':{'norm':LogNorm(vmin=1e-1, vmax=1e2)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[HH0]',
                                    'ytitle':'Intensity',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n 1,1,6 Bragg Peak cut in [HH0]'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_116_HH0'.format(T,E),
                 'Int_ticks':'False'}
    #dsl.append(description)
    
    # 1d cut at 1,1,6 bragg peak along L for bin sizes
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'0.95,1.05',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.05,0.05',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'5.85,0.025,6.15',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-0.1,0.1',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 #'Smoothing':1,
                 #'Plot_parameters':{'norm':LogNorm(vmin=1e-1, vmax=1e2)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[00L]',
                                    'ytitle':'Intensity',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n 1,1,6 Bragg Peak cut in L'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_116_L'.format(T,E),#;-x,y,z;x,y,-z;-x,y,-z',
                 'Int_ticks':'False'}
    #dsl.append(description)
    
    # sweep over cuts in the 00L direction
    for H in np.arange(-5,5+1): # -5,5 from H min to H max for 20 meV coverage in elastic map
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'{},{}'.format(H-0.05,H+0.05),
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-0.05,0.05',
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'-8,0.075,8',
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'{},{},{}'.format(-0.5*E,0.01*E,E),
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=3e-5, vmax=3e-3)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[00L]',
                                        'ytitle':'Energy',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Along [00L] for [{},{},0]'.format(T,E,H,H),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_00LwHH{}'.format(T,E,H),
                     'Int_ticks':'True'}
        #dsl.append(description)
             
    # sweep over cuts in the HH0 direction
    for L in np.arange(-7,7+1): # -7,7 from H min to H max for 20 meV coverage in elastic map
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-5,0.05,5',
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-0.05,0.05',
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'{},{}'.format(L-0.075,L+0.075),
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'{},{},{}'.format(-0.5*E,0.01*E,E),
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=4e-5, vmax=2e-3)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'Energy',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Along [HH0] for [0,0,{}]'.format(T,E,L),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_HH0wL{}'.format(T,E,L),
                     'Int_ticks':'True'}
        #dsl.append(description)
    
    ## Q maps in HH0 and 00L but with various E bounds to look at 1/2 int bragg peaks 
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-6,0.05,6', # good step based on 1D cut of bragg peak width in HH0
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.1,0.1', # good range based on 1D cut of bragg peak width in HH0
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-8,0.075,8', # good step based on 1D cut of bragg peak width in L
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'-1,1', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=2e-3)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[HH0]',
                                    'ytitle':'[00L]',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Q map with E=[-1,1]'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_elasticSatur_Em1to1'.format(T,E),
                 'Int_ticks':'True'}
    #dsl.append(description)

    # 1D cuts to compare directly with IPTS-5307 and look for half int peaks
    for L in np.arange(-1,2+1): # -2,1 from H min to H max for q coverage in elastic map
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'0,0.025,5',
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-0.1,0.1',
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'{},{}'.format(L-0.1,L+0.1),
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'0,5',
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     #'Smoothing':1,
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'Intensity',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Along [HH0] w/ L={} and E=[0,5]'.format(T,E,L),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_HH0wE0to5andL{}'.format(T,E,L),
                     'Int_ticks':'True'}
        #dsl.append(description)

    ## Q maps in HH0 and -HH0 with varied L and E=[0,5] to look at 1/2 int bragg peaks 
    for L in np.arange(6,7+1): # -7,7
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-6,0.02,6', # good step based on 1D cut of bragg peak width in HH0
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-2,0.02,2', # good range based on 1D cut of bragg peak width in HH0
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'{},{}'.format(L-0.1,L+0.1), # good step based on 1D cut of bragg peak width in L
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'0,5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=2e-3)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'[-HH0]',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Q map with E=[0,5] and L={}'.format(T,E,L),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_inelastSatur_HKwL{}'.format(T,E,L),
                     'Int_ticks':'True'}
        #dsl.append(description)
        
    ## same as just above but slight customized specifically to share with Doug 
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-3.65,0.02,3.65', # good step based on 1D cut of bragg peak width in HH0
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-1.35,0.02,0.8', # good range based on 1D cut of bragg peak width in HH0
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'5.9,6.1', # good step based on 1D cut of bragg peak width in L
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'0,5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm(vmin=3e-4, vmax=5e-2)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[HH0]',
                                    'ytitle':'[-HH0]',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Inelastic Q map with L=[5.9,6.1] and E=[0,5]'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_inelastSatur_HKwL6_Doug'.format(T,E,L),
                 'Int_ticks':'True'}
    #dsl.append(description)
        
            ## Q maps in HH0 and 00L but with various E bounds to look at 1/2 int bragg peaks 
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-6,0.05,6', # good step based on 1D cut of bragg peak width in HH0
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'0.9,1.1', # good range based on 1D cut of bragg peak width in HH0
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-8,0.075,8', # good step based on 1D cut of bragg peak width in L
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'1,5', # choose E binning size based on resolution from PyChop (1/2 of resolution)
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm(vmin=5e-4, vmax=2e-3)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[HH0]',
                                    'ytitle':'[00L]',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Q map with E=[1,5]'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_elasticSatur_Em1to5_K1'.format(T,E),
                 'Int_ticks':'True'}
    #dsl.append(description)
    
    ## Q maps in HH0 and -HH0 with L=6 and various E to look for relationship btw smear peak position and E 
    for E_inelast in np.arange(-2,5+1): # 0,5      
        description={'QDimension0':'1,1,0',
                     'QDimension1':'-1,1,0',
                     'QDimension2':'0,0,1',
                     'Dimension0Name':'QDimension0',
                     'Dimension0Binning':'-3.65,0.02,3.65', # good step based on 1D cut of bragg peak width in HH0
                     'Dimension1Name':'QDimension1',
                     'Dimension1Binning':'-1.35,0.04,0.8', # good range based on 1D cut of bragg peak width in HH0
                     'Dimension2Name':'QDimension2',
                     'Dimension2Binning':'5.9,6.1', # good step based on 1D cut of bragg peak width in L
                     'Dimension3Name':'DeltaE',
                     'Dimension3Binning':'{},{}'.format(E_inelast,E_inelast+1), 
                     #'SymmetryOperations':'P 4/nmm',
                     'SymmetryOperations':'x,y,z',
                     'Smoothing':1,
                     'Plot_parameters':{'norm':LogNorm(vmin=1e-4, vmax=1e-2)},
                     'Axes_parameters':{'xrange':None,
                                        'yrange':None,
                                        'xtitle':'[HH0]',
                                        'ytitle':'[-HH0]',
                                        'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n Q map with L=6 and E=[{},{}]'.format(T,E,E_inelast,E_inelast+1),
                                        'aspect_ratio':'auto',
                                        'tight_axes':True},
                     'Name':'FeSi_{}K_{}meV_inelastSatur_HKwL6_E{}to{}'.format(T,E,E_inelast,E_inelast+1),
                     'Int_ticks':'True'}
        #dsl.append(description)
        
    # slice to match sqe sim
    description={'QDimension0':'1,1,0',
                 'QDimension1':'-1,1,0',
                 'QDimension2':'0,0,1',
                 'Dimension0Name':'QDimension0',
                 'Dimension0Binning':'-3.1,-2.9',
                 'Dimension1Name':'QDimension1',
                 'Dimension1Binning':'-0.05,0.05',
                 'Dimension2Name':'QDimension2',
                 'Dimension2Binning':'-1,0.075,1',
                 'Dimension3Name':'DeltaE',
                 'Dimension3Binning':'0,0.5,30',
                 #'SymmetryOperations':'P 4/nmm',
                 'SymmetryOperations':'x,y,z',
                 'Smoothing':1,
                 'Plot_parameters':{'norm':LogNorm(vmin=3e-5, vmax=3e-3)},
                 'Axes_parameters':{'xrange':None,
                                    'yrange':None,
                                    'xtitle':'[00L]',
                                    'ytitle':'Energy',
                                    'title':'FeSi ARCS, T = {} K, Ei = {} meV,\n [-3,-3,L]'.format(T,E),
                                    'aspect_ratio':'auto',
                                    'tight_axes':True},
                 'Name':'FeSi_{}K_{}meV_m3m3L'.format(T,E),
                 'Int_ticks':'True'}
    dsl.append(description)
             
    return dsl

