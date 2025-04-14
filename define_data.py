import os
import numpy as np

########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace 
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(T_Ei_conditions, **kwargs):
    shared_folder='/SNS/ARCS/IPTS-13861/shared/jen/reduction/'
    raw_data_folder='/SNS/ARCS/IPTS-13861/data/'
    mde_folder='/SNS/ARCS/IPTS-13861/shared/Aiden/merged_mde/'

    data_set_list=[]
    for T, E in T_Ei_conditions:
        data_set={'Runs':list(range(69587, 69788+1)), #not accurate since mde are already made
              'BackgroundRuns':None,   # Options: None;list of runs that are added together
              'RawDataFolder':raw_data_folder,      # Options:raw_data_folder string
              'MdeFolder':mde_folder,               # Options:mde_folder string
              'MdeName':'Ge_{}K_{}meV'.format(T,E),   # Options:mde_name string _UBcorrected _1deg_909to910 _5deg_904to914
              'BackgroundMdeName':'',      # Options:None;bg_mde_name string
              'MaskingDataFile':shared_folder+'van67625.nxs',         # Options:None;data_file_name; van_449014_unmasked_negative_detectors.nxs
              'NormalizationDataFile':shared_folder+'van67625.nxs',   # Options:None;data_file_name
              'SampleLogVariables':{'OmegaMotorName':'CCR16Rot','Temperature':T},   # found name by doing hdfview on ARCS_69913_event.nxs and checking under entry->DASlogs
              'UBSetup':{'a':5.66,'b':5.66,'c':5.66,'alpha':90,'beta':90,'gamma':90,'u':'1,0,0','v':'0,1,0'}, # from Dipanshu .m files
               #Data reduction options
              'Ei':None,                             # Options: None;Ei_somehow_determined
              'T0': None,                           # Options: None;T0_determined_from_mantid
              'BadPulsesThreshold':None,            # Options: None;bg_pulses_threshold value
              'TimeIndepBackgroundWindow':None,     # Options: None;[Tib_min,Tib_max]
              'E_min':None,                         # Options: None;Emin
              'E_max':None,                         # Options: None;Emax
              'AdditionalDimensions':None,          # Options: None;list of triplets ("name", min, max)
              'EfCorrectionFunction':None,          # Options:None;'HYSPEC_default_correction';Custom_Ef_Correction_Function_Name
              }
    data_set_list.append(data_set)


    return data_set_list

