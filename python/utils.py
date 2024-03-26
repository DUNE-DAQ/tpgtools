"""
 * @file utils.py
 * @brief Misc python functions that don't fall under any other library.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
from hdf5_converter import convert_tpstream_to_numpy

# This function creates a channel map array, 
# where the first column is the channel number and the second column 
# is the plane view (0 for U, 1 for V, 2 for X).
# TODO understand if this can be avoided by reading from detchannelmaps
def create_channel_map_array(which_channel_map="APA"):
    '''
    :param which_channel_map, options are: 
     - APA for Horizontal Drift 
     - CRP for Vertical Drift (valid also for coldbox)
     - 50L for 50L detector, VD technology
    :return: channel map array
    '''
    # create the channel map array
    channel_map = [] # create empty object for older python versions that can't handle in-scope declaration
    if which_channel_map == "APA":
        channel_map = np.empty((2560, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 2559, 2560, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(800), np.ones(800), np.ones(960)*2))        
    elif which_channel_map == "CRP":
        channel_map = np.empty((3072, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 3071, 3072, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(952), np.ones(952), np.ones(1168)*2))
    elif which_channel_map == "FiftyLChannelMap":
        channel_map = np.empty((128, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 127, 128, dtype=int)
        # channel_map[:, 1] = np.concatenate(np.ones(49)*2, np.zeros(39), np.ones(40))
        channel_map[:, 1] = np.concatenate((np.zeros(39), np.ones(40), np.ones(49)*2))
    else:
        print("Unsupported channel map:", which_channel_map)
        return None
    return channel_map


# Function to save tps in a list
def save_tps_array(filename, max_tps=-1):
    if filename.endswith('.txt'):
        
        # offset to have the first TP at t=0, handle with care
        # time_shift = time_start[0] 
        # time_start -= time_shift
        # time_peak -= time_shift
        # time_start *= 16e-9 # convert to seconds, TODO choose if to have this or not
        # time_peak *= 16e-9 # convert to seconds, TODO choose if to have this or not
            
        # Create a structured array with column names
        dt = np.dtype([('time_start', float), 
                       ('time_over_threshold', float),
                       ('time_peak', float),
                       ('channel', int),
                       ('adc_integral', int),
                       ('adc_peak', int),
                       ('detid', int),
                       ('type', int),
                       ('algorithm', int),
                       ('version', int),
                       ('flag', int)])
        
        # Fill all_tps, that is a numpy array, with the arrays of the variables
        if max_tps == -1:
           print ("The number of TPs to save is not specified, all TPs will be saved") 
        all_tps = np.loadtxt(filename, dtype=dt, max_rows=max_tps, usecols=(0,1,2,3,4,5,6,7,8,9,10))
        
        # sort the list by time_start
        all_tps.sort(order='time_start')
        
        if max_tps == -1:
            max_tps = len(all_tps)
                
        print ("Your ", max_tps, " TPs have been saved in a numpy array and sorted by time start")
        print (" ")
        
        return all_tps
    elif filename.endswith('.hdf5'):
        
        # it's also possible to select a number of records to read, if -1 it will read all but only save max_tps
        all_tps = convert_tpstream_to_numpy(filename=filename, n_tps_to_convert=max_tps, n_records_to_read=-1)
        
        return all_tps
    else:
        print("Unsupported format:", filename.split('.')[-1])
        return None
