"""
 * @file group_maker.py
 * @brief Reads takes tps array and makes groups.
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
def create_channel_map_array(which_detector="APA"):
    '''
    :param which_detector, options are: 
     - APA for Horizontal Drift 
     - CRP for Vertical Drift (valid also for coldbox)
     - 50L for 50L detector, VD technology
    :return: channel map array
    '''
    # create the channel map array
    if which_detector == "APA":
        channel_map = np.empty((2560, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 2559, 2560, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(800), np.ones(800), np.ones(960)*2))        
    elif which_detector == "CRP":
        channel_map = np.empty((3072, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 3071, 3072, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(952), np.ones(952), np.ones(1168)*2))
    elif which_detector == "50L":
        channel_map = np.empty((128, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 127, 128, dtype=int)
        # channel_map[:, 1] = np.concatenate(np.ones(49)*2, np.zeros(39), np.ones(40))
        channel_map[:, 1] = np.concatenate((np.zeros(39), np.ones(40), np.ones(49)*2))
    return channel_map

# Function to save tps in a list
def save_tps_array(filename, max_tps):
    if filename.endswith('.txt'):
        # data = np.genfromtxt(filename, delimiter=' ', max_rows=max_tps)
        # # print (data.shape)
        # # print (data)
        
        # # data = data.transpose()
        # # print (data.shape)
        # # print (data)    
        

        # # TP variables are 11
        # if len(data) == 11:
        #     daq = True
        #     print ("File ", filename, " comes from DAQ, has correct number of variables")
        # else:
        #     print ("WARNING: TPs in this file have an unexpected number of variables:", len(data))
        #     print (" ")
        #     # return # could stop, but for now we continue

        # appo vectors just for clarity, could be avoided
        # THIS IS THE CORRECT ORDER OF THE VARIABLES
        # time_start = data[0]
        # time_over_threshold = data[1]
        # time_peak = data[2] 
        # channel = data[3]
        # adc_integral = data[4]
        # adc_peak = data[5]
        # detid = data[6]
        # type = data[7]
        # algorithm = data[8]
        # version = data[9]
        # flag = data[10]
        
        # del data
    
        # offset to have the first TP at t=0, not ok for new format
        # time_shift = time_start[0] 
        # time_start -= time_shift
        # time_peak -= time_shift
        # time_start *= 16e-9 # convert to seconds, TODO choose if to have this or not
        # time_peak *= 16e-9 # convert to seconds, TODO choose if to have this or not
            
        # Create a structured array with column names
        dt = np.dtype([('time_start', float), ('time_over_threshold', float), ('time_peak', float), ('channel', int), ('adc_integral', int), ('adc_peak', int), ('detid', int), ('type', int), ('algorithm', int), ('version', int), ('flag', int)])
        # Fill all_tps, that is a numpy array, with the arrays of the variables
        # all_tps = np.array([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10]], dtype=dt)
        all_tps = np.loadtxt(filename, dtype=dt, max_rows=max_tps, usecols=(0,1,2,3,4,5,6,7,8,9,10))
        # all_tps = np.ndarray([time_start, time_over_threshold, time_peak, channel, adc_integral, adc_peak, detid, type, algorithm, version, flag], dtype=dt)
        # all_tps = np.rec.fromarrays([time_start, time_over_threshold, time_peak, channel, adc_integral, adc_peak, detid, type, algorithm, version, flag], dtype=dt)
        # print the type of all_tps


        # delete the appo vectors
        # del time_start, time_over_threshold, time_peak, channel, adc_integral, adc_peak, detid, type, algorithm, version, flag
        
        # sort the list by time_start
        all_tps.sort(order='time_start')
        
        print (type(all_tps))
        print (all_tps.shape)
        print (all_tps)
                
        print ("Your ", max_tps, " TPs have been saved in a numpy array and sorted by time start")
        
        return all_tps
    elif filename.endswith('.hdf5'):
        
        # it's also possible to select a number of records to read, if -1 it will read all but only save max_tps
        tps_array = convert_tpstream_to_numpy(filename=filename, n_tps_to_convert=max_tps, n_records_to_read=-1)
        
        print (" ")
        print ("Saved ", len(tps_array), " TPs from file ", filename)
        print (" ")
        
        # move out to a global place? TODO
        dt = np.dtype([('time_start', float), ('time_over_threshold', float), ('time_peak', float), ('channel', int), ('adc_integral', float), ('adc_peak', float), ('detid', int), ('type', int), ('algorithm', int), ('version', int), ('flag', int)])
        
        all_tps = []
        for tps in tps_array:
            # all_tps.append(np.rec.fromarrays([tps_array[i][0], tps_array[i][1], tps_array[i][2], tps_array[i][3], tps_array[i][4], tps_array[i][5], tps_array[i][6], tps_array[i][7], tps_array[i][8], tps_array[i][9], tps_array[i][10]], dtype=dt))
            all_tps.append(np.array(tps, dtype=dt))

        del tps_array
        
        # sort the list by time_start is done in convert_tpstream_to_numpy
        print ("Your ", max_tps, " TPs have been saved in a numpy array, sorted by time start")
        
        return np.array(all_tps)
    else:
        print("Unsupported format:", filename.split('.')[-1])
        return None
