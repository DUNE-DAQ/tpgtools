"""
 * @file hdf5_converter.py
 * @brief Reads hdf5 files and converts them to txt or npy files.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""
import argparse
import sys
import numpy as np

# DAQ libraries
import daqdataformats
import detdataformats
import fddetdataformats
import trgdataformats
from hdf5libs import HDF5RawDataFile

# function to convert a TriggerPrimitive to a numpy array
def tp_to_numpy(tp):
    return np.array([tp.time_start, tp.time_over_threshold, tp.time_peak, tp.channel, tp.adc_integral, tp.adc_peak, tp.detid, tp.type, tp.algorithm, tp.version, tp.flag])

# convert a TriggerPrimitive to a numpy array.
# This function also orders the TPs by time_start
# instead of having them organized in fragments
def convert_tpstream_to_numpy(filename, n_records_to_read=-1, n_tps_to_convert=-1):

    hdf5_file = HDF5RawDataFile(filename)

    # Get all records (TimeSlice)
    records = hdf5_file.get_all_record_ids()
    # Set the number of records to process
    print(f'Input file: {filename}')
    print(f'Total number of records: {len(records)}')
    
    # if argument is not passed, read all records
    if n_records_to_read == -1:
        n_records_to_read = len(records)
        print(f'Number of records to read: {n_records_to_read}')
    elif n_records_to_read > len(records):
        n_records_to_read = len(records)
        print(f'Number of records to process is greater than the number of records in the file. Setting number of records to process to: {n_records_to_read}')
    else:
        print(f'Number of records to read: {n_records_to_read}')

    all_tps = []

    # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
    tp_size = trgdataformats.TriggerPrimitive.sizeof()
    for record in records[:n_records_to_read]:
        print(f'Processing record {record}') # TODO can add verbose option and put all these in it
        # Get all data sources
        datasets = hdf5_file.get_fragment_dataset_paths(record)
        # These are arrays since there can be multiple fragments per record
        n_tps = []
        n_tps_remaining = []
        next_tp_in_datasets = []
        fragments = []
        for dataset in datasets:
            # Get the fragment
            frag = hdf5_file.get_frag(dataset)
            fragments.append(frag)
            # Get the number of TPs
            n_tps_remaining.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            n_tps.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            # Get the data
            next_tp_in_datasets.append(trgdataformats.TriggerPrimitive(frag.get_data(0)).time_start)
        # this is here to order the TPs by time_start and not just as they are written in the hdf5 file
        while len(n_tps_remaining) > 0:
            index = np.argmin(next_tp_in_datasets)
            tp = trgdataformats.TriggerPrimitive(fragments[index].get_data(tp_size*(n_tps[index]-n_tps_remaining[index])))
            all_tps.append(tp_to_numpy(tp))
            if len(all_tps)==n_tps_to_convert:
                all_tps = np.array(all_tps)
                # assigning named dtypes to the columns of the matrix
                dt = np.dtype([('time_start', int), ('time_over_threshold', int), ('time_peak', int), ('channel', int), ('adc_integral', int), 
                               ('adc_peak', int), ('detid', int), ('type', int), ('algorithm', int), ('version', int), ('flag', int)])
                all_tps = np.rec.fromarrays([all_tps[:,0], all_tps[:,1], all_tps[:,2], all_tps[:,3], all_tps[:,4], 
                                             all_tps[:,5], all_tps[:,6], all_tps[:,7], all_tps[:,8], all_tps[:,9], all_tps[:,10]], dtype=dt)

                print(f"Final shape: {all_tps.shape}") # TODO put under verbose mode?
                return all_tps                

            n_tps_remaining[index] -= 1
            if n_tps_remaining[index] == 0:
                del n_tps_remaining[index]
                del next_tp_in_datasets[index]
                del n_tps[index]
                del fragments[index]
            else:
                next_tp_in_datasets[index] = trgdataformats.TriggerPrimitive(fragments[index].get_data(tp_size*(n_tps[index]-n_tps_remaining[index]))).time_start
                


    # all_tps = np.array(all_tps)
    # dt = np.dtype([('time_start', int), ('time_over_threshold', int), ('time_peak', int), ('channel', int), ('adc_integral', int), 
    #   ('adc_peak', int), ('detid', int), ('type', int), ('algorithm', int), ('version', int), ('flag', int)])
    # all_tps = np.rec.fromarrays([all_tps[:,0], all_tps[:,1], all_tps[:,2], all_tps[:,3], all_tps[:,4], 
    #   all_tps[:,5], all_tps[:,6], all_tps[:,7], all_tps[:,8], all_tps[:,9], all_tps[:,10]], dtype=dt)
    # print(f"Final shape: {all_tps.shape}") 
    # return all_tps


# def tpstream_hdf5_to_numpy_as_cpp(filename, n_records_to_read=-1, out_format=["txt"]):

#     hdf5_file = HDF5RawDataFile(filename)

#     # Get all records (TimeSlices)
#     records = hdf5_file.get_all_record_ids()
#     # Set the number of records to process
#     print(f'Input file: {filename}')
#     print(f'Total number of records: {len(records)}')
#     if n_records_to_read == -1:
#         n_records_to_read = len(records)
#         print(f'Number of records to process: {n_records_to_read}')
#     elif n_records_to_read > len(records):
#         n_records_to_read = len(records)
#         print(f'Number of records to process is greater than the number of records in the file. 
#             Setting number of records to process to: {n_records_to_read}')
#     else:
#         print(f'Number of records to process: {n_records_to_read}')

#     all_tps = []

#     # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
#     tp_size = trgdataformats.TriggerPrimitive.sizeof()
#     for record in records[:n_records_to_read]:
#         print(f'Processing record {record}') 
#         # Get all data sources
#         datasets = hdf5_file.get_fragment_dataset_paths(record)

#         for dataset in datasets:
#             # Get the fragment
#             frag = hdf5_file.get_frag(dataset)
#             # Get the number of TPs
#             n_tps=int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof())

#             for i in range(n_tps):
#                 tp = trgdataformats.TriggerPrimitive(frag.get_data(tp_size*i))
#                 # print_tp(tp)
#                 tp_np = tp_to_numpy(tp)
#                 all_tps.append(tp_np)



#     all_tps = np.array(all_tps)
#     print(f"Final shape: {all_tps.shape}")

#     return all_tps
