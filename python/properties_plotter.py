"""
 * @file properties_plotter.py
 * @brief Reads takes tps from file (txt or hdf5) and plots basic quantities.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import matplotlib.pyplot as plt

from utils import save_tps_as_list, create_channel_map_array
from hdf5_converter import convert_tpstream_to_numpy 

# Common options, can move elsewhere TODO
alpha = 0.4 # transparency of the histograms, lower is more opaque
grid_in_superimpose = False
grid_in_not_superimpose = False
image_format='.png'
legend_properties = {'weight': 'bold', 'size': 'small'}


def plotTimePeak(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    time_peak_all_files = []
    for tps_file in tps_lists:
        time_peak_all_files += [tp['time_peak'] - tp['time_start'] for tp in tps_file]
    x_max = np.quantile(time_peak_all_files, quantile)
    
    del time_peak_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        time_peak = []
        for tp in tps_file:
            time_peak.append(tp['time_peak'] - tp['time_start'])
        
        this_filename = file_names[i].split('/')[-1]

        
        label = f"Time Peak, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("Time Peak [ticks]")

        if y_min is not None:
            fig.set_ylim(bottom=y_min)
        if y_max is not None:
            fig.set_ylim(top=y_max)
 
        # bin size is optimized to have a number of bins depending on x_max, thus based on the quantile
        fig.hist(time_peak, bins=np.arange(-0.5, x_max + 0.5), label=label, alpha=alpha, edgecolor='black')

        if not superimpose:
            fig.set_title(f"Time Peak, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plt.savefig(f"{output_folder}{this_filename}_timePeak{image_format}")

    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plt.savefig(f"{output_folder}timePeak{image_format}")
        
    del time_peak # free memory

    return


def  plotTimeOverThreshold(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    time_over_threshold_all_files = []
    for tps_file in tps_lists:
        time_over_threshold_all_files += [tp['time_over_threshold'] for tp in tps_file]
    x_max = np.quantile(time_over_threshold_all_files, quantile)
    
    del time_over_threshold_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        time_over_threshold = [tp['time_over_threshold'] for tp in tps_file]
        this_filename = file_names[i].split('/')[-1]
     
        label = f"Time over Threshold, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("Time over Threshold [ticks]")

        if y_min is not None:
            fig.set_ylim(bottom=y_min)
        if y_max is not None:
            fig.set_ylim(top=y_max)
 
        # bin size is optimized to have a number of bins depending on x_max, thus based on the quantile
        fig.hist(time_over_threshold, bins=np.arange(-0.5, x_max + 0.5), label=label, alpha=alpha, edgecolor='black')

        if not superimpose:
            fig.set_title(f"Time over Threshold, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plt.savefig(f"{output_folder}{this_filename}_timeOverThreshold{image_format}")

    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plt.savefig(f"{output_folder}timeOverThreshold{image_format}")
        
    # free memory
    del time_over_threshold

    return


def plotChannel(tps_lists, file_names, superimpose=False, x_min=0, x_max=None, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    channel_all_files = []
    for tps_file in tps_lists:
        channel_all_files += [tp['channel'] for tp in tps_file] 

    for i, tps_file in enumerate(tps_lists):
        channel = [tp['channel'] for tp in tps_file]
        this_filename = file_names[i].split('/')[-1]
     
        label = f"Channel, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("Channel")

        if y_min is not None:
            fig.set_ylim
        if y_max is not None:
            fig.set_ylim(top=y_max)
        
        if x_max is None:
            fig.set_xlim(right=np.max(channel_all_files))
        else:    
            fig.set_xlim(right=x_max)
            
            
        # bin size is optimized to have a number of bins depending on x_max, thus based on the quantile
        fig.hist(channel, bins=100, label=label, alpha=alpha)
        
        if not superimpose:
            fig.set_title(f"Channel, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plt.savefig(f"{output_folder}{this_filename}_channel{image_format}")
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plt.savefig(f"{output_folder}channel{image_format}")
    
    del channel_all_files # free memory
    del channel # free memory

    return            


def plotADCIntegral(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    adc_integral_all_files = []
    for tps_file in tps_lists:
        adc_integral_all_files += [tp['adc_integral'] for tp in tps_file]
    
    x_max = np.quantile(adc_integral_all_files, quantile)
    
    del adc_integral_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        adc_integral = [tp['adc_integral'] for tp in tps_file]
        this_filename = file_names[i].split('/')[-1]
     
        label = f"ADC Integral, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("ADC Integral")

        if y_min is not None:
            fig.set_ylim(bottom=y_min)
        if y_max is not None:
            fig.set_ylim(top=y_max)
 
        # bin size is optimized to have a number of bins depending on x_max, thus based on the quantile
        # np.arange(-0.5, x_max + 0.5, (x_max+1)/300)
        fig.hist(adc_integral, bins=30, range=(-0.5,x_max+0.5), label=label, alpha=alpha, edgecolor='black')
        
        if not superimpose:
            fig.set_title(f"ADC Integral, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plt.savefig(f"{output_folder}{this_filename}_adcIntegral{image_format}")
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plt.savefig(f"{output_folder}adcIntegral{image_format}")
        
    
    del adc_integral # free memory
    
    return


def plotADCPeak(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    adc_peak_all_files = []
    for tps_file in tps_lists:
        adc_peak_all_files += [tp['adc_peak'] for tp in tps_file]
    x_max = np.quantile(adc_peak_all_files, quantile)
    
    del adc_peak_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        adc_peak = [tp['adc_peak'] for tp in tps_file]
        this_filename = file_names[i].split('/')[-1]
     
        label = f"ADC Peak, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("ADC Peak")

        if y_min is not None:
            fig.set_ylim(bottom=y_min)
        if y_max is not None:
            fig.set_ylim(top=y_max)
 
        # bin size is optimized to have a number of bins depending on x_max, thus based on the quantile
        fig.hist(adc_peak, bins=30, range=(-0.5, x_max + 0.5), label=label, alpha=alpha, edgecolor='black')
        
        if not superimpose:
            fig.set_title(f"ADC Peak, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plt.savefig(f"{output_folder}{this_filename}_adcPeak{image_format}")
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plt.savefig(f"{output_folder}adcPeak{image_format}")
        
    return


def plotDetId(tps_lists, file_names, superimpose=False, y_min=0, y_max=None, output_folder=None, show=False):
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    detid_all_files = []
    for tps_file in tps_lists:
        detid_all_files += [tp['detid'] for tp in tps_file]
    x_max = np.max(detid_all_files)
    
    del detid_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        detid = [tp['detid'] for tp in tps_file]
        this_filename = file_names[i].split('/')[-1]
     
        label = f"DetId, file {this_filename}"
        this_filename = this_filename.split('.')[0]
        
        if not superimpose:
            fig = plt.subplot(111)
            plt.grid(grid_in_not_superimpose) 
        fig.set_xlabel("DetId")
        fig.set_xticks(np.arange(0, x_max + 1, 1))

        if y_min is not None:
            fig.set_ylim(bottom=y_min)
        if y_max is not None:
            fig.set_ylim(top=y_max)
 
        # bin size from 0 to x_max
        fig.hist(detid, bins= x_max +1, label=label, alpha=alpha, edgecolor='black')
        
        if not superimpose:
            fig.set_title(f"DetId, file {this_filename}", fontweight='bold')
            plt.savefig(f"{output_folder}{this_filename}_detid{image_format}")
            if show:
                plt.show()
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        plt.savefig(f"{output_folder}detid{image_format}")
        if show:
            plt.show()
        
    return
