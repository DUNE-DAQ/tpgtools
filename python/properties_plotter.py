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

from utils import save_tps_array, create_channel_map_array
from hdf5_converter import convert_tpstream_to_numpy 

# Common options, can move elsewhere TODO
alpha = 0.4 # transparency of the histograms, lower is more opaque
grid_in_superimpose = False
grid_in_not_superimpose = False
separate_views= True # this creates three separate plots (U, V, X)
image_format='.png'
legend_properties = {'weight': 'bold', 'size': 'small'}


def plotTimePeak(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):
    
    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return
    
    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    time_peak_all_files = []
    for tps_file in tps_lists:
        time_peak_all_files += [tp['time_peak'] - tp['time_start'] for tp in tps_file]
    x_max = np.quantile(time_peak_all_files, quantile)
    # print (f"  Quantile {quantile} of time_peak is {x_max}, set as x_max")
    binsize = x_max/20 # could be handled in a smarter way, TODO
    
    # this is in case, for testing, time_start and time_peak are always the same
    if x_max == 0:
        print ("  All time_peak values are the same, setting x_max to 1")
        x_max = 1
        binsize = 1
    
    del time_peak_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        # if tps_file is empty, skip this file
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
        fig.hist(time_peak, bins=np.arange(-0.5, x_max + 0.5, binsize), label=label, alpha=alpha, edgecolor='black')

        if not superimpose:
            fig.set_title(f"Time Peak, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plot_path = f"{output_folder}{this_filename}_timePeak_{output_suffix}{image_format}"
            plt.savefig(plot_path)

    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plot_path = f"{output_folder}superimposed_timePeak_{output_suffix}{image_format}"
        plt.savefig(plot_path)
        
    del time_peak # free memory
    plt.close()
    
    # return the plot path
    return plot_path


def plotTimeOverThreshold(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):
    
    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return
    
    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    time_over_threshold_all_files = []
    for tps_file in tps_lists:
        time_over_threshold_all_files += [tp['time_over_threshold'] for tp in tps_file]
    x_max = np.quantile(time_over_threshold_all_files, quantile)
    
    del time_over_threshold_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
        fig.hist(time_over_threshold, bins=np.arange(-0.5, x_max + 0.5, x_max/5), label=label, alpha=alpha, edgecolor='black')

        if not superimpose:
            fig.set_title(f"Time over Threshold, file {this_filename}", fontweight='bold')
            if show:
                plt.show()
            plot_path = f"{output_folder}{this_filename}_timeOverThreshold_{output_suffix}{image_format}"
            plt.savefig(plot_path)

    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plot_path = f"{output_folder}superimposed_timeOverThreshold_{output_suffix}{image_format}"
        plt.savefig(plot_path)
        
    # free memory
    del time_over_threshold
    plt.close()
    
    # return the plot path
    return plot_path


def plotChannel(tps_lists, file_names, superimpose=False, x_min=0, x_max=None, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):

    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return
    
    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
       
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    channel_all_files = []
    for tps_file in tps_lists:
        channel_all_files += [tp['channel'] for tp in tps_file] 

    for i, tps_file in enumerate(tps_lists):
        
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
            plot_path = f"{output_folder}{this_filename}_channel_{output_suffix}{image_format}"
            plt.savefig(plot_path)
                        
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plot_path = f"{output_folder}superimposed_channel_{output_suffix}{image_format}"
        plt.savefig(plot_path)
    
    del channel_all_files # free memory
    del channel # free memory
    plt.close()
    
    # return the plot path
    return plot_path


def plotADCIntegral(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):
    
    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return
        
    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
    
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    adc_integral_all_files = []
    for tps_file in tps_lists:
        adc_integral_all_files += [tp['adc_integral'] for tp in tps_file]
    
    x_max = np.quantile(adc_integral_all_files, quantile)
    
    del adc_integral_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
            plot_path = f"{output_folder}{this_filename}_adcIntegral_{output_suffix}{image_format}"
            plt.savefig(plot_path)
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plot_path = f"{output_folder}superimposed_adcIntegral_{output_suffix}{image_format}"
        plt.savefig(plot_path)
        
    
    del adc_integral # free memory
    plt.close()
    
    # return the plot path
    return plot_path


def plotADCPeak(tps_lists, file_names, superimpose=False, quantile=1, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):
    
    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return

    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
        
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    adc_peak_all_files = []
    for tps_file in tps_lists:
        adc_peak_all_files += [tp['adc_peak'] for tp in tps_file]
    x_max = np.quantile(adc_peak_all_files, quantile)
    
    del adc_peak_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
            plot_path = f"{output_folder}{this_filename}_adcPeak_{output_suffix}{image_format}"
            plt.savefig(plot_path)
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        if show:
            plt.show()
        plot_path = f"{output_folder}superimposed_adcPeak_{output_suffix}{image_format}"
        plt.savefig(plot_path)
    
    plt.close()
    # return the plot path
    return plot_path


def plotDetId(tps_lists, file_names, superimpose=False, y_min=0, y_max=None, output_folder=None, output_suffix="", show=False):

    # stop this function if tps_lists is empty
    if len(tps_lists) == 0:
        print("No TPs found, exiting")
        return
    
    # stop if all elements in tps_lists are empty
    if all([len(tps_list) == 0 for tps_list in tps_lists]):
        print("No TPs found, exiting")
        return
        
    plt.figure()
    fig = plt.subplot(111)  # for when superimpose is true
    
    # compute x_max using quantile, considering all the files
    detid_all_files = []
    for tps_file in tps_lists:
        detid_all_files += [tp['detid'] for tp in tps_file]
    x_max = np.max(detid_all_files)
    
    del detid_all_files # free memory

    for i, tps_file in enumerate(tps_lists):
        
        if len(tps_file) == 0:
            print (f"No TPs passed to the plotting function from {file_names[i]}, skipping")
            continue
        
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
            plot_path = f"{output_folder}{this_filename}_detid_{output_suffix}{image_format}"
            plt.savefig(plot_path)
            if show:
                plt.show()
            
    if superimpose:
        fig.legend(prop=legend_properties)
        plt.grid(grid_in_superimpose)
        plot_path = f"{output_folder}superimposed_detid_{output_suffix}{image_format}"
        plt.savefig(plot_path)
        if show:
            plt.show()
    
    plt.close()
    
    # return the plot path
    return plot_path
