"""
 * @file plot_tp_properties.py
 * @brief Plots TP properties from a file or a list of files
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

import sys
sys.path.append("../python/")
from properties_plotter import *
from utils import create_channel_map_array

# parse from command line the args.files and the number of tps to plot
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files",                 nargs="+",              help="files to read")
parser.add_argument("-e", "--output-folder",         type=str,               help="output folder, remember to put / at the end",   default="./")
parser.add_argument("-n", "--number-tps",            type=int,               help="number of tps to plot", default=100)
parser.add_argument("-s", "--superimpose",           action="store_true",    help="superimpose plots")
parser.add_argument("--view",                        nargs="+",              help="default is allViews on same plot, to have separate select one or more among U, V, X", default=["allViews"])
parser.add_argument("--channel-map",                 type=str,               help="channel map to use", default="APA") # thinking of NP04 as default
parser.add_argument("-a", "--all",                   action="store_true",    help="plot all the TP variables",default=False)
parser.add_argument("-t", "--time-peak",             action="store_true",    help="plot time peak")
parser.add_argument("-o", "--time-over-threshold",   action="store_true",    help="plot time over threshold")
parser.add_argument("-i", "--adc-integral",          action="store_true",    help="plot adc integral")
parser.add_argument("-p", "--adc-peak",              action="store_true",    help="plot adc peak")
parser.add_argument("-c", "--channel",               action="store_true",    help="plot channel")
parser.add_argument("-d", "--detid",                 action="store_true",    help="plot detid")
parser.add_argument("--show",                        action="store_true",    help="show plots", default=False)
parser.add_argument("-v", "--verbose",               action="store_true",    help="verbose mode", default=False)
args = parser.parse_args()

if args.all:
    args.time_peak = True
    args.time_over_threshold = True
    args.adc_integral = True
    args.adc_peak = True
    args.channel = True
    args.detid = True

# print recap of options
print("#############################################")
print("Selected options:")
print(" - Files: " + str(args.files))
print(" - Number of TPs: " + str(args.number_tps))
print(" - Superimpose: " + str(args.superimpose))
print(" - Views: " + str(args.view))
print(" - Channel map: " + str(args.channel_map))
print(" - All: " + str(args.all))
print(" - Time peak: " + str(args.time_peak))
print(" - Time over threshold: " + str(args.time_over_threshold))
print(" - ADC integral: " + str(args.adc_integral))
print(" - ADC peak: " + str(args.adc_peak))
print(" - Channel: " + str(args.channel))
print(" - Detid: " + str(args.detid))
print(" - Show: " + str(args.show))
print(" - Verbose: " + str(args.verbose))
print("#############################################")

# read the file(s) and create arrays of TPs, using the class in TriggerPrimitive.py
# the code will deduce if the file is coming from offline or online data

tps_lists = []
for tpFile_path in args.files:

    print("Reading file: " + tpFile_path)

    this_tps_list = save_tps_array(tpFile_path, args.number_tps)
    
    if args.verbose:
        print("TPs list: ", this_tps_list)
    
    tps_lists.append(this_tps_list)
    
print ("Number of TPs lists: ", len(tps_lists))

print (tps_lists[0][0]['channel'])

# creating a channel map for the given detector, to know which channel is in which view    
this_channel_map = create_channel_map_array(args.channel_map)

for view in args.view:
    if view == "allViews":
        print ("Plotting all views on the same plot:")
        # select all channels from map
        channels_to_plot = this_channel_map[:, 0]
    elif view == "U":
        print ("Plotting U view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the U view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 0, 0]  
    elif view == "V":
        print ("Plotting V view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the V view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 1, 0]
    elif view == "X":
        print ("Plotting X view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the X view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 2, 0]
    else :
        print ("View not recognized, please choose from U, V, X or allViews (default)")
        sys.exit(1)  
    
    
    # print ("Plotting channels: ", channels_to_plot)
    
    # Select TPs to plot basing on the selected view and channel map
    # create a new list of lists of TPs, keeping only the TPs that are in the selected view
    # and in the selected channels
    tps_to_plot = []
    for tps_list in tps_lists:
        if args.verbose:
            print ("Number of TPs in this list: ", len(tps_list))
        this_tps_to_plot = []
        for tp in tps_list:
            if tp['channel'] in channels_to_plot:
                this_tps_to_plot.append(tp)
        tps_to_plot.append(this_tps_to_plot)
    
    print("Plotting time peak...")
    if args.time_peak:
        plotTimePeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
    
    print("Plotting time over threshold...")
    if args.time_over_threshold:
        plotTimeOverThreshold(tps_to_plot,
                            args.files,
                            superimpose=args.superimpose,
                            quantile=1,
                            y_min=None, y_max=None,
                            output_folder=args.output_folder,
                            output_suffix=view,
                            show=args.show)

    print("Plotting ADC integral...")
    if args.adc_integral:
        plotADCIntegral(tps_to_plot,
                        args.files,
                        superimpose=args.superimpose,
                        quantile=1,
                        y_min=None, y_max=None,
                        output_folder=args.output_folder,
                        output_suffix=view,
                        show=args.show) 

    print("Plotting ADC peak...")
    if args.adc_peak:
        plotADCPeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print("Plotting channel...")
    if args.channel:
        plotChannel(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    x_min=None, x_max=None,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
                
    print ("Plotting detid...")
    if args.detid:
        plotDetId(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print(" ")
    print ("Finished plotting! Products are in: " + args.output_folder + " and are named after the input file and variables plotted, plus the suffix (if passed).")