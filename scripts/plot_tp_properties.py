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
from matplotlib.backends.backend_pdf import PdfPages
import datetime

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
parser.add_argument("--channel-map",                 type=str,               help="channel map to use. Options are APA, CRP, FiftyLChannelMap", default="APA") # thinking of NP04 as default
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

# # print recap of options
# print("#############################################")
# print("Selected options:")
# print(" - Files: " + str(args.files))
# print(" - Output folder: " + str(args.output_folder))
# print(" - Number of TPs: " + str(args.number_tps))
# print(" - Superimpose: " + str(args.superimpose))
# print(" - Views: " + str(args.view))
# print(" - Channel map: " + str(args.channel_map))
# print(" - All: " + str(args.all))
# print(" - Time peak: " + str(args.time_peak))
# print(" - Time over threshold: " + str(args.time_over_threshold))
# print(" - ADC integral: " + str(args.adc_integral))
# print(" - ADC peak: " + str(args.adc_peak))
# print(" - Channel: " + str(args.channel))
# print(" - Detid: " + str(args.detid))
# print(" - Show: " + str(args.show))
# print(" - Verbose: " + str(args.verbose))
# print("#############################################")

# save the settings printed before in an array of strings 
script_settings = []
script_settings.append(" - Files: " + str(args.files))
script_settings.append(" - Output folder: " + str(args.output_folder))
script_settings.append(" - Number of TPs: " + str(args.number_tps))
script_settings.append(" - Superimpose: " + str(args.superimpose))
script_settings.append(" - Views: " + str(args.view))
script_settings.append(" - Channel map: " + str(args.channel_map))
script_settings.append(" - All: " + str(args.all))
script_settings.append(" - Time peak: " + str(args.time_peak))
script_settings.append(" - Time over threshold: " + str(args.time_over_threshold))
script_settings.append(" - ADC integral: " + str(args.adc_integral))
script_settings.append(" - ADC peak: " + str(args.adc_peak))
script_settings.append(" - Channel: " + str(args.channel))
script_settings.append(" - Detid: " + str(args.detid))
script_settings.append(" - Show: " + str(args.show))
script_settings.append(" - Verbose: " + str(args.verbose))

print("Selected options:")
for i in range(len(script_settings)):
    print (script_settings[i])

# read the file(s) and create arrays of TPs, using the class in TriggerPrimitive.py
# the code will deduce if the file is coming from offline or online data

tps_lists = []
for tpFile_path in args.files:

    print("Reading file: " + tpFile_path)

    this_tps_list = save_tps_array(tpFile_path, args.number_tps)
    
    if args.verbose:
        print("TPs list: ", this_tps_list)
    
    tps_lists.append(this_tps_list)

for i in range(len(tps_lists)):
    print ("Number of TPs in list ", i, " : ", len(tps_lists[i]))

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
    
    print ("Number of channels to plot: ", len(channels_to_plot))
    
    if args.verbose:
        print ("Plotting channels: ", channels_to_plot)
    
    # Select TPs to plot basing on the selected view and channel map
    # create a new list of lists of TPs, keeping only the TPs that are in the selected view
    # and in the selected channels
    tps_to_plot = []
    for tps_list in tps_lists:
        if args.verbose:
            print ("Number of TPs in this list: ", len(tps_list))
        this_tps_to_plot = []
        for tp in tps_list:
            if tp['channel']%(len(this_channel_map[0])) in channels_to_plot:
                this_tps_to_plot.append(tp)
        print ("Number of TPs to plot in this file: ", len(this_tps_to_plot))
        tps_to_plot.append(this_tps_to_plot)
        
    
    print("Plotting time peak...")
    if args.time_peak:
        time_peak_plot_path = plotTimePeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
    
    print("Plotting time over threshold...")
    if args.time_over_threshold:
        time_over_threshold_plot = plotTimeOverThreshold(tps_to_plot,
                            args.files,
                            superimpose=args.superimpose,
                            quantile=1,
                            y_min=None, y_max=None,
                            output_folder=args.output_folder,
                            output_suffix=view,
                            show=args.show)

    print("Plotting ADC integral...")
    if args.adc_integral:
        adc_integral_plot = plotADCIntegral(tps_to_plot,
                        args.files,
                        superimpose=args.superimpose,
                        quantile=1,
                        y_min=None, y_max=None,
                        output_folder=args.output_folder,
                        output_suffix=view,
                        show=args.show) 

    print("Plotting ADC peak...")
    if args.adc_peak:
        adc_peak_plot = plotADCPeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print("Plotting channel...")
    if args.channel:
        channel_plot = plotChannel(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    x_min=None, x_max=None,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
                
    print ("Plotting detid...")
    if args.detid:
        det_id_plot = plotDetId(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print(" ")
    print ("Finished plotting! Products are in: " + args.output_folder + " and are named after the input file and variables plotted, plus the suffix (if passed).")

   
### Create a pdf report with the plots and the settings used

# open pdf file
pdf_report_name = '{}tp_properties_report.pdf'.format(args.output_folder)

with PdfPages(pdf_report_name) as pdf:
    
    # first page contains logo, title and settings
    plt.figure(figsize=(9, 12))
    plt.subplot(2, 1, 1)
    plt.imshow(plt.imread('../tools/dune_logo.jpg')) # could be handled better but this'll do
    plt.axis('off')  # Hide axes
    
    plt.subplot(2, 1, 2)
    plt.axis('off')  # Hide axes
    # set font to bold and size 16
    plt.text(0.2, 1,'TP Properties Report', ha='center', va='center', weight='bold', fontsize=16)
    # Add time of creation
    now = datetime.datetime.now()
    creation_time_text = 'Report created on: {}'.format(now.strftime("%Y-%m-%d %H:%M:%S"))
    plt.text(0.8, 0.9, creation_time_text, ha='center', va='bottom')
    plt.text(0.1, 0.85, 'Selected options:')
    vert_shift = 0
    for setting in script_settings:
        plt.text(0.1, 0.8 + vert_shift, setting, ha='left', va='top')
        vert_shift -= 0.05
    pdf.savefig()  # Save figure to PDF
    plt.close()

    # Add plots
     
    if args.time_peak:
        # Add time_peak_plot_path
        plt.imshow(plt.imread(time_peak_plot_path))
        plt.axis('off')
        pdf.savefig()  # Save figure to PDF
        plt.close()
    
    if args.time_over_threshold:
        # Add time_over_threshold_plot
        plt.imshow(plt.imread(time_over_threshold_plot))
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
    if args.adc_integral:
        # Add adc_integral_plot
        plt.imshow(plt.imread(adc_integral_plot))
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
        
    if args.adc_peak:
        # Add adc_peak_plot
        plt.imshow(plt.imread(adc_peak_plot))
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
        
    if args.channel:
        # Add channel_plot
        plt.imshow(plt.imread(channel_plot))
        plt.axis('off')
        pdf.savefig()
        plt.close()
        
        
    if args.detid:
        # Add det_id_plot
        plt.imshow(plt.imread(det_id_plot))
        plt.axis('off')
        pdf.savefig()
        plt.close()

print ("")
print ("PDF report created: {}".format(pdf_report_name))
