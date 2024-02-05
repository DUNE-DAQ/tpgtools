"""
 * @file image_creator.py
 * @brief Reads array of TPs and converts them to images.
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

# parse from command line the args.files and the number of tps to plot
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files",                nargs="+",      help="files to read")
parser.add_argument("-e", "--output-folder",        type=str,       help="output folder",   default="./")
parser.add_argument("-n", "--number-tps",            type=int,    help="number of tps to plot")
parser.add_argument("-s", "--superimpose",           action="store_true",    help="superimpose plots")
parser.add_argument("-a", "--all",                   action="store_true",    help="plot all the TP variables")
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
    


print("Plotting time peak...")
if args.time_peak:
    plotTimePeak(tps_lists, args.files,
                 superimpose=args.superimpose,
                 quantile=1,
                 y_min=None, y_max=None,
                 output_folder=args.output_folder,
                 show=args.show)

print("Plotting time over threshold...")
if args.time_over_threshold:
    plotTimeOverThreshold(tps_lists, args.files,
                          superimpose=args.superimpose,
                          quantile=1,
                          y_min=None, y_max=None,
                          output_folder=args.output_folder,
                          show=args.show)

print("Plotting ADC integral...")
if args.adc_integral:
    plotADCIntegral(tps_lists, args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    show=args.show)

print("Plotting ADC peak...")
if args.adc_peak:
    plotADCPeak(tps_lists, args.files,
                superimpose=args.superimpose,
                quantile=1,
                y_min=None, y_max=None,
                output_folder=args.output_folder,
                show=args.show)

print("Plotting channel...")
if args.channel:
    plotChannel(tps_lists, args.files,
                superimpose=args.superimpose,
                x_min=None, x_max=None,
                y_min=None, y_max=None,
                output_folder=args.output_folder,
                show=args.show)
            
print ("Plotting detid...")
if args.detid:
    plotDetId(tps_lists, args.files,
                superimpose=args.superimpose,
                y_min=None, y_max=None,
                output_folder=args.output_folder,
                show=args.show)

print(" ")
print ("Finished plotting! Products are in: " + args.output_folder + " and are named after the input file and variables plotted.")