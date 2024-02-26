"""
 * @file create_image.py
 * @brief Creates 2D images of TPs, plotting channel vs time. Colors represent the charge of the TP.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import os
import argparse
import warnings
import gc

import sys
sys.path.append('../python/') 
from utils import save_tps_array, create_channel_map_array
from hdf5_converter import convert_tpstream_to_numpy 
from image_creator import save_image, show_image
from cluster_maker import make_clusters


# parser for the arguments
parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('-i', '--input-file',           type=str,                         help='Input file name')
parser.add_argument('-o', '--output-path',          type=str,    default='./',        help='path to save the image')
parser.add_argument('-n', '--n-tps',                type=int,    default=1,           help='number of tps to process')
parser.add_argument('-t', '--ticks-limit',          type=int,    default=3,           help='closeness in ticks to cluster TPs')
parser.add_argument('-c', '--channel-limit',        type=int,    default=1,           help='closeness in channels to cluster TPs')
parser.add_argument('-m', '--min-tps',              type=int,    default=3,           help='minimum number of TPs to create a cluster and then an image')
parser.add_argument('--channel-map',                type=str,    default="APA",       help='"APA", "CRP" or "50L"')
parser.add_argument('--save-clusters',                action='store_true',              help='write the clusters to a textfile')
parser.add_argument('--show',                       action='store_true',              help='show the image')
parser.add_argument('--fixed-size',                 action='store_true',              help='make the image size fixed')
parser.add_argument('--img-width',                  type=int,    default=200,          help='width of the image, if fixed size')
parser.add_argument('--img-height',                 type=int,    default=300,         help='height of the image, if fixed size')
parser.add_argument('--x-margin',                   type=int,    default=2,           help='margin in x')
parser.add_argument('--y-margin',                   type=int,    default=2,          help='margin in y')
parser.add_argument('-v', '--verbose',              action='store_true',              help='verbose mode')

args = parser.parse_args()
input_file              = args.input_file
output_path             = args.output_path
show                    = args.show
n_tps                   = args.n_tps
ticks_limit             = args.ticks_limit
channel_limit           = args.channel_limit
min_tps                 = args.min_tps
channel_map             = args.channel_map
save_clusters             = args.save_clusters
fixed_size              = args.fixed_size
width                   = args.img_width
height                  = args.img_height
x_margin                = args.x_margin
y_margin                = args.y_margin
verbose                 = args.verbose

# recap of selected options
print ("#############################################")
print("Selected options:")
print(" - Input file: " + str(input_file))
print(" - Output path: " + str(output_path))
print(" - Number of TPs: " + str(n_tps))
print(" - Ticks limit: " + str(ticks_limit))
print(" - Channel limit: " + str(channel_limit))
print(" - Minimum number of TPs to create a cluster: " + str(min_tps))
print(" - Channel Map: " + str(channel_map))
print(" - Save clusters: " + str(save_clusters))
print(" - Show image: " + str(show))
print(" - Fixed size: " + str(fixed_size))
print(" - Image width: " + str(width))
print(" - Image height: " + str(height))
print(" - X margin: " + str(x_margin))
print(" - Y margin: " + str(y_margin))
print(" - Verbose: " + str(verbose))
print ("#############################################")
print (" ")

# reading with this function already orders by time_start
all_TPs = save_tps_array(input_file, max_tps=n_tps) 

if verbose:
    print("TPs read from file: ", input_file)
    print("Number of TPs: ", len(all_TPs))
    print(" ")

# create channel map to distinguish induction and collection
my_channel_map = create_channel_map_array(which_channel_map=channel_map)
print("Channel map created: ", channel_map)

clusters = make_clusters(all_TPs, my_channel_map, 
                     ticks_limit=ticks_limit, 
                     channel_limit=channel_limit, 
                     min_tps_to_cluster=min_tps)

print("Number of clusters: ", len(clusters))

# Check how many views (U, V, X) are in this sample
total_channels = my_channel_map.shape[0]
n_views = np.unique(my_channel_map[all_TPs['channel'] % total_channels, 1]).shape[0]
print("Number of different views in these clusters: ", n_views)
print (" ")

# Prepare the output folder
if not os.path.exists(output_path):
    os.makedirs(output_path)
    print ("Created output folder: " + output_path)


print("Creating event display image...")
# Create the images
save_image(all_TPs, 
           my_channel_map, 
           output_path=output_path, 
           output_name="event_display",
           min_tps_to_create_img=min_tps, 
           make_fixed_size=True, 
           width=width*3, 
           height=height*3, 
           x_margin=x_margin,
           y_margin=y_margin,
           only_collection=False,
           show=show)
print("Done!")
print("Creating images for each cluster...")

for i, cluster in enumerate(clusters):
    save_image(cluster, 
               my_channel_map, 
               output_path=output_path, 
               output_name= "track" + str(i),
               min_tps_to_create_img=min_tps, 
               make_fixed_size=fixed_size, 
               width=width, 
               height=height, 
               x_margin=x_margin,
               y_margin=y_margin,
               only_collection=False,
               show=show)
print("Done!")
print(" ")
    
if save_clusters:
    output_clustersFile = os.path.basename(input_file)
    output_clustersFile = os.path.splitext(output_clustersFile)[0]
    output_clustersFile += '_clusters.txt'    
    print("Writing clusters to text file...")
    with open(output_clustersFile, 'w') as f:
        for i, cluster in enumerate(clusters):
            f.write('cluster ' + str(i) + ':\n')
            for tp in cluster:
                # f.write(str(this_tp) + '\n')
                this_tp = np.array([int(tp["time_start"]), 
                                    int(tp["time_over_threshold"]),
                                    int(tp["time_peak"]),
                                    int(tp["channel"]), 
                                    int(tp["adc_integral"]),
                                    int(tp["adc_peak"]),
                                    int(tp["detid"]), 
                                    int(tp["type"]),
                                    int(tp["algorithm"]),
                                    int(tp["version"]),
                                    int(tp["flag"])])
                f.write (str(this_tp) + '\n')
                # f.write(f"{tp["time_start"]} {tp[1]} {tp[2]} {tp[3]} {tp[4]} {tp[5]} {tp[6]} {tp[7]}\n")
                # f.write(f"{tp[0]} {tp[1]} {tp[2]} {tp[3]} {tp[4]} {tp[5]} {tp[6]} {tp[7]}\n")
            f.write('\n')
    print("Done!")
    print(" ")
