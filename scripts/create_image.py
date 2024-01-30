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
from group_maker import make_groups


# parser for the arguments
parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('-i', '--input-file',           type=str,                         help='Input file name')
parser.add_argument('-o', '--output-path',          type=str,    default='./',        help='path to save the image')
parser.add_argument('-n', '--n-tps',                type=int,    default=1,           help='number of tps to process')
parser.add_argument('-t', '--ticks_limit',          type=int,    default=3,           help='closeness in ticks to group TPs')
parser.add_argument('-c', '--channel_limit',        type=int,    default=1,           help='closeness in channels to group TPs')
parser.add_argument('-m', '--min-tps',              type=int,    default=3,           help='minimum number of TPs to create a group and then an image')
parser.add_argument('--which-detector',             type=str,    default="APA",       help='"APA", "CRP" or "50L"')
parser.add_argument('--save-groups',                action='store_true',              help='write the groups to a textfile')
parser.add_argument('--show',                       action='store_true',              help='show the image')
# parser.add_argument('--save-ds',                  action='store_true',              help='save the dataset') # discuss if to keep it
parser.add_argument('--fixed-size',                 action='store_true',              help='make the image size fixed')
parser.add_argument('--img-width',                  type=int,    default=40,          help='width of the image')
parser.add_argument('--img-height',                 type=int,    default=300,         help='height of the image')
parser.add_argument('--x-margin',                   type=int,    default=5,           help='margin in x')
parser.add_argument('--y-margin',                   type=int,    default=10,          help='margin in y')
# parser.add_argument('--min-tps-to-create-img', type=int,    default=5,          help='minimum number of TPs to create an image')

args = parser.parse_args()
input_file              = args.input_file
output_path             = args.output_path
show                    = args.show
# save_ds               = args.save_ds
n_tps                   = args.n_tps
ticks_limit             = args.ticks_limit
channel_limit           = args.channel_limit
min_tps                 = args.min_tps
which_detector          = args.which_detector
save_groups             = args.save_groups
fixed_size              = args.fixed_size
width                   = args.img_width
height                  = args.img_height
x_margin                = args.x_margin
y_margin                = args.y_margin
# min_tps_to_create_img = args.min_tps_to_create_img

# recap of selected options
print ("#############################################")
print("Selected options:")
print(" - Input file: " + str(input_file))
print(" - Output path: " + str(output_path))
print(" - Number of TPs: " + str(n_tps))
print(" - Ticks limit: " + str(ticks_limit))
print(" - Channel limit: " + str(channel_limit))
print(" - Minimum number of TPs to create a group: " + str(min_tps))
print(" - Which detector: " + str(which_detector))
print(" - Save groups: " + str(save_groups))
print(" - Show image: " + str(show))
# print(" - Save dataset: " + str(save_ds))
print(" - Fixed size: " + str(fixed_size))
print(" - Image width: " + str(width))
print(" - Image height: " + str(height))
print(" - X margin: " + str(x_margin))
print(" - Y margin: " + str(y_margin))
# print(" - Minimum number of TPs to create an image: " + str(min_tps_to_create_img))
print ("#############################################")

# reading with this function already orders by time_start
all_TPs = save_tps_array(input_file, max_tps=n_tps) 

# create channel map to distinguish induction and collection
channel_map = create_channel_map_array(which_detector=which_detector)

groups = make_groups(all_TPs, channel_map, ticks_limit=ticks_limit, channel_limit=channel_limit, min_tps_to_group=min_tps)
print("Number of groups: ", len(groups))

# Check how many views (U, V, X) are in this sample
total_channels = channel_map.shape[0]
n_views = np.unique(channel_map[all_TPs['channel'] % total_channels, 1]).shape[0]
print("Number of different views in these groups: ", n_views)
print (" ")

# Prepare the output folder
if not os.path.exists(output_path):
    os.makedirs(output_path)
    print ("Created output folder: " + output_path)

print("Creating images...")
for i, group in enumerate(groups):
    save_image(np.array(group), channel_map, output_path=output_path, output_name= "track" + str(i), min_tps_to_create_img=min_tps, make_fixed_size=fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, only_collection=False, show=show)
print("Done!")
print(" ")
    
if save_groups:
    print("Writing groups to text file...")
    with open(output_path + 'groups.txt', 'w') as f:
        for i, group in enumerate(groups):
            f.write('Group' + str(i) + ':\n')
            for tp in group:
                f.write()
                this_tp = np.array([tp["time_start"], tp["time_over_threshold"], tp["time_peak"], tp["channel"], tp["adc_integral"], tp["adc_peak"], tp["detid"], tp["type"], tp["algorithm"], tp["version"], tp["flag"]])
                f.write (str(this_tp) + '\n')
                # f.write(f"{tp["time_start"]} {tp[1]} {tp[2]} {tp[3]} {tp[4]} {tp[5]} {tp[6]} {tp[7]}\n")
                # f.write(f"{tp[0]} {tp[1]} {tp[2]} {tp[3]} {tp[4]} {tp[5]} {tp[6]} {tp[7]}\n")
            f.write('\n')
    print("Done!")
    print(" ")

# if save_ds:
#     print("Creating dataset...")
#     if not os.path.exists(output_path+'dataset/'):
#         os.makedirs(output_path+'dataset/')
#     dataset_img, dataset_label = tp2img.create_dataset(groups,  channel_map=channel_map, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, n_views=n_views)
#     print("Dataset shape: ", dataset_img.shape)
#     print("Dataset label shape: ", dataset_label.shape)

#     # print how many entries are different from 0
#     print("Number of non-zero entries: ", np.count_nonzero(dataset_img))

#     np.save(output_path+'dataset/dataset_img.npy', dataset_img)
#     np.save(output_path+'dataset/dataset_label.npy', dataset_label)

#     print(np.unique(dataset_label,return_counts=True))
#     print("Done!")
# print('Done!')

