"""
 * @file image_creator.py
 * @brief Reads array of TPs and converts them to images.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import sys
# import time
import os

import matplotlib.pyplot as plt
# check if ok in daq framework
from mpl_toolkits.axes_grid1 import ImageGrid 
from scipy import sparse


def create_image_one_view(tps_to_draw, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=100, y_min_overall=-1, y_max_overall=-1):
    '''
    :param tps_to_draw: all trigger primitives to draw
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image
    :param height: height of the image
    :param x_margin: margin on the x axis
    :param y_margin: margin on the y axis
    :param y_min_overall: minimum y value of the image
    :param y_max_overall: maximum y value of the image
    :return: image
    '''

    if y_min_overall != -1:
        t_start = y_min_overall
    else:
        t_start = tps_to_draw[0]["time_start"] 

    if y_max_overall != -1:
        t_end = y_max_overall
    else:
        t_end = np.max(tps_to_draw["time_start"] + tps_to_draw["time_over_threshold"])
    
    x_max = (tps_to_draw["channel"].max())
    x_min = (tps_to_draw["channel"].min())

    x_range = x_max - x_min 
    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_width = width
        img_height = height
    else:
        img_width =  x_range + 2*x_margin
        img_height =  y_range + 2*y_margin
    img = np.zeros((img_height, img_width))
    
    # fill image
    if (not make_fixed_size):
        for tp in tps_to_draw:
            x = (tp["channel"] - x_min) + x_margin
            y_start = (tp["time_start"] - t_start) + y_margin
            y_end = (tp["time_start"] + tp["time_over_threshold"] - t_start) + y_margin
            img[int(y_start)-1:int(y_end), int(x)-1] = tp["adc_integral"]/(y_end - y_start)
    else:
    # We stretch the image inwards if needed but we do not upscale it. In this second case we build a padded image
        if img_width < x_range:
            stretch_x = True
            print('Warning: image width is smaller than the range of the TPs. The image will be stretched inwards.')
        else:
            x_margin = (img_width - x_range)/2
            stretch_x = False
            
        if img_height < y_range:
            stretch_y = True
            print('Warning: image height is smaller than the range of the TPs. The image will be stretched inwards.')
        else:
            y_margin = (img_height - y_range)/2
            stretch_y = False

        if stretch_x & stretch_y:
            for tp in tps_to_draw:
                x=(tp["channel"] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp["time_start"] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp["time_start"] + tp["time_over_threshold"] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp["adc_integral"]/(y_end - y_start)
        elif stretch_x:
            for tp in tps_to_draw:                
                x=(tp["channel"] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp["time_start"] - t_start) + y_margin
                y_end = (tp["time_start"] + tp["time_over_threshold"] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp["adc_integral"]/(y_end - y_start)
        elif stretch_y:
            for tp in tps_to_draw:
                x = (tp["channel"] - x_min) + x_margin
                y_start = (tp["time_start"] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp["time_start"] + tp["time_over_threshold"] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start):int(y_end), int(x)-1] = tp["adc_integral"]/(y_end - y_start)
        else:
            for tp in tps_to_draw:
                x = (tp["channel"] - x_min) + x_margin
                y_start = (tp["time_start"] - t_start) + y_margin
                y_end = (tp["time_start"] + tp["time_over_threshold"] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp["adc_integral"]/(y_end - y_start)
   
    return img

def create_images(tps_to_draw, channel_map, min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200, only_collection=True):
    '''
    :param tps_to_draw: all trigger primitives to draw
    :param channel_map: channel map array
    :param min_tps_to_create_img: minimum number of TPs to create an image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image
    :param height: height of the image
    :param x_margin: margin on the x axis
    :param y_margin: margin on the y axis
    :return: images or -1 if there are not enough TPs
    '''

    y_min_overall = tps_to_draw[0]["time_start"]
    y_max_overall = np.max(tps_to_draw["time_start"] + tps_to_draw["time_over_threshold"])
    total_channels = channel_map.shape[0]
    img_u, img_v, img_x = np.array([[-1]]), np.array([[-1]]), np.array([[-1]])

    # X plane, take only the tps where the corrisponding position in the channel map is 2
    tps_x = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 2)]
    if tps_x.shape[0] >= min_tps_to_create_img:
        img_x = create_image_one_view(tps_x, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    if only_collection:
        return img_x # calling here to avoid wasting execution time
    
    # U plane, take only the tps where the corrisponding position in the channel map is 0      
    tps_u = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 0)]
    if tps_u.shape[0] >= min_tps_to_create_img:
        img_u = create_image_one_view(tps_u, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    
    # V plane, take only the tps where the corrisponding position in the channel map is 1
    tps_v = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 1)]
    if tps_v.shape[0] >= min_tps_to_create_img: 
        img_v = create_image_one_view(tps_v, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    
    return img_u, img_v, img_x

def show_image(tps_to_draw, channel_map, min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200, img_u, img_v, img_x):
    '''
    :param tps_to_draw: all trigger primitives in the group
    :param channel_map: channel map
    :param min_tps: minimum number of TPs to create the image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image if make_fixed_size is True
    :param height: height of the image if make_fixed_size is True
    :param x_margin: margin on the x axis if make_fixed_size is False
    :param y_margin: margin on the y axis if make_fixed_size is False
    '''
    img_u = []
    img_v = []
    img_x = []
    
    # If not given as argument, create images
    if img_u[0, 0] == -1 or img_v[0, 0] == -1 or img_x[0, 0] == -1:
        img_u, img_v, img_x = create_images(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)

    n_views = 0
    #show images
    if img_u[0, 0] != -1:
        n_views += 1
        plt.figure(figsize=(8, 20))
        plt.imshow(img_u)
        plt.show()
    if img_v[0, 0] != -1:
        n_views += 1
        plt.figure(figsize=(8, 20))
        plt.imshow(img_v)
        plt.show()
    if img_x[0, 0] != -1:
        n_views += 1
        plt.figure(figsize=(8, 20))
        plt.imshow(img_x)
        plt.show()

    if n_views > 1:
        fig = plt.figure()
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(1,3),
                        axes_pad=0.15,
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="7%",
                        cbar_pad=0.15,
                        )
        if img_u[0, 0] != -1:
            im = grid[0].imshow(img_u)
            grid[0].set_title('U plane')
        if img_v[0, 0] != -1:
            im = grid[1].imshow(img_v)
            grid[1].set_title('V plane')
        if img_x[0, 0] != -1:
            im = grid[2].imshow(img_x)
            grid[2].set_title('X plane')
        grid.cbar_axes[0].colorbar(im)
        grid.axes_llc.set_yticks(np.arange(0, img_u.shape[0], 100))
        plt.show()

def save_image(tps_to_draw, channel_map, output_path, output_name='test', min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200, show_image=False):
    '''
    :param tps_to_draw: all trigger primitives in the event
    :param channel_map: channel map
    :param output_path: path where to save the images
    :param output_name: base name for the images
    :param min_tps: minimum number of TPs to create the image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param width: width of the image if make_fixed_size is True
    :param height: height of the image if make_fixed_size is True
    :param x_margin: margin on the x axis if make_fixed_size is False
    :param y_margin: margin on the y axis if make_fixed_size is False
    '''
    total_channels = channel_map.shape[0]

    #create images
    img_u, img_v, img_x = create_images(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
    max_pixel_value_overall = np.max([np.max(img_u), np.max(img_v), np.max(img_x)])

    #save images
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    x_max = (tps_to_draw["channel"].max())
    x_min = (tps_to_draw["channel"].min())

    x_range = x_max - x_min 
    t_start = tps_to_draw[0]["time_start"]

    t_end = np.max(tps_to_draw["time_start"] + tps_to_draw["time_over_threshold"])
    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_width = width
        img_height = height
        if img_width > x_range:
            x_margin = (img_width - x_range)/2    
        if img_height > y_range:
            y_margin = (img_height - y_range)/2
    else:
        img_width =  x_range + 2*x_margin
        img_height =  y_range + 2*y_margin
    yticks_labels = [t_start-y_margin + i*(y_range + 2*y_margin)//10 for i in range(10)]

    n_views = 0

    if img_u[0, 0] != -1:
        tps_u = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 0)]
        
        x_min_u = (tps_u["channel"].min())
        x_max_u = (tps_u["channel"].max())
        x_range_u = x_max_u - x_min_u
        x_margin_u = x_margin

        if make_fixed_size:
            if img_width > x_range_u:
                x_margin_u = (img_width - x_range_u)/2    
        xticks_labels_u = [x_min_u-x_margin_u + i*(x_range_u + 2*x_margin_u)//2 for i in range(2)]

        n_views += 1
        plt.figure(figsize=(10, 26))
        plt.title('U plane')
        plt.imshow(img_u)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_u.shape[0], img_u.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_u.shape[1], img_u.shape[1]/2), labels=xticks_labels_u)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'u_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_v[0, 0] != -1:
        tps_v = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 1)]

        x_min_v = (tps_v["channel"].min())
        x_max_v = (tps_v["channel"].max())
        x_range_v = x_max_v - x_min_v

        x_margin_v = x_margin

        if make_fixed_size:
            if img_width > x_range_v:
                x_margin_v = (img_width - x_range_v)/2
        xticks_labels_v = [x_min_v-x_margin_v + i*(x_range_v + 2*x_margin_v)//2 for i in range(2)]
        
        n_views += 1
        plt.figure(figsize=(10, 26))    
        plt.title('V plane')
        plt.imshow(img_v)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")

        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_v.shape[1], img_v.shape[1]/2), labels=xticks_labels_v)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'v_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_x[0, 0] != -1:
        tps_x = tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 2)]

        x_min_x = (tps_x["channel"].min())
        x_max_x = (tps_x["channel"].max())
        x_range_x = x_max_x - x_min_x

        x_margin_x = x_margin

        if make_fixed_size:
            if img_width > x_range_x:
                x_margin_x = (img_width - x_range_x)/2
        xticks_labels_x = [x_min_x-x_margin_x + i*(x_range_x + 2*x_margin_x)//2 for i in range(2)]

        n_views += 1
        plt.figure(figsize=(10, 26))
        plt.title('X plane')
        plt.imshow(img_x)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel("Channel")
        plt.ylabel("Time (ticks)")
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_x.shape[0], img_x.shape[0]/10), labels=yticks_labels)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_x.shape[1], img_x.shape[1]/2), labels=xticks_labels_x)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'x_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()
    if n_views == 0:
        print(f'No images saved! Do a better grouping algorithm! You had {tps_to_draw.shape[0]} tps and {min_tps_to_create_img} as min_tps_to_create_img.' )
        print(f'You have {tps_to_draw[np.where(channel_map[tps_to_draw[:, 3]% total_channels, 1] == 0)].shape[0]} U tps, {tps_to_draw[np.where(channel_map[tps_to_draw[:, 3]% total_channels, 1] == 1)].shape[0]} V tps and {tps_to_draw[np.where(channel_map[tps_to_draw[:, 3]% total_channels, 1] == 2)].shape[0]} Z tps.')


    if n_views > 1:
        fig = plt.figure(figsize=(10, 26))
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(1,3),
                        axes_pad=0.5,
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="30%",
                        cbar_pad=0.25,
                        )   


        if img_u[0, 0] != -1:
            im = grid[0].imshow(img_u, vmin=0, vmax=max_pixel_value_overall)
            grid[0].set_title('U plane')
            # grid[0].set_xticks(np.arange(0, img_u.shape[1], img_u.shape[1]/2))
            # grid[0].set_xticklabels(xticks_labels_u)
        if img_v[0, 0] != -1:
            im = grid[1].imshow(img_v, vmin=0, vmax=max_pixel_value_overall)
            grid[1].set_title('V plane')
            # grid[1].set_xticks(np.arange(0, img_v.shape[1], img_v.shape[1]/2))
            # grid[1].set_xticklabels(xticks_labels_v)
        if img_x[0, 0] != -1:
            im = grid[2].imshow(img_x, vmin=0, vmax=max_pixel_value_overall)
            grid[2].set_title('X plane')
            # grid[2].set_xticks(np.arange(0, img_x.shape[1], img_x.shape[1]/2))
            # grid[2].set_xticklabels(xticks_labels_x)

        grid.cbar_axes[0].colorbar(im)
        # grid.axes_llc.set_yticks(yticks_labels)
        # use the same yticks_labels for all the images
        grid.axes_llc.set_yticks(np.arange(0, img_v.shape[0], img_v.shape[0]/10))
        grid.axes_llc.set_yticklabels(yticks_labels)
        
        # save the image
        plt.savefig(output_path+ 'multiview_' + os.path.basename(output_name) + '.png')
        plt.close()
    
    # just a more compact way of calling this instead of having to repeat the arguments
    if show_image:
        show_image(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, img_u=img_u, img_v=img_v, img_x=img_x)

# def create_dataset(groups, channel_map, make_fixed_size=True, width=70, height=1000, x_margin=5, y_margin=50, n_views=3, use_sparse=False, unknown_label=99, idx=7, dict_lab=None):
#     '''
#     :param groups: list of groups
#     :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
#     :param width: width of the image
#     :param height: height of the image
#     :param x_margin: margin on the x axis
#     :param y_margin: margin on the y axis
#     :return: dataset [[img],[label]] in numpy array format
#     '''
#     if not use_sparse:
#         # Each pixel must have a value between 0 and 255. Maybe problematic for high ADC values. I can't afford havier data types for the moment
#         dataset_img = np.zeros((len(groups), height, width, n_views), dtype=np.uint8) 
#         dataset_label = np.empty((len(groups), 1), dtype=np.uint8)
#         i=0
#         for group in (groups):

#             # create the label. I have to do it this way because the label is not the same for all the datasets
#             label = label_generator_snana(group, unknown_label=unknown_label, idx=idx, dict_lab=dict_lab)
#             # append to the dataset as an array of arrays
#             if n_views > 1:
#                 img_u, img_v, img_x = all_views_img_maker(np.array(group), channel_map, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
#                 if img_u[0, 0] != -1:
#                     dataset_img[i, :, :, 0] = img_u
#                 if img_v[0, 0] != -1:
#                     dataset_img[i, :, :, 1] = img_v
#                 if img_x[0, 0] != -1:
#                     dataset_img[i, :, :, 2] = img_x 
#             else:  
#                 img = from_tp_to_imgs(np.array(group), make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
#                 if img[0, 0] != -1:
#                     dataset_img[i, :, :, 0] = img

#             dataset_label[i] = [label]            
#             i+=1
#     else:
#         dataset_img = []
#         dataset_label = np.empty((len(groups), 1), dtype=np.uint8)
#         i=0
#         for group in (groups):

#             # create the label. I have to do it this way because the label is not the same for all the datasets
#             label = label_generator_snana(group, unknown_label=unknown_label, idx=idx, dict_lab=dict_lab)
#             # append to the dataset as an array of arrays
#             if n_views > 1:
#                 img_u, img_v, img_x = all_views_img_maker(np.array(group), channel_map, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
#                 if img_u[0, 0] != -1:
#                     img_u = sparse.csr_matrix(img_u)
#                 if img_v[0, 0] != -1:
#                     img_v = sparse.csr_matrix(img_v)
#                 if img_x[0, 0] != -1:
#                     img_x = sparse.csr_matrix(img_x)
#                 dataset_img.append([img_u, img_v, img_x])
#             else:  
#                 img = from_tp_to_imgs(np.array(group), make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin)
#                 if img[0, 0] != -1:
#                     img = sparse.csr_matrix(img)
#                 dataset_img.append([img.data, img.indices, img.indptr])
#             dataset_label[i] = [label]            
#             i+=1
        
#         dataset_img = np.array(dataset_img, dtype=object)

#     return (dataset_img, dataset_label)

# def label_generator_snana(group,idx=7, unknown_label=10, dict_lab=None):

#     label = group[0]["truth"]
#     if np.all(group["truth"] == label):
#         if dict_lab is None:
#             return label
#         else:
#             return dict_lab[label]
#     else:
#         return unknown_label
