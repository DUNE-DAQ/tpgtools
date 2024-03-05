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
# from scipy import sparse

from utils import save_tps_array, create_channel_map_array


def create_image_one_view(tps_to_draw, make_fixed_size=False, width=100, height=100, x_margin=10, y_margin=100, y_min_overall=-1, y_max_overall=-1, verbose=False):
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
        t_start = tps_to_draw[0]['time_start'] 

    if y_max_overall != -1:
        t_end = y_max_overall
    else:
        t_end = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])
    
    x_max = (tps_to_draw['channel'].max())
    x_min = (tps_to_draw['channel'].min())

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
            x = (tp['channel'] - x_min) + x_margin
            y_start = (tp['time_start'] - t_start) + y_margin
            y_end = (y_start + tp['time_over_threshold'])
            
            # print all the values used in the next line 
            if verbose:
                print(f'x: {x}, y_start: {y_start}, y_end: {y_end}, tp["adc_integral"]: {tp["adc_integral"]}, y_end - y_start: {y_end - y_start}')
            
            img[int(y_start)-1:int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
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
                x=(tp['channel'] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp['time_start'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
        elif stretch_x:
            for tp in tps_to_draw:                
                x=(tp['channel'] - x_min)/x_range * (img_width - 2*x_margin) + x_margin
                y_start = (tp['time_start'] - t_start) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
        elif stretch_y:
            for tp in tps_to_draw:
                x = (tp['channel'] - x_min) + x_margin
                y_start = (tp['time_start'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start):int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
        else:
            for tp in tps_to_draw:
                x = (tp['channel'] - x_min) + x_margin
                y_start = (tp['time_start'] - t_start) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
   
    return img

def create_images(tps_to_draw, channel_map, min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200, only_collection=False, verbose=False):
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
    print ("Creating images...")
    y_min_overall = tps_to_draw[0]["time_start"]
    y_max_overall = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])
    total_channels = channel_map.shape[0]
    img_u, img_v, img_x = np.array([[-1]]), np.array([[-1]]), np.array([[-1]])

    # X plane, take only the tps where the corrisponding position in the channel map is 2
    print ("Creating X plane images...")
    tps_x = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 2)]
    if tps_x.shape[0] >= min_tps_to_create_img:
        img_x = create_image_one_view(tps_x, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    if only_collection:
        print (" ")
        return img_u, img_v, img_x # calling here to avoid wasting execution time. U and V will be empty
    
    # U plane, take only the tps where the corrisponding position in the channel map is 0      
    print ("Creating U plane images...")
    tps_u = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 0)]
    if tps_u.shape[0] >= min_tps_to_create_img:
        img_u = create_image_one_view(tps_u, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    
    # V plane, take only the tps where the corrisponding position in the channel map is 1
    print ("Creating V plane images...")
    tps_v = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 1)]
    if tps_v.shape[0] >= min_tps_to_create_img: 
        img_v = create_image_one_view(tps_v, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, y_min_overall=y_min_overall, y_max_overall=y_max_overall)
    
    print (" ")
    return img_u, img_v, img_x

def show_image(tps_to_draw, channel_map, min_tps_to_create_img=2, make_fixed_size=False, width=500, height=1000, x_margin=10, y_margin=200, only_collection=False, img_u=None, img_v=None, img_x=None):
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
   
    # If not given as argument, create images
    if img_u[0, 0] == -1 or img_v[0, 0] == -1 or img_x[0, 0] == -1:
        img_u, img_v, img_x = create_images(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, only_collection=only_collection)

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

def save_image(tps_to_draw, channel_map, output_path, output_name='test', min_tps_to_create_img=2, make_fixed_size=False, width=200, height=300, x_margin=10, y_margin=200, only_collection=False, show=False):
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
    n_x_ticks = 5
    fontsize = 30
    t_shift = tps_to_draw[0]['time_start']
    tps_to_draw['time_start'] = tps_to_draw['time_start'] - t_shift 
    tps_to_draw['time_peak'] = tps_to_draw['time_peak'] - t_shift

    total_channels = channel_map.shape[0]

    #create images
    img_u, img_v, img_x = create_images(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, only_collection=only_collection)
    max_pixel_value_overall = np.max([np.max(img_u), np.max(img_v), np.max(img_x)])
    
    #save images
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    x_max = (tps_to_draw['channel'].max())
    x_min = (tps_to_draw['channel'].min())

    x_range = x_max - x_min 
    t_start = tps_to_draw[0]['time_start']

    t_end = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])
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
        tps_u = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 0)]
        
        x_min_u = (tps_u['channel'].min())
        x_max_u = (tps_u['channel'].max())
        x_range_u = x_max_u - x_min_u
        x_margin_u = x_margin

        if make_fixed_size:
            if img_width > x_range_u:
                x_margin_u = int((img_width - x_range_u)/2)    
        xticks_labels_u = [x_min_u-x_margin_u + i*(x_range_u + 2*x_margin_u)//n_x_ticks for i in range(n_x_ticks)]

        n_views += 1
        plt.figure(figsize=(20, 100))
        plt.title('U plane', fontsize=fontsize)
        plt.imshow(img_u)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_u.shape[0], img_u.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_u.shape[1], img_u.shape[1]/n_x_ticks), labels=xticks_labels_u, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'u_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_v[0, 0] != -1:
        tps_v = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 1)]

        x_min_v = (tps_v['channel'].min())
        x_max_v = (tps_v['channel'].max())
        x_range_v = x_max_v - x_min_v

        x_margin_v = x_margin

        if make_fixed_size:
            if img_width > x_range_v:
                x_margin_v = int((img_width - x_range_v)/2)
        xticks_labels_v = [x_min_v-x_margin_v + i*(x_range_v + 2*x_margin_v)//n_x_ticks for i in range(n_x_ticks)]
        
        n_views += 1
        plt.figure(figsize=(20, 100))    
        plt.title('V plane', fontsize=fontsize)
        plt.imshow(img_v)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)

        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_v.shape[1], img_v.shape[1]/n_x_ticks), labels=xticks_labels_v, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'v_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_x[0, 0] != -1:
        tps_x = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 2)]

        x_min_x = (tps_x['channel'].min())
        x_max_x = (tps_x['channel'].max())
        x_range_x = x_max_x - x_min_x

        x_margin_x = x_margin

        if make_fixed_size:
            if img_width > x_range_x:
                x_margin_x = int((img_width - x_range_x)/2)
        xticks_labels_x = [x_min_x-x_margin_x + i*(x_range_x + 2*x_margin_x)//n_x_ticks for i in range(n_x_ticks)]

        n_views += 1
        plt.figure(figsize=(20, 100))
        plt.title('X plane', fontsize=fontsize)
        plt.imshow(img_x)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_x.shape[0], img_x.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_x.shape[1], img_x.shape[1]/n_x_ticks), labels=xticks_labels_x, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'x_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()
    if n_views == 0:
        print(f'No images saved, if this is unexpected you might want to change some parameters.')
        print(f'You selected {tps_to_draw.shape[0]} tps and {min_tps_to_create_img} as minimum tps to create an image.' )
        print(f'  You have {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 0)].shape[0]} U tps,')
        print(f'  {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 1)].shape[0]} V tps')
        print(f'  and {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 2)].shape[0]} Z tps.')


    if n_views > 1:
        fig = plt.figure(figsize=(60, 40))
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(1,3),
                        axes_pad=0.5,
                        share_all=False,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="10%",
                        cbar_pad=0.25,
                        )   
        
        if img_u[0, 0] != -1:
            im = grid[0].imshow(img_u, vmin=0, vmax=max_pixel_value_overall)
            grid[0].set_title('U plane', fontsize=fontsize)
            grid[0].set_xticks(np.arange(0, img_u.shape[1], img_u.shape[1]/n_x_ticks))
            xticks_labels_u = [x_min_u-x_margin_u + i*(x_range_u + 2*x_margin_u)//n_x_ticks for i in range(n_x_ticks)]
            grid[0].set_xticklabels(xticks_labels_u, fontsize=fontsize)
        if img_v[0, 0] != -1:
            im = grid[1].imshow(img_v, vmin=0, vmax=max_pixel_value_overall)
            grid[1].set_title('V plane', fontsize=fontsize)
            grid[1].set_xticks(np.arange(0, img_v.shape[1], img_v.shape[1]/n_x_ticks))
            xticks_labels_v = [x_min_v-x_margin_v + i*(x_range_v + 2*x_margin_v)//n_x_ticks for i in range(n_x_ticks)]
            grid[1].set_xticklabels(xticks_labels_v, fontsize=fontsize)
        if img_x[0, 0] != -1:
            im = grid[2].imshow(img_x, vmin=0, vmax=max_pixel_value_overall)
            grid[2].set_title('X plane', fontsize=fontsize)
            grid[2].set_xticks(np.arange(0, img_x.shape[1], img_x.shape[1]/n_x_ticks))
            xticks_labels_x = [x_min_x-x_margin_x + i*(x_range_x + 2*x_margin_x)//n_x_ticks for i in range(n_x_ticks)]
            grid[2].set_xticklabels(xticks_labels_x, fontsize=fontsize)

        grid.cbar_axes[0].colorbar(im)
        # grid.cbar.ax.tick_params(labelsize=fontsize) 
        grid.cbar_axes[0].set_yticklabels(np.arange(0, max_pixel_value_overall, max_pixel_value_overall/10, dtype=int), fontsize=fontsize, )
        # grid.axes_llc.set_yticks(yticks_labels)
        # use the same yticks_labels for all the images
        grid.axes_llc.set_yticks(np.arange(0, img_v.shape[0], img_v.shape[0]/10))
        grid.axes_llc.set_yticklabels(yticks_labels, fontsize=fontsize)
        
        # save the image
        plt.savefig(output_path+ 'multiview_' + os.path.basename(output_name) + '.png')
        plt.close()
    
    # just a more compact way of calling this instead of having to repeat the arguments
    if show:
        show_image(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, width=width, height=height, x_margin=x_margin, y_margin=y_margin, only_collection=only_collection, img_u=img_u, img_v=img_v, img_x=img_x)

def create_event_display(tps_to_draw, channel_map, min_tps_to_create_img, output_name='test', make_fixed_size=False, height=300, only_collection=False, verbose=False, y_margin=100):
    '''
    :param tps_to_draw: all trigger primitives in the event
    :param channel_map: channel map
    :param output_name: base name for the images
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param height: height of the image if make_fixed_size is True
    '''
    print ("Creating images...")
    y_min_overall = tps_to_draw[0]["time_start"]
    y_max_overall = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])
    total_channels = channel_map.shape[0]
    img_u, img_v, img_x = np.array([[-1]]), np.array([[-1]]), np.array([[-1]])

    # X plane, take only the tps where the corrisponding position in the channel map is 2
    print ("Creating X plane images...")
    tps_x = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 2)]
    if tps_x.shape[0] >= min_tps_to_create_img:
        img_x = create_event_display_one_view(tps_x, channel_map=channel_map, make_fixed_size=make_fixed_size, height=height, y_min_overall=y_min_overall, y_max_overall=y_max_overall, view=2)
    if only_collection:
        print (" ")
        return img_u, img_v, img_x # calling here to avoid wasting execution time. U and V will be empty
    
    # U plane, take only the tps where the corrisponding position in the channel map is 0      
    print ("Creating U plane images...")
    tps_u = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 0)]
    if tps_u.shape[0] >= min_tps_to_create_img:
        img_u = create_event_display_one_view(tps_u, channel_map=channel_map, make_fixed_size=make_fixed_size, height=height, y_min_overall=y_min_overall, y_max_overall=y_max_overall, view=0)
    
    # V plane, take only the tps where the corrisponding position in the channel map is 1
    print ("Creating V plane images...")
    tps_v = tps_to_draw[np.where(channel_map[tps_to_draw['channel']% total_channels, 1] == 1)]
    if tps_v.shape[0] >= min_tps_to_create_img: 
        img_v = create_event_display_one_view(tps_v, channel_map=channel_map, make_fixed_size=make_fixed_size, height=height, y_min_overall=y_min_overall, y_max_overall=y_max_overall, view=1)
    
    print (" ")
    return img_u, img_v, img_x

def create_event_display_one_view(tps_to_draw, channel_map, make_fixed_size=False, height=1000, y_min_overall=-1, y_max_overall=-1, view=2, y_margin=100):
    '''
    :param tps_to_draw: all trigger primitives to draw
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param height: height of the image
    :param y_min_overall: minimum y value of the image
    :param y_max_overall: maximum y value of the image
    :return: image
    '''

    if y_min_overall != -1:
        t_start = y_min_overall
    else:
        t_start = tps_to_draw[0]['time_start'] 

    if y_max_overall != -1:
        t_end = y_max_overall
    else:
        t_end = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])

    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_height = height
    else:
        img_height =  y_range + 2*y_margin
    
    these_channels = np.where(channel_map[:, 1] == view)
    x_min = these_channels[0].min()
    x_max = these_channels[0].max()
    x_range = x_max - x_min

    img = np.zeros((img_height, x_range))

    # fill image
    if (not make_fixed_size):
        for tp in tps_to_draw:
            x = tp['channel']
            y_start = (tp['time_start'] - t_start) + y_margin
            y_end = (y_start + tp['time_over_threshold'])
            img[int(y_start):int(y_end), int(x)] = tp['adc_integral']/(y_end - y_start)
    else:
    # We stretch the image inwards if needed but we do not upscale it. In this second case we build a padded image
        if img_height < y_range:
            stretch_y = True
            print('Warning: image height is smaller than the range of the TPs. The image will be stretched inwards.')
        else:
            y_margin = (img_height - y_range)/2
            stretch_y = False

        if stretch_y:
            for tp in tps_to_draw:
                x = (tp['channel'] - x_min)
                y_start = (tp['time_start'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start)/y_range * (img_height - 2*y_margin) + y_margin
                img[int(y_start):int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)
        else:
            for tp in tps_to_draw:
                x = (tp['channel'] - x_min)
                y_start = (tp['time_start'] - t_start) + y_margin
                y_end = (tp['time_start'] + tp['time_over_threshold'] - t_start) + y_margin
                img[int(y_start)-1:int(y_end), int(x)-1] = tp['adc_integral']/(y_end - y_start)

 
    return img

def save_single_event_display(tps_to_draw, channel_map, min_tps_to_create_img, output_path, output_name='test', make_fixed_size=False, height=300, only_collection=False, verbose=False, y_margin=100):
    '''
    :param tps_to_draw: all trigger primitives in the event
    :param channel_map: channel map
    :param output_path: path where to save the images
    :param output_name: base name for the images
    :param min_tps: minimum number of TPs to create the image
    :param make_fixed_size: if True, the image will have fixed size, otherwise it will be as big as the TPs
    :param height: height of the image if make_fixed_size is True
    '''

    n_x_ticks = 4

    fontsize = 30
    t_shift = tps_to_draw[0]['time_start']
    tps_to_draw['time_start'] = tps_to_draw['time_start'] - t_shift 
    tps_to_draw['time_peak'] = tps_to_draw['time_peak'] - t_shift

    total_channels = channel_map.shape[0]

    #save images
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    #create images
    img_u, img_v, img_x = create_event_display(tps_to_draw, channel_map, min_tps_to_create_img=min_tps_to_create_img, make_fixed_size=make_fixed_size, height=height, only_collection=only_collection, y_margin=y_margin)
    max_pixel_value_overall = np.max([np.max(img_u), np.max(img_v), np.max(img_x)])

    t_start = tps_to_draw[0]['time_start'] 
    t_end = np.max(tps_to_draw['time_start'] + tps_to_draw['time_over_threshold'])

    y_range = int((t_end - t_start))

    # create the image
    if make_fixed_size:
        img_height = height
        if img_height > y_range:
            y_margin = (img_height - y_range)/2
    else:
        img_height =  y_range + 2*y_margin
    yticks_labels = [t_start-y_margin + i*(y_range + 2*y_margin)//10 for i in range(10)]

    n_views = 0
    xticks_labels_u = []
    xticks_labels_v = []
    xticks_labels_x = []

    if img_u[0, 0] != -1:
        n_views += 1
        max_value = np.where(channel_map[:, 1] == 0)[0].max()
        min_value = np.where(channel_map[:, 1] == 0)[0].min()
        x_range = max_value - min_value
        xticks_labels_u = [min_value + i*x_range//n_x_ticks for i in range(n_x_ticks)]

        plt.figure(figsize=(20, 100))
        plt.title('U plane', fontsize=fontsize)
        plt.imshow(img_u)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_u.shape[0], img_u.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_u.shape[1], img_u.shape[1]/n_x_ticks), labels=xticks_labels_u, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'u_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()
    
    if img_v[0, 0] != -1:
        n_views += 1
        max_value = np.where(channel_map[:, 1] == 1)[0].max()
        min_value = np.where(channel_map[:, 1] == 1)[0].min()
        x_range = max_value - min_value
        xticks_labels_v = [min_value + i*x_range//n_x_ticks for i in range(n_x_ticks)]

        plt.figure(figsize=(20, 100))    
        plt.title('V plane', fontsize=fontsize)
        plt.imshow(img_v)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)

        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_v.shape[0], img_v.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_v.shape[1], img_v.shape[1]/n_x_ticks), labels=xticks_labels_v, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'v_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if img_x[0, 0] != -1:
        n_views += 1

        max_value = np.where(channel_map[:, 1] == 2)[0].max()
        min_value = np.where(channel_map[:, 1] == 2)[0].min()
        x_range = max_value - min_value
        xticks_labels_x = [min_value + i*x_range//n_x_ticks for i in range(n_x_ticks)]

        plt.figure(figsize=(20, 100))
        plt.title('X plane', fontsize=fontsize)
        plt.imshow(img_x)
        # plt.colorbar()
        # add x and y labels
        plt.xlabel('channel', fontsize=fontsize)
        plt.ylabel("Time (ticks)", fontsize=fontsize)
        # set y axis ticks
        plt.yticks(ticks=np.arange(0, img_x.shape[0], img_x.shape[0]/10), labels=yticks_labels, fontsize=fontsize)
        # set x axis ticks
        plt.xticks(ticks=np.arange(0, img_x.shape[1], img_x.shape[1]/n_x_ticks), labels=xticks_labels_x, fontsize=fontsize)

        # save the image, with a bbox in inches smaller than the default but bigger than tight
        plt.savefig(output_path+ 'x_' + os.path.basename(output_name) + '.png', bbox_inches='tight', pad_inches=1)
        plt.close()

    if n_views == 0:
        print(f'No images saved, if this is unexpected you might want to change some parameters.')
        print(f'You selected {tps_to_draw.shape[0]} tps and {min_tps_to_create_img} as minimum tps to create an image.' )
        print(f'  You have {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 0)].shape[0]} U tps,')
        print(f'  {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 1)].shape[0]} V tps')
        print(f'  and {tps_to_draw[np.where(channel_map[tps_to_draw["channel"]% total_channels, 1] == 2)].shape[0]} Z tps.')


    if n_views > 1:
        fig = plt.figure(figsize=(60, 40))
        grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                        nrows_ncols=(1,3),
                        axes_pad=0.5,
                        share_all=False,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="10%",
                        cbar_pad=0.25,
                        )   
        
        if img_u[0, 0] != -1:
            im = grid[0].imshow(img_u, vmin=0, vmax=max_pixel_value_overall)
            grid[0].set_title('U plane', fontsize=fontsize)
            grid[0].set_xticks(np.arange(0, img_u.shape[1], img_u.shape[1]/n_x_ticks))
            grid[0].set_xticklabels(xticks_labels_u, fontsize=fontsize)
        if img_v[0, 0] != -1:
            im = grid[1].imshow(img_v, vmin=0, vmax=max_pixel_value_overall)
            grid[1].set_title('V plane', fontsize=fontsize)
            grid[1].set_xticks(np.arange(0, img_v.shape[1], img_v.shape[1]/n_x_ticks))
            grid[1].set_xticklabels(xticks_labels_v, fontsize=fontsize)
        if img_x[0, 0] != -1:
            im = grid[2].imshow(img_x, vmin=0, vmax=max_pixel_value_overall)
            grid[2].set_title('X plane', fontsize=fontsize)
            grid[2].set_xticks(np.arange(0, img_x.shape[1], img_x.shape[1]/n_x_ticks))
            grid[2].set_xticklabels(xticks_labels_x, fontsize=fontsize)
            
        grid.cbar_axes[0].colorbar(im)
        # grid.cbar.ax.tick_params(labelsize=fontsize) 
        grid.cbar_axes[0].set_yticklabels(np.arange(0, max_pixel_value_overall, max_pixel_value_overall/10, dtype=int), fontsize=fontsize, )
        # grid.axes_llc.set_yticks(yticks_labels)
        # use the same yticks_labels for all the images
        grid.axes_llc.set_yticks(np.arange(0, img_v.shape[0], img_v.shape[0]/10))
        grid.axes_llc.set_yticklabels(yticks_labels, fontsize=fontsize)
        
        # save the image
        plt.savefig(output_path+ 'multiview_' + os.path.basename(output_name) + '.png')
        plt.close()

def save_event_display(tps_to_draw, channel_map, min_tps_to_create_img, output_path, output_name='test', make_fixed_size=False, height=300, only_collection=False, verbose=False, y_margin=100):

    n_channels_per_apa = channel_map.shape[0]
    apas = np.unique(tps_to_draw['channel']//n_channels_per_apa)    
    for apa in apas:
        tps_apa = tps_to_draw[np.where(tps_to_draw['channel']//n_channels_per_apa == apa)]
        tps_apa['channel'] = tps_apa['channel']%n_channels_per_apa
        save_single_event_display(tps_apa, channel_map, min_tps_to_create_img, output_path, output_name=f'{output_name}_apa{apa}', make_fixed_size=make_fixed_size, height=height, only_collection=only_collection, verbose=verbose)

