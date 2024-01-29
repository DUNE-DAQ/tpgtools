"""
 * @file group_maker.py
 * @brief Reads takes tps array and makes groups.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import sys
# import time
import os

# import detchannelmaps

# This function creates a channel map array, 
# where the first column is the channel number and the second column 
# is the plane view (0 for U, 1 for V, 2 for X).
# TODO understand if this can be avoided by reading from detchannelmaps
def create_channel_map_array(which_detector="APA"):
    '''
    :param which_detector, options are: 
     - APA for Horizontal Drift 
     - CRP for Vertical Drift (valid also for coldbox)
     - 50L for 50L detector, VD technology
    :return: channel map array
    '''
    # create the channel map array
    if which_detector == "APA":
        channel_map = np.empty((2560, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 2559, 2560, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(800), np.ones(800), np.ones(960)*2))        
    elif which_detector == "CRP":
        channel_map = np.empty((3072, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 3071, 3072, dtype=int)
        channel_map[:, 1] = np.concatenate((np.zeros(952), np.ones(952), np.ones(1168)*2))
    elif which_detector == "50L":
        channel_map = np.empty((128, 2), dtype=int)
        channel_map[:, 0] = np.linspace(0, 127, 128, dtype=int)
        # channel_map[:, 1] = np.concatenate(np.ones(49)*2, np.zeros(39), np.ones(40))
        channel_map[:, 1] = np.concatenate((np.zeros(39), np.ones(40), np.ones(49)*2))
    return channel_map

# This function creates the groups basing on channel and time proximity
def make_groups(all_tps, channel_map, ticks_limit=3, channel_limit=1, min_tps_to_group=2, group_induction_view=False):
    '''
    :param all_tps: all trigger primitives in the event
    :param channel_map: channel map
    :param ticks_limit: maximum time window to consider, in ticks
    :param channel_limit: maximum channel distance to consider in the same group
    :param min_tps_to_group: minimum number of hits to consider to form a group
    :return: list of groups
    '''
    # If I have tps from different planes, grouping only basing on time.
    # To avoid this, you can remove all the tps from different planes from the all_tps array
    total_channels = channel_map.shape[0] # 2560 for APA, 3072 for CRP, 128 for 50L
    if np.unique(channel_map[all_tps["channel"] % total_channels, 1]).shape[0] != 1:
        if group_induction_view == False:
            # exclude all induction TPs from the all_tps array
            all_tps = all_tps[np.where(channel_map[all_tps["channel"] % total_channels, 1] != 2)]            
        else:
            print('Warning: induction plane included in the list of TPs. Grouping only by time.')
            return group_maker_only_by_time(all_tps, channel_map, ticks_limit=ticks_limit, channel_limit=channel_limit, min_tps_to_group=min_tps_to_group)

    # When only collection plane is present, we can group by time and channel proximity
    groups = []
    buffer = []
    # Loop over the TPs and create groups.
    # The idea is that TPs are iterated over the time and are place in a buffer.
    # If there is a TP close in terms of channel, it is added to the buffer.
    # When there are no more close TPs in terms of time (valid since they should come ordered), 
    # the group is added to the list of groups or discarded if there are not enough TPs in it.
    for tp in all_tps:
        if len(buffer)==0:
            buffer.append([tp])
        else:
            buffer_copy = buffer.copy()
            buffer = []
            appended = False
            idx = 0
            for candidate in (buffer_copy):
                time_cond = (tp["time_start"] - np.max(np.array(candidate)["time_start"]+np.array(candidate)["time_over_threshold"])) <= ticks_limit
                if time_cond:
                    chan_cond = np.min(np.abs(tp["channel"] - np.array(candidate)["channel"])) <= channel_limit
                    if chan_cond:
                        if not appended:
                            candidate.append(tp)
                            buffer.append(candidate)
                            idx_appended = idx
                            appended = True
                        else:
                            for tp2 in candidate:
                                buffer[idx_appended].append(tp2)
                    else:
                        buffer.append(candidate)
                        idx += 1
                else:
                    if len(candidate) >= min_tps_to_group:
                        groups.append(np.array(candidate, dtype=int))
            if not appended:
                buffer.append([tp])
    if len(buffer) > 0:
        for candidate in buffer:
            if len(candidate) >= min_tps_to_group:
                groups.append(np.array(candidate, dtype=int))

    groups = np.array(groups, dtype=object)
    return groups

# This function creates groups basing only on time proximity, used when including induction view
def make_groups_only_by_time(all_tps, channel_map, ticks_limit=5, min_tps_to_group=4):
    groups = []
    current_group = np.array([all_tps[0]])
    for tp in all_tps[1:]:
        if tp["time_start"] - np.max(current_group["time_start"]+current_group["time_over_threshold"]) <= ticks_limit:
            current_group = np.append(current_group, tp)
        else:
            if len(current_group) >= min_tps_to_group:
                groups.append(current_group)
            current_group = np.array([tp])

    return groups
