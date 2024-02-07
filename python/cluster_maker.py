"""
 * @file cluster_maker.py
 * @brief Reads takes tps array and makes clusters.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import sys
# import time
import os

from utils import save_tps_array, create_channel_map_array

# import detchannelmaps

# This function creates the clusters basing on channel and time proximity
def make_clusters(all_tps, channel_map, ticks_limit=3, channel_limit=1, min_tps_to_cluster=2, cluster_induction_view=False):
    '''
    :param all_tps: all trigger primitives in the event
    :param channel_map: channel map
    :param ticks_limit: maximum time window to consider, in ticks
    :param channel_limit: maximum channel distance to consider in the same cluster
    :param min_tps_to_cluster: minimum number of hits to consider to form a cluster
    :return: list of clusters
    '''
    # If I have tps from different planes, clustering only basing on time.
    # To avoid this, you can remove all the tps from different planes from the all_tps array
    total_channels = channel_map.shape[0] # 2560 for APA, 3072 for CRP, 128 for 50L
    if np.unique(channel_map[all_tps["channel"] % total_channels, 1]).shape[0] != 1:
        if cluster_induction_view == False:
            # exclude all induction TPs from the all_tps array
            all_tps = all_tps[np.where(channel_map[all_tps["channel"] % total_channels, 1] != 2)]            
        else:
            print('Warning: induction plane included in the list of TPs. clustering only by time.')
            return cluster_maker_only_by_time(all_tps, channel_map, ticks_limit=ticks_limit, channel_limit=channel_limit, min_tps_to_cluster=min_tps_to_cluster)

    # When only collection plane is present, we can cluster by time and channel proximity
    clusters = []
    buffer = []
    # Loop over the TPs and create clusters.
    # The idea is that TPs are iterated over the time and are place in a buffer.
    # If there is a TP close in terms of channel, it is added to the buffer.
    # When there are no more close TPs in terms of time (valid since they should come ordered), 
    # the cluster is added to the list of clusters or discarded if there are not enough TPs in it.
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
                    if len(candidate) >= min_tps_to_cluster:
                        clusters.append(np.array(candidate))
            if not appended:
                buffer.append([tp])
    if len(buffer) > 0:
        for candidate in buffer:
            if len(candidate) >= min_tps_to_cluster:
                clusters.append(np.array(candidate))

    return clusters

# This function creates clusters basing only on time proximity, used when including induction view
def make_clusters_only_by_time(all_tps, channel_map, ticks_limit=5, min_tps_to_cluster=4):
    clusters = []
    current_cluster = [all_tps[0]]
    # print names of dict
    print ("all_tps.dtype.names")
    print(all_tps.dtype.names)
    for tp in all_tps[1:]:
        if tp["time_start"] - np.max(current_cluster["time_start"]+current_cluster["time_over_threshold"]) <= ticks_limit:
            current_cluster = np.append(current_cluster, tp)
        else:
            if len(current_cluster) >= min_tps_to_cluster:
                clusters.append(current_cluster)
            current_cluster = [tp]
 
    return clusters
