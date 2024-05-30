import numpy as np
import pandas as pd
import hdf5libs
import justintime.utils.rawdataunpacker as rdu
from scipy.spatial import ConvexHull
import time
import warnings
# import logging

# class UnpackerServiceFilter(rdu.UnpackerService):
#     """
#     Extends the standard justintime unpacker service to allow a filter
#     on time channel before storing a file
    
#     New method:
#     unpack_filtered(raw_data_file, tr_id, time_min, time_max, seq_id=0)
#     """
#     def unpack_filtered(
#             self,
#             raw_data_file, tr_id: int,
#             time_min: int, time_max: int,
#             seq_id: int=0) -> dict:
#         res = {}

#         trh = raw_data_file.get_trh((tr_id, seq_id))
#         tr_source_ids = raw_data_file.get_source_ids((tr_id, seq_id))

#         for sid in tr_source_ids:
#             frag = raw_data_file.get_frag((tr_id, seq_id),sid)
            
#             if not (frag.get_window_begin() <= time_max
#                     and frag.get_window_end() >= time_min):
#                 continue

#             for n,up in self.fragment_unpackers.items():
#                 if not up.match(frag.get_fragment_type(), sid.subsystem):
#                     # logging.debug(f"fragment {sid} (type {frag.get_fragment_type()}) and unpacker {n} - no match")
#                     continue
                
#                 logging.debug(f"[{n}] Unpacking Subsys={sid.subsystem}, id={sid.id}")                
#                 r = up.unpack(frag)
#                 logging.debug(f"[{n}] Unpacking Subsys={sid.subsystem}, id={sid.id} completed ({len(r) if r is not None else 0})")
#                 res.setdefault(n,{})[sid.id] = r

#         return res


def read_hdf5_time_slices_jit(file_path, init_slice=0, n_slices=1, detector="VDColdbox", verbosity=2):
    """
    Create a pandas dataframe from Time Slices using justintime.

    Parameters
    ----------
    file_path : str
        Path to hdf5 file to read.
    init_slice : int, optional
        Index of the first Time Slice to start reading from. Default
        is 0.
    n_slices : int, optional
        Number of Time Slices to read. If -1, read out all Time Slices
        after init_slice. Default is 1.
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "VDColdbox".
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 1.
    
    Returns
    -------
    pd.Dataframe
        Datafram comtaining TPs in the hdf5 file.
        Columns are:
        - time_start
        - time_peak
        - time_over_threshold
        - channel
        - adc_integral
        - adc_peak
        - flag
        - plane
    """
    tp_up = rdu.TPFragmentPandasUnpacker(detector)
    up = rdu.UnpackerService()
    up.add_unpacker('tp', tp_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    tss = [ i for i,_ in rdf.get_all_timeslice_ids()]
    # process only the first TR
    slice_end = init_slice + n_slices if n_slices != -1 else -1
    tss = tss[init_slice:slice_end]

    df_tps = []
    if verbosity >= 4:
        t_start = time.time()
    for ts in tss:
        if verbosity >= 2:
            print(f"--- Reading Time Slice {ts} ---")

        unpacked_ts = up.unpack(rdf, ts)

        if 'tp' in unpacked_ts:
            if verbosity >= 3:
                print("Assembling TPs")
            df_tp = pd.concat(unpacked_ts['tp'].values())
            df_tps.append(df_tp.sort_values(by=['time_start', 'channel']))
            if verbosity >= 3:
                print(f"TPs dataframe assembled: {len(df_tp)}")
    df_tps = pd.concat(df_tps, ignore_index=True)
    df_tps = df_tps.reindex(sorted(df_tps.columns), axis=1)
    if verbosity >= 1:
        print(f"TPC adcs dataframe assembled {len(df_tps)}x{len(df_tps.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load TPs: {time.time() - t_start:.2f}s")
    return df_tps

def read_hdf5_record_tps_jit(file_path, init_record=0, n_records=1, detector="VDColdbox", verbosity=2):
    """
    Create a pandas dataframe from Records using justintime.

    Parameters
    ----------
    file_path : str
        Path to hdf5 file to read.
    init_record : int, optional
        Index of the first record to start reading from. Default is 0.
    n_records : int, optional
        Number of records to read. If -1, read out all records after
        init_record. Default is 1.
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "VDColdbox".
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 1.
    
    Returns
    -------
    pd.Dataframe
        Datafram comtaining TPs in the hdf5 file.
        Columns are:
        - time_start
        - time_peak
        - time_over_threshold
        - channel
        - adc_integral
        - adc_peak
        - flag
        - plane
    """
    tp_up = rdu.TPFragmentPandasUnpacker(detector)
    up = rdu.UnpackerService()
    up.add_unpacker('tp', tp_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    trs = [ i for i,_ in rdf.get_all_record_ids()]
    # process only the first TR
    record_end = init_record + n_records if n_records != -1 else -1
    trs = trs[init_record:record_end]

    df_tps = []
    if verbosity >= 4:
        t_start = time.time()
    for tr in trs:
        if verbosity >= 2:
            print(f"--- Reading Record {tr} ---")

        unpacked_tr = up.unpack(rdf, tr)

        if 'tp' in unpacked_tr:
            if verbosity >= 3:
                print("Assembling TPs")
            df_tp = pd.concat(unpacked_tr['tp'].values())
            # Record the record of the event for easier access in future
            df_tp = df_tp.assign(record_id=tr)
            df_tps.append(df_tp.sort_values(by=['time_start', 'channel']))
            if verbosity >= 3:
                print(f"TPs dataframe assembled: {len(df_tp)}")
    df_tps = pd.concat(df_tps, ignore_index=True)
    df_tps = df_tps.reindex(sorted(df_tps.columns), axis=1)
    if verbosity >= 1:
        print(f"TPC adcs dataframe assembled {len(df_tps)}x{len(df_tps.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load TPs: {time.time() - t_start:.2f}s")
    return df_tps

def read_hdf5_raw_jit(
        file_path,
        records_list=None, init_record=0, n_records=1,
        detector="VDColdbox",
        verbosity=2):
    """
    Create a pandas dataframe contain raw data.

    Parameters
    ----------
    file_path : str
        Path to hdf5 file to read.
    init_record : int, optional
        Index of the first Record to start reading from. Default is 0.
    n_records : int, optional
        Number of Recordss to read. If -1, read out all Records after
        init_record. Default is 1.
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "VDColdbox".
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 1.
    
    Returns
    -------
    pd.Dataframe
        Dataframe with channels as columns and timestamps as the index
        indicating the recorded ADC of the channel at that time.
    """
    wethf_up = rdu.WIBEthFragmentPandasUnpacker(detector)
    up = rdu.UnpackerService()
    up.add_unpacker('bde_eth', wethf_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    if records_list is None:
        records_list = [ i for i,_ in rdf.get_all_record_ids()]
        # process only the first TR
        record_end = init_record + n_records if n_records != -1 else -1
        records_list = records_list[init_record:record_end]

    df_tpcs = []
    if verbosity >= 4:
        t_start = time.time()
    for tr in records_list:
        if verbosity >= 2:
            print(f"--- Reading Record {tr} ---")

        unpacked_tr = up.unpack(rdf, tr)
        if 'bde_eth' in unpacked_tr:
            dfs_bde = {k:v for k,v in unpacked_tr['bde_eth'].items() if not v is None}
            if verbosity >= 3:
                print(f"Assembling WIBEth Frames {len(dfs_bde)}")

            idx = pd.Index([], dtype='uint64')
            for df in dfs_bde.values():
                idx = idx.union(df.index)

            df_tpc = pd.DataFrame(index=idx, dtype='uint16')

            for df in dfs_bde.values():
                df_tpc = df_tpc.join(df)

            df_tpcs.append(df_tpc.reindex(sorted(df_tpc.columns), axis=1))

            if verbosity >= 3:
                print(f"TPC ADCs dataframe assembled {len(df_tpc)}x{len(df_tpc.columns)}")
        # if 'tp' in unpacked_tr:
        #     if verbosity >= 3:
        #         print("Assembling TPs")
        #     df_tp = pd.concat(unpacked_tr['tp'].values())
        #     df_tps.append(df_tp.sort_values(by=['time_start', 'channel']))
        #     if verbosity >= 3:
        #         print(f"TPs dataframe assembled: {len(df_tp)}")
    # Must not ignore index here - index stores the timestamp
    df_tpcs = pd.concat(df_tpcs)
    df_tpcs = df_tpcs.reindex(sorted(df_tpcs.columns), axis=1)
    if verbosity >= 1:
        print(f"Full TPC ADCs dataframe assembled {len(df_tpcs)}x{len(df_tpcs.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load ADCs: {time.time() - t_start:.2f}s")
    return df_tpcs


def read_roi_from_raw(
        file_path,
        cluster_row,
        buffer_channel=0, buffer_time=0,
        detector="VDColdbox",
        verbosity=0):
    """
    Get the raw data for the region containing the cluster referenced
    by in the input `cluster_row`.

    Extra area can be added to the region using the `buffer_channel`
    and `buffer_time` parameters.

    TODO buffer_channel could be smarter. Options should be:
     - Readout all channels
     - Readout all of the plane the cluster is on
     - Not go beyond the current plane unless specified
    This will likely need a new `multi_plane` parameter to be
    introduced.

    Parameters
    ----------
    file_path : str
        Location of the file contain the raw records to be read.
    cluster_row : pd.Series
        The row from a clusters dataframe (generated by
        `cluster_finder.create_basic_dataframe()`) as created by
        `df.loc[0]`. This should contain the range in channel and time,
        as well as a list of record ids which contain the raw
        records.
    buffer_channel : int, optional
        Number of extra channels to add either side of the region
        containing the cluster.
    buffer_time : int, optional
        Amount of extra time in ns to include either side of the region
        containing the cluster.
    
    Returns
    -------
    pd.DataFrame
        Dataframe with channels as columns and timestamps as the index
        indicating the recorded ADC of the channel at that time.
    """
    channel_min = cluster_row["channel_min"]
    channel_max = cluster_row["channel_max"]
    time_min = cluster_row["time_min"]
    time_max = cluster_row["time_max"]

    # TODO determine column name
    if "record_ids" not in cluster_row.keys():
        records_list = None
    else:
        records_list = cluster_row["record_ids"]

    return read_roi_from_raw_manual(
        file_path,
        channel_min, channel_max,
        time_min, time_max,
        records_list=records_list,
        buffer_channel=buffer_channel,
        buffer_time=buffer_time,
        detector=detector,
        verbosity=verbosity)


def read_roi_from_raw_manual(
        file_path,
        channel_min, channel_max,
        time_min, time_max,
        records_list=None,
        buffer_channel=0, buffer_time=0,
        detector="VDColdbox",
        verbosity=0):
    """TODO, as above, but manual max/min time/channel inputs."""
    # <0 alright since we simply use >= operator
    channel_min -= buffer_channel
    channel_max += buffer_channel
    time_min -= buffer_time
    time_max += buffer_time

    wethf_up = rdu.WIBEthFragmentPandasUnpacker(detector)
    up = rdu.UnpackerService()
    up.add_unpacker('bde_eth', wethf_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    if records_list is None:
        records_list = [ i for i,_ in rdf.get_all_record_ids()]

    df_tpcs = []
    if verbosity >= 4:
        t_start = time.time()
    for tr in records_list:
        if verbosity >= 2:
            print(f"--- Reading Record {tr} ---")
        
        # # This filtered unpacker ensures we only load frag data if it's
        # #   in the desired time window
        # unpacked_tr = up.unpack_filtered(rdf, tr, time_min, time_max)
        unpacked_tr = up.unpack(rdf, tr)

        if 'bde_eth' in unpacked_tr:
            dfs_bde = {k:v for k,v in unpacked_tr['bde_eth'].items() if not v is None}
            if verbosity >= 3:
                print(f"Assembling WIBEth Frames {len(dfs_bde)}")

            idx = pd.Index([], dtype='uint64')
            for df in dfs_bde.values():
                idx = idx.union(df.index)
            # Filters only events within desired time range
            time_filter = channel_filter = np.logical_and(
                idx >= time_min, idx <= time_max)
            if not np.any(time_filter):
                if verbosity >= 3:
                    print("No data within desired range found, skipping...")
                continue
            df_tpc = pd.DataFrame(index=idx[time_filter], dtype='uint16')

            for df in dfs_bde.values():
                # Filters data only in desired channel range
                # NOT SMART! Can only filter concurrence
                channel_filter = np.logical_and(
                    df.columns >= channel_min, df.columns <= channel_max)
                if np.any(channel_filter):
                    df_tpc = df_tpc.join(df.loc[time_filter, channel_filter])
            
            df_tpcs.append(df_tpc.reindex(sorted(df_tpc.columns), axis=1))

            if verbosity >= 3:
                print(f"TPC ADCs dataframe assembled {len(df_tpc)}x{len(df_tpc.columns)}")
        # if 'tp' in unpacked_tr:
        #     if verbosity >= 3:
        #         print("Assembling TPs")
        #     df_tp = pd.concat(unpacked_tr['tp'].values())
        #     df_tps.append(df_tp.sort_values(by=['time_start', 'channel']))
        #     if verbosity >= 3:
        #         print(f"TPs dataframe assembled: {len(df_tp)}")
    # Must not ignore index here - index stores the timestamp
    df_tpcs = pd.concat(df_tpcs)
    df_tpcs = df_tpcs.reindex(sorted(df_tpcs.columns), axis=1)
    if verbosity >= 1:
        print(f"Full TPC ADCs dataframe assembled {len(df_tpcs)}x{len(df_tpcs.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load ADCs: {time.time() - t_start:.2f}s")
    return df_tpcs