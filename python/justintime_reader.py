import numpy as np
import pandas as pd
import hdf5libs
import justintime.utils.rawdataunpacker as rdu
import tpgsandbox.emulation.algos as tpalgos
from scipy.spatial import ConvexHull
import time
import warnings


def _unpacking_core(
        file, unpacker,
        frags_list=None, init_frag=0, n_frags=1,
        channel_bounds=None, time_bounds=None,
        verbosity=2):
    """
    Create a pandas dataframe contain raw data.
    NOTE: only one unpacker allowed in the unpacker ("tp" or
    "wib_eth").

    Parameters
    ----------
    file : HDF5RawDataFile
        Path to hdf5 file to read.
    unpacker : UnpackerService
        Unpacker service with unpackers already added.
    frags_list : list [int], optional
        If specified, read only the frags in the frags_list, else read
        all frags. Default is None
    init_frag : int, optional
        Index of the first frag to start reading from. Only used if
        frags_list is not specified. Default is 0.
    n_frags : int, optional
        Number of frags to read. If -1, read out all frags after
        init_frag. Only used if frags_list is not specified. Default is
        1.
    channel_bounds : tuple[int, int], optional
        (minimum, maximum) channels to be read out, inclusive of
        minimum and maximum. If None, read all channels. Default is
        None.
    time_bounds : tuple[int, int], optional
        (minimum, maximum) time ticks to be read out, inclusive of
        minimum and maximum. If None, read all ticks. Default is None.
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 1.
    
    Returns
    -------
    pd.Dataframe
        DataFrame containing the unpacked data from the file.
    """
    record_type = file.is_trigger_record_type()

    if frags_list is None:
        if record_type:
            frags_list = [ i for i,_ in file.get_all_record_ids()]
        else:
            frags_list = [ i for i,_ in file.get_all_timeslice_ids()]
        # process only the first TR
        frag_end = init_frag + n_frags if n_frags != -1 else -1
        frags_list = frags_list[init_frag:frag_end]

    type_str = "Record" if record_type else "TimeSlice"
    if "tp" in unpacker.fragment_unpackers.keys():
        return _tps_unpacker(
            file, unpacker,
            frags_list, type_str,
            verbosity=verbosity)
    elif 'wib_eth' in unpacked_trig:
        channel_fn = _make_bounds_func(channel_bounds)
        time_fn = _make_bounds_func(time_bounds)
        return _raw_unpacker(
            file, unpacker, frags_list, type_str,
            channel_filter_fn=channel_fn, time_filter_fn=time_fn,
            verbosity=verbosity)
    else:
        raise ValueError(
            "Cannot find either tp or wib_eth data in the unpacker.")

def _make_bounds_func(bounds):
    """Create a function which selects data within the bounds."""
    if bounds is None:
        return bounds
    def func(data):
        return np.logical_and(
            data >= bounds[0], data <= bounds[1])

def _tps_unpacker(file, unpacker, frags_list, type_str, verbosity=2):
    """Core unpacker for TP type data"""
    df_data = []
    if verbosity >= 4:
        t_start = time.time()
    for trig in frags_list:
        if verbosity >= 2:
            print(f"--- Reading {type_str} {trig} ---")

        unpacked_trig = unpacker.unpack(file, trig)
        if verbosity >= 3:
            print("Assembling TPs")
        df_tp = pd.concat(unpacked_trig['tp'].values())
        df_data.append(df_tp.sort_values(by=['time_start', 'channel']))
        if verbosity >= 3:
            print(f"TPs dataframe assembled: {len(df_tp)}")
    
    df_data = pd.concat(df_data, ignore_index=True)
    df_data = df_data.reindex(sorted(df_data.columns), axis=1)
    if verbosity >= 1:
        print(f"Full TPs dataframe assembled {len(df_data)}x{len(df_data.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load TPs: {time.time() - t_start:.2f}s")
    return df_data

def _raw_unpacker(
        file, unpacker, frags_list, type_str,
        channel_filter_fn=None, time_filter_fn=None,
        verbosity=2):
    """
    Core unpacker for raw data.

    Includes time and channel filters, which are passed as None for no
    filters, or a tuple (minimum, maximum) to select only time
    ticks/channels in range [minimum, maximum].
    """
    df_data = []
    if verbosity >= 4:
        t_start = time.time()
    for trig in frags_list:
        dfs_bde = {k:v for k,v in unpacked_trig['wib_eth'].items() if not v is None}
        if verbosity >= 3:
            print(f"Assembling WIBEth Frames {len(dfs_bde)}")

        idx = pd.Index([], dtype='uint64')
        for df in dfs_bde.values():
            idx = idx.union(df.index)
        # Filters only events within desired time range
        if time_filter_fn is not None:
            time_filter = time_filter_fn(idx)
            if not np.any(time_filter):
                if verbosity >= 3:
                    print("No data within desired range found, skipping...")
                continue
        else:
            time_filter = slice(None)

        df_adc = pd.DataFrame(index=idx, dtype='uint16')

        for df in dfs_bde.values():
            if channel_filter_fn is not None:
                channel_filter = channel_filter_fn(df.columns)
                if np.any(channel_filter):
                    df = df.loc[time_filter, channel_filter]
            else:
                df = df.loc[time_filter]
            df_adc = df_adc.join(df)

        df_data.append(df_adc.reindex(sorted(df_adc.columns), axis=1))

        if verbosity >= 3:
            print(f"TPC ADCs dataframe assembled {len(df_adc)}x{len(df_adc.columns)}")
    df_data = pd.concat(df_data)
    df_data = df_data.reindex(sorted(df_data.columns), axis=1)
    if verbosity >= 1:
        print(f"Full ADCs dataframe assembled {len(df_data)}x{len(df_data.columns)}")
    if verbosity >= 4:
        print(f"Time taken to load ADCs: {time.time() - t_start:.2f}s")
    return df_data

def read_hdf5_tps_jit(file_path, init_slice=0, n_slices=1, detector="PD2HD", verbosity=2):
    """
    Create a pandas dataframe of TPs using justintime.

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
        by the Channel Map. Default is "PD2HD".
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

    return _unpacking_core(rdf, up,
                           init_frag=init_slice, n_frags=n_slices,
                           verbosity=verbosity)

def read_hdf5_raw_jit(
        file_path,
        records_list=None, init_record=0, n_records=1,
        detector="PD2HD",
        verbosity=1):
    """
    Create a pandas dataframe contain raw data.

    Parameters
    ----------
    file_path : str
        Path to hdf5 file to read.
    records_list : list[int], optional
        List of record IDs to be read out. Default is None.
    init_record : int, optional
        Index of the first Record to start reading from. Default is 0.
    n_records : int, optional
        Number of Recordss to read. If -1, read out all Records after
        init_record. Default is 1.
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "PD2HD".
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
    up.add_unpacker('wib_eth', wethf_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    return _unpacking_core(rdf, up,
                           frags_list=records_list,
                           init_frag=init_record, n_frags=n_records,
                           verbosity=verbosity)

def read_roi_from_raw(
        file_path,
        cluster_row,
        buffer_channel=0, buffer_time=0,
        detector="PD2HD",
        verbosity=0):
    """
    Get the raw data for the region containing the cluster referenced
    by in the input `cluster_row`.

    Extra area can be added to the region using the `buffer_channel`
    and `buffer_time` parameters.

    NOTE buffer_channel could be smarter. Options should be:
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
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "PD2HD".
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 0.
    
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
        (channel_min, channel_max),
        (time_min, time_max),
        records_list=records_list,
        buffer_channel=buffer_channel,
        buffer_time=buffer_time,
        detector=detector,
        verbosity=verbosity)


def read_roi_from_raw_manual(
        file_path,
        channel_bounds,
        time_bounds,
        records_list=None,
        buffer_channel=0, buffer_time=0,
        detector="VDColdbox",
        verbosity=0):
    """
    Get the raw data for the region containing the cluster referenced
    by in the input `cluster_row`.

    Extra area can be added to the region using the `buffer_channel`
    and `buffer_time` parameters.

    Parameters
    ----------
    file_path : str
        Location of the file contain the raw records to be read.
    records_list : list[int], optional
        List of record IDs to be read out. Default is None.
    channel_bounds : tuple[int, int], optional
        (minimum, maximum) channels to be read out, inclusive of
        minimum and maximum. If None, read all channels. Default is
        None.
    time_bounds : tuple[int, int], optional
        (minimum, maximum) time ticks to be read out, inclusive of
        minimum and maximum. If None, read all ticks. Default is None.
    buffer_channel : int, optional
        Number of extra channels to add either side of the region
        containing the cluster.
    buffer_time : int, optional
        Amount of extra time in ns to include either side of the region
        containing the cluster.
    detector : str, optional
        Name of the detector the hdf5 file corresponds to, referenced
        by the Channel Map. Default is "PD2HD".
    verbosity : int, optional
        Amount of text to print, 0 for no outputs. Default is 0.
    
    Returns
    -------
    pd.DataFrame
        Dataframe with channels as columns and timestamps as the index
        indicating the recorded ADC of the channel at that time.
    """
    # <0 alright since we simply use >= operator
    channel_min = channel_bounds[0] - buffer_channel
    channel_max = channel_bounds[1] + buffer_channel
    time_min = time_bounds[0] - buffer_time
    time_max = time_bounds[1] + buffer_time
    channel_bounds = (channel_min, channel_max)
    time_bounds = (time_min, time_max)

    wethf_up = rdu.WIBEthFragmentPandasUnpacker(detector)
    up = rdu.UnpackerService()
    up.add_unpacker('wib_eth', wethf_up)
    if verbosity >= 1:
        print(f"Opening {file_path}")
    rdf = hdf5libs.HDF5RawDataFile(file_path)

    return _unpacking_core(
        rdf, up,
        frags_list=records_list, n_frags=-1,
        channel_bounds=channel_bounds, time_bounds=time_bounds,
        verbosity=verbosity)

def get_frag_types_in_file(file_path):
    """Returns the set of types contained within the passed hdf5 file"""
    rdf = hdf5libs.HDF5RawDataFile(file_path)
    results = []
    for id in rdf.get_all_record_ids():
        for path in rdf.get_fragment_dataset_paths(id):
            results += [rdf.get_frag(path).get_fragment_type()]
    return set(results)

def emulate_tps_from_raw(
        raw_df,
        running_sum_ratio, threshold,
        ped_acc_limit=10,
        init_ped_range=100,
        init_ped_algo="mode", running_sum_algo="standard",
        input_stage=None, return_stage="tps",
        detector="PD2HD", offset_times=False):
    """
    Takes a DataFrame of raw ADC values and creates a DataFrame of
    emulated TPs.

    TPs are emulated via the following procedure:
    `"format"`:
        Format the DataFrame to use 16 bit integer values for ADCs, and
        64 bit integer values for timestamps. If `offset_times` is set
        as True, the lowest timestamp is set as 0.
    `"ped"`:
        The channels are pedestal subtracted (high pass filter). This
        is done using an accumulator which will update the pedestal
        when the positive or negative `ped_acc_limit` is reached,
        attemping to estimate the median.
        The initial expected median is estimated by a separate
        `init_ped_algo` alogorithm (`"mode"` or `"mean"`), based on the
        initial `init_ped_algo` time ticks of each channel.
    `"sum"`:
        A running sum is applied to the channels (low pass filter). The
        decay coefficient (ratio of previous value to keep) is
        configurable as `running_sum_ratio`. The running sum algorithm
        is configurable. Known sums are: `"standard"`.
    `"tps"`:
        TPs are generated by searching for ADC values above a
        threshold, configurable as `threshold`.
    
    Which stage the input data starts from, and the stage at which to
    output data is configurable via the `input_stage` and
    `return_stage` arguments, by passing the strings refering to each
    stage above (or None for raw input data without any processing).

    Parameters
    ----------
    raw_df : pd.DataFrame
        DataFrame of raw ADC values with columns as channel numbers and
        index as timestamps.
    running_sum_ratio : float
        Number between 0. and 1. indicating the fraction of the
        previous running sum value to be kept.
    threshold : int
        ADC count threshold which must be reached to generate a TP.
    ped_acc_limit : int, optional
        Number of ticks above or below the current median required to
        cause an update to the pedestal. Default is 10.
    init_ped_range : int, optional
        Number of ticks to use to estimate the initial pedestal value.
        Default is 100.
    init_ped_algo : str ["mode", "mean"], optional
        Method to generate initial pedestal estimate. Default is
        "mode".
    running_sum_algo : str ["standard"], optional
        Running sum method to use. Default is "standard".
    intput_stage : str [None, "format", "ped", "sum"], optional
        What stage the input DataFrame is at (i.e. what was the last
        step applied to it). If None, DataFrame is assumed to be raw
        values from an hdf5 file with no processing applied. Default is
        None.
    return_stage : str ["format", "ped", "sum", "tps"], optional
        What stage to return the DataFrame. Can be used to investigate
        the TP generation at iterim stages. Default is "tps".
    detector : str, optional
        Which detector the data is from (to generate channel map). The
        channel map should be discoverable as
        `detector + "ChannelMap"`. Default is "PD2HD".
    offset_times : bool, optional
        If True, offset the lowest timestamp of the data to 0, else
        keep true timestamps. Default is False.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing the data at the stage specified by
        `return_stage`. By default, this is a DataFrame containing TPs.
    """
    return_stages = ["format", "ped", "sum", "tps"]
    if return_stage not in return_stages:
        raise ValueError(
            f"return_stage: {return_stage}, must be one of {return_stages}")
    input_stages = [None] + input_stages[:-1]
    if input_stage not in return_stages:
        raise ValueError(
            f"input_stage: {input_stage}, must be one of {input_stages}")
    init_ped_algos = ["mode", "mean"]
    if init_ped_algo not in init_ped_algos:
        raise ValueError(
            f"init_ped_algo: {init_ped_algo}, "
            + f"must be one of {init_ped_algos}")
    running_sum_algos = ["standard"]
    if running_sum_algo not in running_sum_algos:
        raise ValueError(
            f"running_sum_algo: {running_sum_algo}, "
            + f"must be one of {running_sum_algos}")
    input_ind = input_stages.index(input_stage)
    return_ind = return_stages.index(return_stage)
    if input_ind >= return_ind:
        raise ValueError(f"input_stage {input_stage} must be earlier "
                         + f"than the return_stage {return_stage}")

    if input_stage <= 0:
        formatted_raw = raw_df.astype('int16')
        t0 = formatted_raw.index.min() if offset_times else 0
        formatted_raw.index = formatted_raw.index.astype('int64') - t0
        if return_stage == "format":
            return formatted_raw

    if input_stage <=1: 
        df_ped, df_ped_var = tpalgos.emulate_ped(
            formatted_raw, limit=ped_acc_limit, init_ped_range=init_ped_range)
        df_adc = formatted_raw-df_ped
        if return_stage == "ped":
            return df_adc

    if input_stage <= 2:
        match running_sum_algo:
            case "standard":
                df_rs_adc = tpalgos.emulate_running_sum(df_adc, running_sum_ratio)
        if return_stage == "sum":
            return df_rs_adc

    channel_map = detchannelmaps.make_map(f"{detector}ChannelMap")
    df_emu_tps = tpalgos.generate_tps(df_adc, threshold, channel_map)
    return df_emu_tps
