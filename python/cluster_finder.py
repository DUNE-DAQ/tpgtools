import trgtools
import numpy as np
import awkward as ak
from sklearn.cluster import dbscan
import pandas as pd
from scipy.spatial import ConvexHull

# UTILS, the following 4 functions ought to be moved to their own module...

def make_ak_slicer(array):
    """
    Given an input awkward array of integers, return a function which
    acts on a numpy array and returns an awkward array. The returned
    awkward array will be the same shape as the input awkward array to
    this function, where the values are replaced with the value at the
    corresponding index of the numpy array acted upon.

    This acts analogously to slicing a numpy aray using an array of
    integers.

    Example
    -------
    ```
    >>> ak_indicies = ak.Array([[2], [], [3,0]])
    >>> np_array = np.array([0, 10, 20, 30, 40])
    >>> slicer = make_ak_slicer(ak_indicies)
    >>> slicer(np_array)
    <Array [[20], [], [30, 0]] type='3 * var * int64'>

    Parameters
    ----------
    array : ak.Array
        Array of integer indicies for slicing.
    
    Returns
    -------
    function
        Function which takes a numpy input and returns an awkward array
        based on the slicing array.
    """
    n_points = ak.num(array, axis=1)
    flat_inds = ak.flatten(array)
    def slicer(arr):
        return ak.unflatten(arr[flat_inds], n_points)
    return slicer

def slice_to_ak(np_array, ak_indicies):
    """
    Given an input awkward array of integers, return an awkward array
    with the same shape as `ak_indicies`, where the values are
    replaced with th evalue at the value at the corresponding index of
    `np_array`.

    This acts analogously to slicing a numpy aray using an array of
    integers.

    Example
    -------
    ```
    >>> ak_indicies = ak.Array([[2], [], [3,0]])
    >>> np_array = np.array([0, 10, 20, 30, 40])
    >>> slice_to_ak(np_array, ak_indicies)
    <Array [[20], [], [30, 0]] type='3 * var * int64'>

    Parameters
    ----------
    np_array : np.ndarray
        Array to be sliced.
    ak_indicies : ak.Array
        Awkward array of integer indicies for slicing.
    
    Returns
    -------
    ak.Array
        Sliced view of `np_array`.
    """
    n_points = ak.num(ak_indicies, axis=1)
    flat_inds = ak.flatten(ak_indicies)
    return ak.unflatten(np_array[flat_inds], n_points)

def add_ak_event_offset(vals, return_n_bits=False):
    """
    Generates event uniqueness by bitshifting integer values to leave
    room for, and set the first n bits to reference the event number.

    If insufficient bits are present in the data type, the code will
    attempt to increase the size of the integer to hold the new bits.

    WARNING: does not necessaryily preserve datatype

    Parameters
    ----------
    vals : ak.Array
        Awkward array of some integer values.
    
    Returns
    -------
    ak.Array
        Awkward array of same shape of vals, where the index of the
        outer layer is refernce as the lowest n bits of the new values.
    """
    n_events = ak.num(vals, axis=0)
    n_extra_bits = len(bin(n_events)) - 2
    n_mother_bits = len(bin(ak.max(vals))) - 2
    offset_vals = vals
    if (n_extra_bits + n_mother_bits) > 31:
        offset_vals = ak.values_astype(offset_vals, np.int64)
    if (n_extra_bits + n_mother_bits) > 63:
        offset_vals = ak.values_astype(offset_vals, np.int128)
    if (n_extra_bits + n_mother_bits) > 127:
        raise ValueError("Overflow - numbers too large")
    new_array = (offset_vals << n_extra_bits) + np.arange(n_events)
    if return_n_bits:
        return new_array, n_extra_bits
    return new_array

def ak_unique_along_final_axis(array):
    """Returns unique values along final axis of awkward array."""
    offset_arr, n_bits = add_ak_event_offset(array, return_n_bits=True)
    unique_vals = np.unique(ak.ravel(offset_arr))
    cluster_nums = unique_vals % (2**n_bits)
    sorting = np.argsort(cluster_nums)
    # timeit implies this is ~3x faster than using
    #   np.unique(cluster_nums, return_counts=True)
    run_lengths = ak.run_lengths(cluster_nums[sorting])
    # Recover the initial values, sorted by event they appear in
    values = unique_vals[sorting] >> n_bits
    return ak.unflatten(values, run_lengths)


# Actual functions

def get_positions_array(
        tp_data,
        channel_factor=1.,
        time_factor=100.,
        time="time_peak"):
    """Return an array of channel/time coordinates with dtype np.uint"""
    return np.array([
            tp_data["channel"] / channel_factor,
            tp_data[time] / time_factor],
        dtype=np.uint).T

def db_cluster_tps(tp_data, epsilon=20, min_hits=7, channel_factor=1., time_factor=100.):
    """Create labels from DB clustering (sklearn)"""
    c_t_positions = get_positions_array(tp_data,
                                        channel_factor=channel_factor,
                                        time_factor=time_factor)
    _, all_labels = dbscan(c_t_positions, eps=epsilon, min_samples=min_hits)
    return all_labels

def create_clusters_array(cluster_labels):
    """
    Create an awkward array containg the indicies of the clusters
    indexed by `cluster_labels`. Axis 0 indexes which cluster, running
    from [0, max(cluster_labels)]. Axis 1 is a ragged axis containing
    the indicies of each TP in the cluster as integers.

    Parameters
    ----------
    cluster_labels : np.ndarray
        An array with the length of `tp_data` where each index is an
        integer label representing which cluster a TP belongs to (such
        that the number of clusters = 1 + max(cluster_labels)). -1 is
        used to index TP which do not belong in a cluster.
    
    Returns
    -------
    ak.Array
        Awkward array containing indicies for each cluster.
    """
    index, positions = np.unique(cluster_labels, return_inverse=True)
    
    # Unclustered vals are given label -1
    index = index[1:]
    small_mask = positions != 0 # 0 in position corresponds to -1 in index
    small_pos = positions[small_mask]
    small_tps = np.arange(positions.size)[small_mask]

    # Perform one sort as the complex step
    # This method is relatively slow for small counts,
    # but has ~ n log n scaling of quicksort for large n
    sort_args = np.argsort(small_pos)
    ak_size = ak.run_lengths(small_pos[sort_args])
    return ak.unflatten(small_tps[sort_args], ak_size)

# def decorate_ak_comprehension(func):
#     def new_func(tps, indicies):
#         pass
#     return new_func            

def create_basic_dataframe(tp_data, cluster_labels, extra_columns=[]):
    """
    Create a pandas dataframe containing the indicies of tp withing
    each cluster defined by `cluster_labels`.

    In addition, the dataframe contains the following properties:
    - `n_hits`: number of hits in the cluster.
    - `channel_min`: Lowest channel recorded in the cluster.
    - `channel_max`: Highest channel recorded in the cluster.
    - `time_max`: Minimum time tick recorded in the cluster. Calculated
        from `time_start`.
    - `time_max`: Maximum time tick recorded in the cluster. Calculated
        from `time_start` + `time_over_threshold`.

    Additional columns can be added using the `extra_columns`. This
    shold be a dictionary containing functions which act upon an array
    which is the TP data of all points in the cluster.

    NOTE on channels - currently this isn't smart. If a cluster manages
    to go over multiple APAs, potentially there will be a load of extra
    channels inbetween. Would need a smarter way of referncing the
    channels.

    Parameters
    ----------
    tp_data : np.ndarray or pd.Dataframe
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`, or generated by justintime)
    cluster_labels : np.ndarray
        An array with the length of `tp_data` where each index is an
        integer label representing which cluster a TP belongs to (such
        that the number of clusters = 1 + max(cluster_labels)). -1 is
        used to index TP which do not belong in a cluster.
    extra_columns : list [ function ], optional
        List containing as series of functions which add additional
        cluster properties to the dataframe columns.
        
        Fucntions should take an array of TPs, and awkward array of
        cluser indicies, and return a dictionary of type
        {str : array-like} containing the column name, and values for
        each cluster.
        
        A function which looks at one dimensions arrays of TPs can be
        converted to this format by decorating with the
        `decorate_ak_comprehension` function in this module.
    
    Returns
    -------
    pd.DataFrame
        Dataframe containing cluster indicies and properties.
    """
    clusters = create_clusters_array(cluster_labels)
    
    n_hits = ak.num(clusters, axis=1)

    data_indicies = ak.flatten(clusters)

    tp_peaks = ak.unflatten(tp_data["time_peak"][data_indicies], n_hits)
    tp_start_times = ak.unflatten(tp_data["time_start"][data_indicies], n_hits)
    tp_end_times = tp_start_times + ak.unflatten(tp_data["time_over_threshold"][data_indicies], n_hits)
    tp_channels = ak.unflatten(tp_data["channel"][data_indicies], n_hits)

    cols_dict = {
        "tps": clusters.to_list(),
        "n_hits": n_hits,
        "channel_min": ak.min(tp_channels, axis=1),
        "channel_max": ak.max(tp_channels, axis=1),
        "time_min": ak.min(tp_start_times, axis=1),
        "time_max": ak.max(tp_end_times, axis=1)}

    try:
        record_ids = ak.unflatten(tp_data["record_id"][data_indicies], n_hits)
        unique_records = ak_unique_along_final_axis(record_ids)
        cols_dict.update({"record_ids": unique_records.to_list()})
    except (KeyError, ValueError) as e:
        pass

    for func in extra_columns:
        cols_dict.update(func(tp_data, clusters))

    return pd.DataFrame(cols_dict)


def axes_finder(tps, indicies, channel_scale=1., time_scale=1.):
    """
    From a set of TPs and awkward array of cluster indicies, create a
    dictionary containing the major and minor axis sizes of all
    clusters based on the channel numebr and peak times.
    
    Parameters
    ----------
    tp_data : np.ndarray or pd.Dataframe
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`, or generated by justintime)
    indicies : ak.Array
        Array of indicies which form clusters.
    channel_scale : float, optional
        Scaling parameter applied to the channel, allows for moving the
        axis sizes to real space.
    time_scale : float, optional
        Scaling parameter applied to the times, allows for moving the
        axis sizes to real space.

    Returns
    -------
    dict {"major_axis": np.array, "minor_axis": np.array}
        Dictionary containg the major and minor axis values for all
        clusters.
    """
    slicer = make_ak_slicer(indicies)
    c = slicer(tps["channel"])
    t = slicer(tps["time_peak"])
    c = c - c[:, 0] - ak.mean(c - c[:, 0], axis=1)
    t = t - t[:, 0] - ak.mean(t - t[:, 0], axis=1)
    
    minor, major = fast_axis_approx(c, t)
    return {"minor_axis": minor, "major_axis": major}

def fast_axis_approx(x, y):
    """Get a fast approximation of the major/minor axis dimensions of a set of points"""
    # Given a convex hull, major and minor axes are the eignvalues
    # of the covariance of the hull points:
    # https://math.stackexchange.com/questions/207685/how-to-find-the-minimal-axis-parallel-ellipse-enclosing-a-set-of-points
    # This is an approximation by taking all points, strictly it ought to be just the outer points
    # See calc_axes_accurate for a precise version

    # This is the calculation of the covariance matrix
    xx = ak.var(x, axis=-1)
    yy = ak.var(y, axis=-1)
    xy = ak.covar(x, y, axis=-1)
    
    # And calculate eigenvalues using eigenvalue equation
    # [[xx-e, xy], [xy, yy-e]] => e^2 - (xx+yy)e + xx * yy - xy^2
    mid = (xx + yy)
    offset = np.sqrt(mid**2 - 4*(xx * yy - xy **2))
    minor = np.sqrt((mid - offset)/2)
    major = np.sqrt((mid + offset)/2)
    return minor, major

def axes_finder_accurate(tps, indicies, channel_scale=1., time_scale=1.):
    """
    From a set of TPs and awkward array of cluster indicies, create a
    dictionary containing the major and minor axis sizes of all
    clusters based on the channel numebr and peak times.
    
    Parameters
    ----------
    tp_data : np.ndarray or pd.Dataframe
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`, or generated by justintime)
    indicies : ak.Array
        Array of indicies which form clusters.
    channel_scale : float, optional
        Scaling parameter applied to the channel, allows for moving the
        axis sizes to real space.
    time_scale : float, optional
        Scaling parameter applied to the times, allows for moving the
        axis sizes to real space.

    Returns
    -------
    dict {"major_axis": np.array, "minor_axis": np.array}
        Dictionary containg the major and minor axis values for all
        clusters.
    """
    slicer = make_ak_slicer(indicies)
    c = ak.unflatten(slicer(tps["channel"]), 1, axis=1)
    t = ak.unflatten(slicer(tps["time_peak"]), 1, axis=1)
    c = c - c[:, 0] - ak.mean(c - c[:, 0], axis=1)
    t = t - t[:, 0] - ak.mean(t - t[:, 0], axis=1)
    coords = ak.concatenate((c*channel_scale, t*time_scale), axis=2)

    axes = ak.from_iter(map(calc_axes_accurate, coords))
    return {"minor_axis": axes[:,0], "major_axis": axes[:,1]}

def calc_axes_accurate(points):
    """Gets the major/minor axis dimensions of a set of points"""
    # TODO improve this!
    # The QJ option "joggles" the points by some ammount (up to 1??) to
    # ensure the convex hull can be formed, but this probably adds some
    # extra error - try to find a more exact (more efficient too?) way.
    hull = points[ConvexHull(points, qhull_options='QJ').vertices]
    # Given a convex hull, major and minor axes are the eignvalues
    # of the covariance of the hull points:
    # https://math.stackexchange.com/questions/207685/how-to-find-the-minimal-axis-parallel-ellipse-enclosing-a-set-of-points
    covar = np.cov(hull, rowvar=False, bias=True)
    eigvals = np.linalg.eigvals(covar)
    return np.sort(np.sqrt(eigvals))

def get_ind_matches(coll_ranges, starts, ends, ind_list, window):
    """ N.B requires TPs entirely within window"""
    match_ind_arr = np.arange(coll_ranges.shape[0])
    curr_group = 0
    groups = [np.array([])]*coll_ranges.shape[0]
    unmatched = []
    curr_inds = []
    init_i = 0
    for i, (start, end) in enumerate(zip(starts, ends)):
        if end > (coll_ranges[curr_group, 1] + window):
            # Useful printouts for debugging, uncomment if wanted.
            #print(f"Next group: {end} > {(coll_ranges[curr_group, 1] + window)}")
            #print(f"Added group: ind_list[{init_i}:{i}]")
            #print(ind_list[init_i:i])
            groups[curr_group] = ind_list[init_i:i]
            # Could do something mildly fancy here to avoid checking groups < curr_group
            curr_group = match_ind_arr[(coll_ranges[:, 1] + window) >= end][0]
            init_i = i
        if start >= (coll_ranges[curr_group, 0] - window):
            #print(f"Added to current group: {end} > {(coll_ranges[curr_group, 1] + window)} and {start} >= {(coll_ranges[curr_group, 0] - window)}")
            continue
        # Probably could remove this too with appropriate ordering of if statements,
        # or it may just require a "last tp not in range" flag, which I can't be bothered to add...
        #print(f"Unmatched: {start} >= {(coll_ranges[curr_group, 0] - window)}")
        unmatched += [ind_list[i]]
        init_i = i+1
    return ak.Array(groups), np.array(unmatched)

def match_all_tps_in_windows(
        coll_tps, u_tps, v_tps,
        coll_gap, u_window, v_window,
        u_offset=375, v_offset=188,
        return_ranges=False):
    """
    Match together TPs on all planes based on closeness in time.
    
    A group of collection TPs is constructed if they are separated by
    at most `coll_gap` ticks. Then, the same window is applied to the
    TPs on the induction planes, offset by the expected number of ticks
    for charge to drift to the collection plane. This window can be
    slightly widened to account for uncertainty in drift time using the
    `u_window` and `v_window` arguments. If an induction TP is
    contained completely withing this window, it is added into the
    group.

    BEWARE: it does not consider any spatial coordinates. Matching is
    done purely on time. This is appropriate for small clusters of TPs
    when the supplied tps exist only on a few neighbouring channels. To
    try to match a cluster, the `coll_tps` would need to contain only
    the TPs contained in that cluster as determined by some other
    algorithmn.

    Parameters
    ----------
    coll_tps : pd.DataFrame
        DataFrame containing TPs from the collection plane which should
        be matched to induction planes. Collection TPs are grouped
        together if the diffence between the last seen tps <=
        `coll_gap`, regardless of wire. Ensure this does not contain
        TPs appear at the same time, but with large channel
        separations.
    u_tps : pd.DataFrame
        DataFrame of TPs generated from the U plane (plane 0).
    v_tps : pd.DataFrame
        DataFrame of TPs generated from the V plane (plane 1).
    coll_gap : int
        Number of ticks (16ns) that must have passed since the previous
        collection TP to create a new group. For saftey, this should be
        at least `2 * max(u_window, v_window)`
    u_window : int
        Number of extra ticks (16ns) leeway added to the U plane window
        when deciding if a TP should belong to the group. 1 tick should
        be sufficient in most cases. If there are signifcant difference
        in grouping with different values, it may indicate the
        `u_offset` is incorrect.
    u_window : int
        Number of extra ticks (16ns) leeway added to the V plane window
        when deciding if a TP should belong to the group. 1 tick should
        be sufficient in most cases. If there are signifcant difference
        in grouping with different values, it may indicate the
        `v_offset` is incorrect.
    u_offset : int, optional
        Number of ticks (16ns) expected for charge to drift from the
        U-plane to the collection plane. PD-II horizontal drift, this
        is expected to be 375 ticks. Default is 375.
    u_offset : int, optional
        Number of ticks (16ns) expected for charge to drift from the
        V-plane to the collection plane. PD-II horizontal drift, this
        is expected to be 188 ticks. Default is 188.
    return_ranges : bool, optional
        If true, return the time ranges of the collection planes TPs
        used in the groupings as an `(n, 2)`-shape array, for `n`
        groups found. This could be manually calculated per group as
        `[min(tps.time_start),
         max(tps.time_start + tps.time_over_threshold)]`. Dfeault is
        False.
    
    Returns
    -------
    `(coll_groups, u_matched, v_matched, u_unmatched, v_unmatched)` or
    `coll_groups, u_matched, v_matched, u_unmatched, v_unmatched,
     matches)` if `return_ranges` is True.
    
    coll_groups : ak.Array
        Awkward array of integers which pick out the TPs in `coll_tps`
        that form each group. Each entry in the 0 axis corresponds to
        one TP group.
    u_matched : ak.Array
        Awkward array of integers which pick out the TPs in `u_tps`
        that form each group. Each entry in the 0 axis corresponds to
        one TP group, matching with the same outer index of
        `coll_groups`.
    v_matched : ak.Array
        Awkward array of integers which pick out the TPs in `v_tps`
        that form each group. Each entry in the 0 axis corresponds to
        one TP group, matching with the same outer index of
        `coll_groups`.
    u_unmatched : np.ndarray
        Numnpy array of integers which pick out the TPs in `u_tps`
        which are not matched to any group.
    v_unmatched : np.ndarray
        Numnpy array of integers which pick out the TPs in `v_tps`
        which are not matched to any group.
    matches : np.ndarray, optional
        Only returned if `return_ranges` is True. Shape is `(n, 2)`,
        for `n` groups found. The first column of axis 1 is the
        minimum time start of the group, the second column of axis 1 is
        the maximum time end of the group. Axis 0 matches to axis 0 of
        the `coll_groups` array.
    """
    used_collections = np.full(len(coll_tps), False)
    coll_starts = coll_tps["time_start"]
    u_starts = u_tps["time_start"]
    v_starts = v_tps["time_start"]
    min_ts = min(np.min(coll_tps["time_start"]), np.min(u_tps["time_start"]), np.min(v_tps["time_start"]))
    def ticks(times):
        return (times.to_numpy() - min_ts)/16
    coll_sort = np.argsort(coll_starts)
    coll_starts = ticks(coll_starts)[coll_sort]
    coll_ends = ticks(coll_tps["time_start"] + coll_tps["time_over_threshold"])
    coll_list = np.arange(coll_starts.size)[coll_sort]
    u_sort = np.argsort(u_starts)
    u_starts = ticks(u_starts)[u_sort] - u_offset
    u_ends = ticks(u_tps["time_start"] + u_tps["time_over_threshold"])[u_sort] - u_offset
    u_list = np.arange(u_starts.size)[u_sort]
    v_sort = np.argsort(v_starts)
    v_starts = ticks(v_starts)[v_sort] - v_offset
    v_ends = ticks(v_tps["time_start"] + v_tps["time_over_threshold"])[v_sort] - v_offset
    v_list = np.arange(v_starts.size)[v_sort]

    matches = []
    coll_groups = []
    init_i = 0
    curr_start = coll_starts[0]
    last_end = coll_starts[0]
    for i, (start, end) in enumerate(zip(coll_starts, coll_ends)):
        if start - last_end > coll_gap:
            matches += [np.array([curr_start, end])]
            coll_groups += [coll_list[init_i:i]]
            curr_start = start
            init_i = i
        last_end = end
    matches += [np.array([curr_start, coll_ends[-1]])]
    matches = np.array(matches)
    coll_groups += [coll_list[init_i:i]]
    coll_groups = ak.Array(coll_groups)
    
    u_matched, u_unmatched = get_ind_matches(matches, u_starts, u_ends, u_list, u_window)
    v_matched, v_unmatched = get_ind_matches(matches, v_starts, v_ends, v_list, v_window)

    if return_ranges:
        return coll_groups, u_matched, v_matched, u_unmatched, v_unmatched, matches
    return coll_groups, u_matched, v_matched, u_unmatched, v_unmatched
