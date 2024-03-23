import trgtools
import numpy as np
import awkward as ak
from sklearn.cluster import dbscan
import pandas as pd
from scipy.spatial import ConvexHull


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
    n_points = ak.num(ak_slices, axis=1)
    flat_inds = ak.flatten(ak_slices)
    return ak.unflatten(np_array[flat_inds], n_points)

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

def decorate_ak_comprehension(func):
    def new_func(tps, indicies):
        pass
    return new_func

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

    Parameters
    ----------
    tp_data : np.ndarray
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`)
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
    tp_channels = ak.unflatten(tp_data["time_peak"][data_indicies], n_hits)

    cols_dict = {
        "tps": clusters.to_list(),
        "n_hits": n_hits,
        "channel_min": ak.min(clusters, axis=1),
        "channel_max": ak.max(clusters, axis=1),
        "time_min": ak.min(tp_start_times, axis=1),
        "time_max": ak.max(tp_end_times, axis=1)}

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
    tp_data : np.ndarray
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`)
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
    tp_data : np.ndarray
        Structed array containing the TP data (i.e. as generated by
        `trgtools.TPReader`)
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