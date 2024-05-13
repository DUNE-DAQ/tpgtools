#!/usr/bin/env python
"""
For a given TPStream, plot the channel histogram and print
the channels that are above a given limit and Tuckey's fence.
"""


from trgtools import TPReader

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os


# As of 05-13-2024.
DETECTOR_CHANNELMAPS = {
            "HD": os.path.join(os.environ["DETCHANNELMAPS_SHARE"], "config/pd2hd/PD2HDChannelMap_v6.txt"),
            "VD": os.path.join(os.environ["DETCHANNELMAPS_SHARE"], "config/pd2vd/PD2VDBottomTPCChannelMap_v1.txt"),
            "VDCB": os.path.join(os.environ["DETCHANNELMAPS_SHARE"], "config/vdcoldbox/vdcbce_chanmap_v4.txt")
        }

# Colors.
BLUE = "#63ACBE"
RED  = "#EE442F"


def get_detector_element_channel_hist_map(map_name: str) -> dict[str, int]:
    """
    Get a map with the detector elements as the key and the channels (as histogram edges)
    as the value.

    Parameters:
        map_name (str): Channel map name to load from.

    Returns a dict with the detector element name and channel histogram edges.
    """
    # Only care about the offline channels (0) and detector element name (2) columns.
    df = pd.read_table(DETECTOR_CHANNELMAPS[map_name], sep="\s+", usecols=[0, 2], header=None)
    named_series = pd.Series(df[0].values, index=df[2])  # Many channels will belong to the same det_elem.
    unique_det_elems = named_series.keys().unique()      # So the det_elem will repeat many times.

    det_elem_edges = {}
    for det_elem in unique_det_elems:
        min_channel = named_series[det_elem].values.min()
        max_channel = named_series[det_elem].values.max()
        det_elem_edges[det_elem] = np.arange(min_channel-0.5, max_channel+1)  # Need padding for hist edges.

    return det_elem_edges


def read_requested_fragments(reader: TPReader, fragment_index: int, num_fragments: int) -> None:
    """
    Read the requested fragments and make changes if necessary.

    Parameters:
        reader          (TPReader): TPReader object to read fragments from.
        fragment_index  (int)     : First fragment path to read from.
        num_fragments   (int)     : Number of fragments to read.

    Returns nothing. Mutates :reader:.
    """
    len_paths = len(reader.get_fragment_paths())
    if fragment_index > len_paths:
        raise IndexError(f"There are only {len_paths} fragment paths. Choose a smaller index.")
    if fragment_index + num_fragments > len_paths:
        num_fragments = len_paths - fragment_index
        print( "WARNING: fragment_index + num_fragments > the number of fragments in the file.")
        print(f"         Using num_fragments = {num_fragments} instead.")

    reader.set_fragment_paths(reader.get_fragment_paths()[fragment_index:fragment_index+num_fragments])
    reader.read_all_fragments()
    return


def get_tukey_upper_fence(data: np.ndarray, k: float = 1.5) -> np.float64:
    """
    Calculate the upper Tukey's fence for the given k.

    Parameters:
        data (np.ndarray): Data sample to calculate the fence for.
        k    (float)     : Fence extension value. k = 1.5 gives Tukey's proposed 'outlier' limit.
    """
    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25
    return k * iqr + q75


def plot_channel_histogram(hist: np.ndarray,
                           bins: np.ndarray,
                           limit: int,
                           tukeys_fence: float,
                           run_id: int,
                           file_index: int,
                           det_elem: str) -> None:
    """
    Plots the very fine bin sized channel histogram outlier limit lines.

    Parameters:
        hist         (np.ndarray): Calculated histogram counts to plot.
        bins         (np.ndarray): Calculated histogram bin edges to plot.
        limit        (int)       : The user-given limit.
        tukeys_fence (float)     : The Tukey's upper fence limit.
        run_id       (int)       : HDF5 run ID.
        file_index   (int)       : HDF5 file index.
        det_elem     (str)       : String of the detector element this is plotting.

    Returns nothing. Saves a plot in the current directory. Overwrites plots
    generated from the same HDF5.
    """
    plt.figure(figsize=(6, 4), dpi=200)  # Needs to be quick, so hard setting PNGs.
    ax = plt.gca()

    # Plot in linear style.
    ax.stairs(hist, edges=bins, color=BLUE, alpha=0.6, label="Linear", fill=True)
    ax.set_yscale("linear")

    # Plot in log style.
    ax2 = ax.twinx()
    ax2.stairs(hist, edges=bins, color=RED, alpha=0.6, label="Log", fill=True)
    ax2.set_yscale("log")

    # Plot the limits. Choosing log to show this. Maybe should be configurable.
    ax2.axhline(limit, color=RED, label=f"User Limit (Log): {limit}")
    ax2.axhline(tukeys_fence, ls=(0, (5, 5)), color=RED, label=f"Tukey's Fence (Log): {tukeys_fence}")

    # Set axis and tick colors.
    ax.spines['left'].set_color(BLUE)
    ax.yaxis.label.set_color(BLUE)
    ax.tick_params('y', colors=BLUE)

    ax.spines['right'].set_color(RED)  # Actually belongs to ax and not ax2.
    ax2.yaxis.label.set_color(RED)
    ax2.tick_params('y', colors=RED)

    # Set the plot order.
    ax.set_zorder(2)
    ax.patch.set_visible(False)
    ax2.set_zorder(1)

    # Set the legend.
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles + handles2
    labels = labels + labels2
    plt.legend(handles=handles, labels=labels)

    # Title and axes labels.
    plt.title(f"Overactive TP Channels In {det_elem}")
    ax.set_xlabel("Channel Number")
    ax.set_ylabel("TriggerPrimitive Count (Linear)")
    ax2.set_ylabel("TriggerPrimitive Count (Log)")

    # Finally save and close.
    plt.tight_layout()
    plt.savefig(f"tpg_overactive_channels-{run_id}.{file_index:04}-{det_elem}.png")
    plt.close()
    return


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument("file", type=click.Path(exists=True, readable=True))
@click.option('-d', "--detector", type=click.Choice(["HD", "VD", "VDCB"], case_sensitive=False),
              help="Select the detector to process.", required=True, prompt=True)
@click.option('-i', "--fragment-index", type=click.INT, default=10,
              help="Fragment index to start processing on.")
@click.option('-k', "--outlier-multiplier", "tukey_k", type=click.FLOAT, default=1.5,
              help="Outlier multiplier for Tukey's fences.")
@click.option('-l', "--limit", type=click.INT, default=1000,
              help="User limit to mark overactive channels.")
@click.option('-n', "--num-fragments", type=click.INT, default=1,
              help="The number of fragments to process.")
@click.option('-v', "--verbose", count=True,
              help="Be more verbose about fragment reading.")
def main(file, detector, fragment_index, num_fragments, limit, tukey_k, verbose):
    tp_reader = TPReader(file, verbose)
    channel_hist_map = get_detector_element_channel_hist_map(detector)

    # Fragment processing logics.
    read_requested_fragments(tp_reader, fragment_index, num_fragments)

    # Process the channel histogram.
    channels = tp_reader.tp_data['channel']

    for det_elem, bins in channel_hist_map.items():
        print("Detector Element:", det_elem)
        hist, _ = np.histogram(channels, bins=bins)
        if np.all(hist == 0):  # No activity in this det_elem for the loaded fragments.
            print("No activity! Try a different fragment or investigate more thoroughly.")
            print("=" * 60)
            continue

        print(f"Channels Above Limit ({limit}):\n",
              np.array2string(np.where(hist > limit)[0], separator=',', threshold=np.inf))

        tukeys_fence = get_tukey_upper_fence(hist, tukey_k)
        print(f"Channels Above Tukey's Fence ({tukeys_fence}):\n",
              np.array2string(np.where(hist > tukeys_fence)[0], separator=',', threshold=np.inf))

        plot_channel_histogram(hist, bins, limit, tukeys_fence, tp_reader.run_id, tp_reader.file_index, det_elem)

        print("=" * 60)
    return


if __name__ == "__main__":
    main()
