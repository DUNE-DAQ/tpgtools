# tpgtools 
Here is a short summary of the applications and scripts available in `tpgtools` 

## Emulator

`wibeth_tpg_algorithms_emulator` is an emulator application for TPG algorithms. It takes as input a Trigger Record file (`*.hdf5` file) and it will execute the selected TPG algorithm on the Trigger Record data. The application is single threaded, pinned to core 0. The core number is configurable.   

To use the tool use the following:
```sh
$ wibeth_tpg_algorithms_emulator --help 
TPG algorithms emulator using input from Trigger Record files
Usage: wibeth_tpg_algorithms_emulator [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--file-path-input TEXT   Path to the input file
  -a,--algorithm TEXT         TPG Algorithm (SimpleThreshold / AbsRS)
  -i,--implementation TEXT    TPG implementation (AVX / NAIVE). Default: AVX
  -m,--channel-map TEXT       Select a valid channel map: None, VDColdboxChannelMap, ProtoDUNESP1ChannelMap, PD2HDChannelMap, HDColdboxChannelMap, FiftyLChannelMap
  -n,--num-TR-to-read INT     Number of Trigger Records to read. Default: select all TRs.
  -t,--tpg-threshold INT      Value of the TPG threshold. Default value is 500.
  -c,--core INT               Set core number of the executing TPG thread. Default value is 0.
  --save-adc-data             Save ADC data (first frame only)
  --save-trigprim             Save trigger primitive data
  --parse_trigger_primitive   Parse Trigger Primitive records
  -s,--num-frames-to-save INT Set the number of frames of ADC data from a TR to save: -1 (all) or 1. Default: 1 (first frame only).
```

The command line option `save_adc_data` allows to save the raw ADC values in a txt file after the 14-bit to 16-bit expansion. The command line option `save_trigprim`  allows to save the in a file the Trigger Primitive object information in a txt file. 

Example of usage: 
```sh
$ wibeth_tpg_algorithms_emulator -f swtest_run000035_0000_dataflow0_datawriter_0_20231102T083908.hdf5  -a SimpleThreshold -m PD2HDChannelMap -t 500 --save-trigprim --parse_trigger_primitive
$ wibeth_tpg_algorithms_emulator -f swtest_run000035_0000_dataflow0_datawriter_0_20231102T083908.hdf5  -a AbsRS -m PD2HDChannelMap -t 500 --save-adc-data  -n 5 
```


## Utility tools and scripts

* `wibeth_tpg_workload_emulator` is a simple emulator for TPG algorithm in either a naive or in AVX implementation. The application allows to emulate the workload when running a TPG algorithm and therefore monitor performance metrics. It requires an input binary frame file (check assets-list for valid input files) and it will execute the desired TPG algorithm for a configurable duration (default value is 120 seconds). The application is single threaded, pinned to core 0 (configurable). Check the helper page for more details. 

Example of usage: `wibeth_tpg_workload_emulator -f wibeth_frame_file.bin` 

* `wibeth_tpg_validation` is a simple emulator for validating different TPG algorithms, either in naive or in AVX implementation. Check the helper page for more details. Usage: `wibeth_tpg_validation -f wibeth_frame_file.bin` 

* `wibeth_binary_frame_reader`: reads a WIBEth frame file (`.bin` file) and prints all the ADC values on screen. Usage `wibeth_binary_frame_reader <input_file_name>`.  

* `wibeth_binary_frame_modifier` is used to create a custom WIBEth frame file suitable for testing different patterns. The application will produce an output file `wibeth_output.bin`. There are no command line options, please refer to the code for further details (e.g. what ADC value to set, which time frame to use, etc.). 

* `streamed_TPs_to_text` tool used to streamed TPs from an HDF5 file (e.g. TPSTream recording) into the text format. Example of usage: 
```sh
streamed_TPs_to_text -i INPUT_TPSTREAM.hdf5  -o OUTPUT.txt
```

* `check_fragment_TPs` tool for reading TP fragments from file and check that they have start times within the request window

## Python utility tools for TP studies

### Python libraries
In the `python` directory, there are some libraries with tools to perform simple analyses on TPs, convert between formats, create images. 
Plotting functions make use of `matplotlib`; to know how to load it in your environment, see the section below.

* `hdf5_converter.py` is a library to read the hdf5 TPstream, using directly the DAQ hdf5 tools.

* `utils.py` is a library with some utility functions, like the one to save TPs into a numpy array (`save_tps_array`).

* `properties_plotter.py` contains some functions to plot the properties of the TPs, like the start time, time over threshold, etc. 
They have several parameters to customize the plots, in particular in some cases the range is handled through the quantile of the distribution. 
Default value is 1, but to get rid of outliers it can be set to 0.9, for example.

* `cluster_maker.py` is a library to cluster TPs basing on channel and/or time. 
Each cluster should then be a track, if the conditions are set to be strict enough (deafult is channel limit 1 and time limit 3 ticks).

* `image_creator.py` contains functions to create images from TPs, that will be 2D histograms (channel vs time). 
Also in this case, many parameters can be set to customize the images, for example to choose if to fix the size or if it will depend on the cluster case by case.


### Python scripts
In the `scripts` directory, there are some example scripts that make useof the python libraries.
To avoid unwanted files in the repo, I added `*.png` to the `.gitignore` file.
Taking these as inspiration, you can create your own scripts to perform the analyses you need.
They have a verbose option (`-v`) that prints out the TPs that are being processed.

#### `plot_tp_properties.py` 

Script to plot all the properties of the TPs: start time, time over threshold, channel, detid, peak, peak time, integral. 
It can accept one or more files are input, both in hdf5 or text format. 
The flag `--superimpose` can be used to plot all the input files on the same plot. 
Default output path is `./`, output files will be named `input_file_timeStart.png` and so on.
You can see all options with 
```sh
python plot_tp_properties.py --help
```

Here an example of usage. `--all` is used to plot all the properties, but you can also choose to plot only one or some of them.
```sh
python plot_tp_properties.py -f INPUT_TPSTREAM.hdf5 INPUT_TPS.txt -n 1000 -e my/output/folder/ --superimpose --all --view U V X --channel-map CRP
python3.10 plot_tp_properties.py -f INPUT_TPSTREAM.hdf5 INPUT_TPS.txt  -e ./ --superimpose --all  --channel-map FiftyLChannelMap
python3.10 plot_tp_properties.py -f INPUT_TPSTREAM.hdf5 INPUT_TPS.txt  -e ./ --superimpose --all --view U V X --channel-map FiftyLChannelMap
```


#### `create_images.py` 
Script to create images from TP clusters, ie one or more tracks. 
It accepts one input file (hdf5 or text) and will produce images with increasing numbering in the output directory (default is `./`).
Default detector is "APA", but it can be changed with the flag `--channel-map` to "CRP" or "50L".
The images are by default `png`, named `<view>_track<number>.png`, e. g. `u_track12.png`.
The clusters can be saved to a text file to take a look at them in a readable format, with the flag `--save-clusters`.

You can see full usage with `--help`, here an example:
```sh
python create_images.py -i INPUT_TPSTREAM.hdf5 -n 1000 -o my/output/folder/ --ticks-limit 5 --channel-limit 2 --min-tps 3 
```


#### Setup DAQ environment on lxplus or NP04 machines (e.g. `np04-srv-019`)
To use the tools and scripts in this repository, the DUNE-DAQ software environment must be setup. The following commands are valid for lxplus machines and NP04 machines (e.g. `np04-srv-019`). 
```sh
source /cvmfs/dunedaq.opensciencegrid.org/setup_dunedaq.sh
setup_dbt latest
dbt-setup-release fddaq-v4.2.0
```
The version of the DUNE-DAQ software can also be a different one than `v4.2.0`.


#### Setup matplotlib on lxplus or NP04 machines (e.g. `np04-srv-019`)
To use the `matplotlib` python module run the following command on a console where the DUNE-DAQ software area has not been sourced:
```sh
export PREFIX_PATH=/your/chosen/path/matplotlib/
pip install --prefix=$PREFIX_PATH matplotlib
export PYTHONPATH=$PREFIX_PATH/lib/python3.10/site-packages/:$PYTHONPATH
```

## Validation tools

### Pattern generator

`wibeth_tpg_pattern_generator` is a pattern generator application for validating TPG algorithms. It takes as input an existing WIBEth frame file (`*.bin` file) and it will (a) set ADC values for a given channel and set of time ticks according to the desired pattern, and (b) reset all other ADC values to 0. 
To use the tool run the following:
```sh
$ wibeth_tpg_pattern_generator --help
TPG pattern generator using as input and storing output as raw WIBEth ADC binary files
Usage: wibeth_tpg_pattern_generator [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--file-path-input TEXT   Path to the input file
  -o,--path-output TEXT       Path to the output directory. Default: . (i.e. pwd)
  -n,--num-frames-to-read INT Number of WIBEth frames to read. Default: select all frames.
  -i,--input_channel UINT     Input channel number for adding fake hit. Default: 0. (max:63)
  -t,--tpg-threshold INT      Value of the TPG threshold. Default value is 500.
  --save-adc-data             Save ADC data (first frame only)
  --save-trigprim             Save trigger primitive data
  -c,--clock-tick-offset INT  Time tick of pattern start. Default: 1 (max:63).
  -s,--out_suffix TEXT        Append string (suffix) to output hit file name
  -p,--select-pattern TEXT    Test pattern name (patt_golden, patt_pulse, patt_square, patt_square_left, patt_square_right). Default: patt_golden.
  -w,--overwrite-wibeth-header
                              Overwrite crate, slot, stream IDs (needed for offline channel map). Default: false.
  -v,--verbose                Printout additional information while the application is running. Default: false.
```

Example of usage:
```sh
rm patt_golden_1_wibeth_output.bin patt_golden_1_wibeth_output_pedsub.bin patt_golden_1_wibeth_output_pedsub_hits.txt
wibeth_tpg_pattern_generator -f wibeth_output_all_zeros.bin -n 2 -i 0 -t 499 -o 1 -p patt_golden --save-trigprim -s __1 -v -w 
wibeth_tpg_pattern_generator -f wibeth_output_all_zeros.bin -o . --save-trigprim -w -n 2 -t 64 -i 0 -c 63 -p patt_golden -s __63
```

### Workload emulator

`wibeth_tpg_workload_emulator` is an emulator application for TPG algorithms. It takes as input a WIBEth frame file (`*.bin` file) and it will execute the selected TPG algorithm on the WIBEth frame data. This application is used for comparing NAIVE and AVX implementations of a given TPG algorithm using a well known input test pattern file. 

To use the tool run the following:
```sh
$ wibeth_tpg_workload_emulator --help
Generate workload for TPG algorithms using binary frame files.
Usage: wibeth_tpg_workload_emulator [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--file-path-input TEXT   Path to the input file
  -o,--path-output TEXT       Path to the output directory. Default: .
  -a,--algorithm TEXT         TPG Algorithm (SimpleThreshold / AbsRS). Default: SimpleThreshold
  -i,--implementation TEXT    TPG implementation (AVX / NAIVE). Default: AVX
  -m,--channel-map TEXT       Select a valid channel map: None, VDColdboxChannelMap, ProtoDUNESP1ChannelMap, PD2HDChannelMap, HDColdboxChannelMap, FiftyLChannelMap
  -d,--duration-test INT      Duration (in seconds) to run the test. Default value is 120.
  -t,--tpg-threshold INT      Value of the TPG threshold. Default value is 500.
  -c,--core INT               Set core number of the executing TPG thread. Default value is 0.
  --save-adc-data             Save ADC data (first frame only)
  --save-trigprim             Save trigger primitive data
  -r,--repeat-timer BOOLEAN   Repeat frame processing until certain time elapsed (true/false). Default: true.
  -s,--out-suffix TEXT        Append string to output hit file name (e.g. __1). Default: empty string).
  -n,--num-frames-to-read INT Number of frames to read. Default: -1 (select all frames).
```

Example of usage: 
```sh
$ wibeth_tpg_pattern_generator -f /cvmfs/dunedaq.opensciencegrid.org/assets/files/d/d/1/wibeth_output_all_zeros.bin -o . --save-trigprim -w -n 2 -t 64 -i 0 -c 63 -p patt_golden -s __63
$ tpg_workload_emulator -f patt_golden_chan_0_tick_63_wibeth_output.bin -r false -a SimpleThreshold -i NAIVE  -n 2 -t 499 -m ProtoDUNESP1ChannelMap
```

Please note, when using `wibeth_output_all_zeros.bin` input file from the asset repository, the `-w` option is needed to overwrite the header information. The generated pattern file, `patt_golden_chan_0_tick_63_wibeth_output.bin`, is then used as input to the `tpg_workload_emulator` app.  

More examples of usage:
```sh
$ wibeth_tpg_workload_emulator -f patt_golden_chan_0_tick_1_wibeth_output.bin -r false -a SimpleThreshold -i NAIVE -n 2 -t 64  --save-trigprim -s __1
$ wibeth_tpg_workload_emulator -f patt_golden_chan_0_tick_1_wibeth_output.bin -r false -a SimpleThreshold -i AVX -n 2 -t 64  --save-trigprim -s __1
wibeth_tpg_workload_emulator -f patt_golden_chan_0_tick_63_wibeth_output.bin -r false -m VDColdboxChannelMap --save-trigprim -n 2 -t 64 -c 63 -a AbsRS -i AVX
```

### Running `pytest`

The pytest framework is used to streamline the validation process and to create an automated workflow using the validation tools described in this section. 

One of the standard type of tests performed during TPG algorithm development is the comparison of the NAIVE and AVX implementation of the algorithm under study. However, existing algorithms should also be tested regularly to ensure readout and other code changes do not break the TPG functionality. 

The integration of `pytest` in the tpgtools repositor has the following directory structure:
```sh
tpgtools/
pytest.ini
conftest.py
tests/
    __init__.py
    test_tpg_implementation.py
scripts/
    test_tpg_implementation_bundle.sh
```

The tested parameters are provided on the command line in the following order (separated with a dash): 
```sh
pattern_name-channel_number-clock_tick-threshold-algorithm e.g. patt_golden-0-1-64-2-SimpleThreshold
where
    pattern_name   = patt_golden, patt_square, patt_square_left, patt_square_right
    channel_number = [0, 63]
    clock_tick     = [1, 63]
    threshold      = 10, 64, 499
    algorithm      = SimpleThreshold, AbsRS
```

Regarding the patterns, For example, the `patt_square` which creates a constant-ADC pattern, few ticks wide, on a given channel, has a fixed amplitude of 500. For `patt_golden`, the sequence of ADC values is: 500, 502, 504, 505, 506, 505, 504, 502, 500.


To run the full suite of tests use:
```sh
addopts = --quiet -ra  (in pytest.ini)
pytest -k "test_tuple_params" tests
```
As the `pytest` script currently scans/tests all parameters sequentially, it takes too long to complete the test job. Therefore the parameters are provided explicitly in order to run only a subset of all tests. The examples below show how to start a test from the command line or to run the default subset of tests.


#### Example to run a single test
```sh
cd sourcecode/tpgtools
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
pytest -k "test_all_params[patt_golden-0-1-64-2-SimpleThreshold]" tests
pytest --printout=more -k "test_all_params[patt_golden-0-1-64-2-AbsRS]" tests
```

This test gives immediate feedback as to whether the test passed or failed. The output from the test is stored in `/tmp/pytest-of-<user>` on the machine where the test was run. 

#### Example to run the default test bundle 
```sh
cd sourcecode/tpgtools/scripts
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
./test_tpg_implementation_bundle.sh -d | tee pytest_log.txt
```

This test takes a about 15 minutes to complete and provides a reasonable coverage of the parameters. In a next version we aim to reduce the time and provide a one-line summary of all tests. 

## Notes
- The tools and scripts developed have been used for TPG related activities. They have not been generalized to cover all use-cases. If there is a need or feature request, ask mainteners of the repository.  
- The repository also contains tools for WIB2 frames but they are not kept up to date. Please refer to the code or ask mainteners of the repository for help. 

