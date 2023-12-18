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
```

The command line option `save_adc_data` allows to save the raw ADC values in a txt file after the 14-bit to 16-bit expansion. The command line option `save_trigprim`  allows to save the in a file the Trigger Primitive object information in a txt file. 

Example of usage: 
```sh
$ wibeth_tpg_algorithms_emulator -f swtest_run000035_0000_dataflow0_datawriter_0_20231102T083908.hdf5  -a SimpleThreshold -m PD2HDChannelMap -t 500 --save-trigprim --parse_trigger_primitive
$ wibeth_tpg_algorithms_emulator -f swtest_run000035_0000_dataflow0_datawriter_0_20231102T083908.hdf5  -a AbsRS -m PD2HDChannelMap -t 500 --save-adc-data  -n 5 
```


## Utility tools and scripts

* `wibeth_tpg_workload_emulator` is a simple emulator for TPG algorithm in either a naive or in AVX implementation. The application allows to emulate the workload when running a TPG algorithm and therefore monitor performance metrics. It requires an input binary frame file (check assets-list for valid input files) and it will execute the desired TPG algorithm for a configurable duration (default value is 120 seconds). The application is single threaded, pinned to core 0 (configurable). Check the helper page for more details. Example of usage: `wibeth_tpg_workload_emulator -f wibeth_frame_file.bin` 

* `wibeth_tpg_validation` is a simple emulator for validating different TPG algorithms, either in naive or in AVX implementation. Check the helper page for more details. Usage: `wibeth_tpg_validation -f wibeth_frame_file.bin` 

* `wibeth_binary_frame_reader`: reads a WIBEth frame file (`.bin` file) and prints all the ADC values on screen. Usage `wibeth_binary_frame_reader <input_file_name>`.  

* `wibeth_binary_frame_modifier` is used to create a custom WIBEth frame file suitable for testing different patterns. The application will produce an output file `wibeth_output.bin`. There are no command line options, please refer to the code for further details (e.g. what ADC value to set, which time frame to use, etc.). 

* `plot_trigprim_output_data.py` plots the Trigger Primitive output file obtained through `wibeth_tpg_algorithms_emulator` (when `save_trigprim` flag is enabled) and produces a plot called `output_trigger_primitives.png` . This script requires the use of `matplotlib`. To use the script run the following command: 
```sh
python3 plot_trigprim_output_data.py  -f TP_OUTPUT.TXT
```

* `streamed_TPs_to_text` tool used to streamed TPs from an HDF5 file (e.g. TPSTream recording) into the text format. Example of usage: 
```sh
streamed_TPs_to_text -i INPUT_TPSTREAM.hdf5  -o OUTPUT.txt
```

* `check_fragment_TPs` tool for reading TP fragments from file and check that they have start times within the request window


#### Setup matplotlib on NP04 machines (e.g. `np04-srv-019`)
To use the `matplotlib` python module run the following command on a console where the DUNE-DAQ software area has not been sourced:
```sh
export PREFIX_PATH=$HOME
pip install --prefix=$PREFIX_PATH matplotlib
export PYTHONPATH=$HOME/lib/python3.10/site-packages/:$PYTHONPATH
```


## Notes
- The tools and scripts developed have been used for TPG related activities. They have not been generalized to cover all use-cases. If there is a need or feature request, ask mainteners of the repository.  
- The repository also contains tools for WIB2 frames but they are not kept up to date. Please refer to the code or ask mainteners of the repository for help. 

