"""
 * @file plot_tp_properties.py
 * @brief Plots TP properties from a file or a list of files
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2024.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import datetime

import sys
sys.path.append("../python/")
from properties_plotter import *
from utils import create_channel_map_array

# parse from command line the args.files and the number of tps to plot
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files",                 nargs="+",              help="files to read")
parser.add_argument("-e", "--output-folder",         type=str,               help="output folder, remember to put / at the end",   default="./")
parser.add_argument("-n", "--number-tps",            type=int,               help="number of tps to plot", default=100)
parser.add_argument("-s", "--superimpose",           action="store_true",    help="superimpose plots")
parser.add_argument("--view",                        nargs="+",              help="default is allViews on same plot, to have separate select one or more among U, V, X", default=["allViews"])
parser.add_argument("--channel-map",                 type=str,               help="channel map to use. Options are APA, CRP, FiftyLChannelMap", default="APA") # thinking of NP04 as default
parser.add_argument("-a", "--all",                   action="store_true",    help="plot all the TP variables",default=False)
parser.add_argument("-t", "--time-peak",             action="store_true",    help="plot time peak")
parser.add_argument("-o", "--time-over-threshold",   action="store_true",    help="plot time over threshold")
parser.add_argument("-i", "--adc-integral",          action="store_true",    help="plot adc integral")
parser.add_argument("-p", "--adc-peak",              action="store_true",    help="plot adc peak")
parser.add_argument("-c", "--channel",               action="store_true",    help="plot channel")
parser.add_argument("-d", "--detid",                 action="store_true",    help="plot detid")
parser.add_argument("--show",                        action="store_true",    help="show plots", default=False)
parser.add_argument("-v", "--verbose",               action="store_true",    help="verbose mode", default=False)
args = parser.parse_args()

if args.all:
    args.time_peak = True
    args.time_over_threshold = True
    args.adc_integral = True
    args.adc_peak = True
    args.channel = True
    args.detid = True

# # print recap of options
# print("#############################################")
# print("Selected options:")
# print(" - Files: " + str(args.files))
# print(" - Output folder: " + str(args.output_folder))
# print(" - Number of TPs: " + str(args.number_tps))
# print(" - Superimpose: " + str(args.superimpose))
# print(" - Views: " + str(args.view))
# print(" - Channel map: " + str(args.channel_map))
# print(" - All: " + str(args.all))
# print(" - Time peak: " + str(args.time_peak))
# print(" - Time over threshold: " + str(args.time_over_threshold))
# print(" - ADC integral: " + str(args.adc_integral))
# print(" - ADC peak: " + str(args.adc_peak))
# print(" - Channel: " + str(args.channel))
# print(" - Detid: " + str(args.detid))
# print(" - Show: " + str(args.show))
# print(" - Verbose: " + str(args.verbose))
# print("#############################################")

# save the settings printed before in an array of strings 
script_settings = []
script_settings.append(" - Files: " + str(args.files))
script_settings.append(" - Output folder: " + str(args.output_folder))
script_settings.append(" - Number of TPs: " + str(args.number_tps))
script_settings.append(" - Superimpose: " + str(args.superimpose))
script_settings.append(" - Views: " + str(args.view))
script_settings.append(" - Channel map: " + str(args.channel_map))
script_settings.append(" - All: " + str(args.all))
script_settings.append(" - Time peak: " + str(args.time_peak))
script_settings.append(" - Time over threshold: " + str(args.time_over_threshold))
script_settings.append(" - ADC integral: " + str(args.adc_integral))
script_settings.append(" - ADC peak: " + str(args.adc_peak))
script_settings.append(" - Channel: " + str(args.channel))
script_settings.append(" - Detid: " + str(args.detid))
script_settings.append(" - Show: " + str(args.show))
script_settings.append(" - Verbose: " + str(args.verbose))

print("Selected options:")
for i in range(len(script_settings)):
    print (script_settings[i])

# read the file(s) and create arrays of TPs, using the class in TriggerPrimitive.py
# the code will deduce if the file is coming from offline or online data

tps_lists = []
for tpFile_path in args.files:

    print("Reading file: " + tpFile_path)

    this_tps_list = save_tps_array(tpFile_path, args.number_tps)
    
    if args.verbose:
        print("TPs list: ", this_tps_list)
    
    tps_lists.append(this_tps_list)

for i in range(len(tps_lists)):
    print ("Number of TPs in list ", i, " : ", len(tps_lists[i]))

print (tps_lists[0][0]['channel'])

# creating a channel map for the given detector, to know which channel is in which view    
this_channel_map = create_channel_map_array(args.channel_map)

for view in args.view:
    if view == "allViews":
        print ("Plotting all views on the same plot:")
        # select all channels from map
        channels_to_plot = this_channel_map[:, 0]
    elif view == "U":
        print ("Plotting U view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the U view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 0, 0]  
    elif view == "V":
        print ("Plotting V view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the V view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 1, 0]
    elif view == "X":
        print ("Plotting X view:")
        # the first element of channel map is the channel number, the second is the view
        # we want to plot all the TPs in the X view
        channels_to_plot = this_channel_map[this_channel_map[:, 1] == 2, 0]
    else :
        print ("View not recognized, please choose from U, V, X or allViews (default)")
        sys.exit(1)  
    
    print ("Number of channels to plot: ", len(channels_to_plot))
    
    if args.verbose:
        print ("Plotting channels: ", channels_to_plot)
    
    # Select TPs to plot basing on the selected view and channel map
    # create a new list of lists of TPs, keeping only the TPs that are in the selected view
    # and in the selected channels
    tps_to_plot = []
    for tps_list in tps_lists:
        if args.verbose:
            print ("Number of TPs in this list: ", len(tps_list))
        this_tps_to_plot = []
        for tp in tps_list:
            if tp['channel']%(len(this_channel_map[0])) in channels_to_plot:
                this_tps_to_plot.append(tp)
        print ("Number of TPs to plot in this file: ", len(this_tps_to_plot))
        tps_to_plot.append(this_tps_to_plot)
        
    
    print("Plotting time peak...")
    if args.time_peak:
        time_peak_plot_path = plotTimePeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
    
    print("Plotting time over threshold...")
    if args.time_over_threshold:
        plotTimeOverThreshold(tps_to_plot,
                            args.files,
                            superimpose=args.superimpose,
                            quantile=1,
                            y_min=None, y_max=None,
                            output_folder=args.output_folder,
                            output_suffix=view,
                            show=args.show)

    print("Plotting ADC integral...")
    if args.adc_integral:
        plotADCIntegral(tps_to_plot,
                        args.files,
                        superimpose=args.superimpose,
                        quantile=1,
                        y_min=None, y_max=None,
                        output_folder=args.output_folder,
                        output_suffix=view,
                        show=args.show) 

    print("Plotting ADC peak...")
    if args.adc_peak:
        plotADCPeak(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    quantile=1,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print("Plotting channel...")
    if args.channel:
        plotChannel(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    x_min=None, x_max=None,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)
                
    print ("Plotting detid...")
    if args.detid:
        plotDetId(tps_to_plot,
                    args.files,
                    superimpose=args.superimpose,
                    y_min=None, y_max=None,
                    output_folder=args.output_folder,
                    output_suffix=view,
                    show=args.show)

    print(" ")
    print ("Finished plotting! Products are in: " + args.output_folder + " and are named after the input file and variables plotted, plus the suffix (if passed).")

   
### Create a pdf report

# open pdf file


pdf_report_name = '{}/tp_properties_report.pdf'.format(args.output_folder)

with PdfPages(pdf_report_name) as pdf:
    # Example:
    plt.figure(figsize=(8, 6))
    plt.imshow(plt.imread('../tools/dune_logo.jpg'))
    plt.axis('off')  # Hide axis
    # plt.title('DUNE Logo')
    pdf.savefig()  # Save figure to PDF
    plt.close()

    # Add title of the report
    pdf_report_name = pdf_report_name.replace(args.output_folder, '')  # Remove folder path
    # pdf.subtitle('Performance Report: {}'.format(pdf_report_name))
    
    # pdf.set_font('Times', 'B', 16)
    plt.text(0.1, 0.9, 'Selected options:')
    plt.text(0.1, 0.9, script_settings, ha='left', va='top')
    pdf.savefig()  # Save figure to PDF
    plt.close()

    # Add plots. Add time_peak_plot_path
    plt.imshow(plt.imread(time_peak_plot_path))
    pdf.savefig()  # Save figure to PDF
    # clean page
    plt.close()

    # Add time of creation
    now = datetime.datetime.now()
    creation_time_text = 'Report created on: {}'.format(now.strftime("%Y-%m-%d %H:%M:%S"))
    plt.text(0.1, 0.9, creation_time_text, ha='left', va='top')
    pdf.savefig()  # Save figure to PDF

print("PDF report created: {}".format(pdf_report_name))


# pdf.ln(1)
# pdf.image('../tools/dune_logo.jpg', w=180)
# pdf.ln(2)
# pdf.set_font('Times', 'B', 16)
# pdf.cell(40,10,'TP Properties Report')
# pdf.ln(10)

# # script_settings
# pdf.set_font('Times', '', 12)
# pdf.write(5, 'Selected options:')
# pdf.ln(5)
# for i in range(len(script_settings)):
#     pdf.write(5, script_settings[i])
#     pdf.ln(5)
    
# pdf.ln(10)

# # add image 
# pdf.image(time_peak_plot_path, w=180)

# pdf.write(5, 'Report made on {}'.format(current_time()))

# print('The report was created and saved to {}'.format(pdf_name))


# create a pdf report with the plots of the TPs
# def create_report_performance(input_dir, output_dir, daqconfs_cpupins_folder_parent_dir, process_pcm_files=False, process_uprof_files=False, print_info=True, streams='8, 16, 24, 32, 40, and 48', pdf_name='performance_report', repin_threads_file=None, comment=['TBA']):    
# directory([input_dir, output_dir])
# pcm_file, uprof_file, core_utilization_file, reformated_uprof_file, reformated_core_utilization_file, all_file, all_plots_file = make_name_list(input_dir)



# # Processing the data first
# if process_pcm_files:
#     for i, file_pcm_i in enumerate(pcm_file):
#         add_new_time_format(input_dir, file_pcm_i)

# if process_uprof_files:
#     for i, file_uprof_i in enumerate(uprof_file):
#         uprof_pcm_formatter(input_dir, file_uprof_i)
#         add_new_time_format(input_dir, 'reformatter_{}'.format(file_uprof_i))

# cpupins_utilazation_reformatter(input_dir)
# for i, file_core_i in enumerate(core_utilization_file):
#         add_new_time_format_utilization(input_dir, 'reformatter_{}'.format(file_core_i))
    
# if process_pcm_files or process_uprof_files:
#     print('Finish the processing of the data.')

# print(all_file[0])
# info_pcm_basic = break_file_name(all_file[0])

# # creating report
# pdf.set_font('Times', '', 10)
# pdf.write(5, 'The tests were run for the WIB{} data format for {} streams. The Figures 1 and 2 show the results of the tests ran (Table1) using the different metrics. \n'.format(info_pcm_basic[4], streams))
# pdf.write(5, '    * L2-hits is the fraction of requests that make it to L2 at all. Similar for L3. \n')
# pdf.write(5, '    * L2-misses is the fraction of requests that make it to L2 at all and then miss in L2. Similar for L3. \n')
# pdf.ln(10)

# #-------------------------------------------TABLE-----------------------------------------------
# # Data to tabular
# rows_data = []
# headers = ['Test', 'Readout SRV', 'dunedaq', 'OS', 'NODE', 'General comments']
# rows_data.append(headers)

# line_height = pdf.font_size * 2
# col_width = [pdf.epw/3.8, pdf.epw/8, pdf.epw/7, pdf.epw/12, pdf.epw/12, pdf.epw/5]  
# lh_list = [] #list with proper line_height for each row

# for i, file_i in enumerate(all_file):
#     info = break_file_name(file_i)
#     test_info = re.sub('_', ' ', info[5])
#     line = [test_info, info[2], info[1], check_OS(info[2]), info[3], comment[i]]
#     rows_data.append(line)

# # Determine line heights based on the number of words in each cell
# for row in rows_data:
#     max_lines = 1  # Initialize with a minimum of 1 line
#     for datum in row:
#         lines_needed = len(str(datum).split('\n'))  # Count the number of lines
#         max_lines = max(max_lines, lines_needed)

#     lh_list.append(line_height * max_lines)
    
# # Add table rows with word wrapping and dynamic line heights
# for j, row in enumerate(rows_data):
#     line_height_table = lh_list[j] 
#     for k, datum in enumerate(row):
#         pdf.multi_cell(col_width[k], line_height_table, datum, border=1, align='L', new_x=XPos.RIGHT, new_y=YPos.TOP, max_line_height=pdf.font_size)
        
#     pdf.ln(line_height_table)
    
# pdf.write(5, 'Table 1. Summary of the tests ran. \n')    
# pdf.ln(10)

# #--------------------------------------------FIGURES------------------------------------------------
# plot_vars_comparison(input_dir, output_dir, pdf_name)

# if info[3] == '0' or info[3] == '01':
#     pdf.image('{}/{}_results_{}_socket0.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#     pdf.write(5, 'Figure 1. Socket0 results of the tests ran using the metrics CPU Utilization (%), Memory Bandwidth (GB/sec), Instructions Per Cycle, Instructions Retired Any (Million).')
#     pdf.ln(10)
#     pdf.image('{}/{}_results_cache_{}_socket0.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#     pdf.write(5, 'Figure 2. Socket0 results of the tests ran using the metrics L2 Cache Misses (Million), L2 Cache [Misses/Hits] (%), L3 Cache Misses (Million), and L3 Cache [Misses/Hits] (%).')
#     pdf.ln(10)
    
#     if info[3] == '01':
#         pdf.image('{}/{}_results_{}_socket1.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#         pdf.write(5, 'Figure 3. Socket1 results of the tests ran using the metrics CPU Utilization (%), Memory Bandwidth (GB/sec), Instructions Per Cycle, Instructions Retired Any (Million).')
#         pdf.ln(10)
#         pdf.image('{}/{}_results_cache_{}_socket1.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#         pdf.write(5, 'Figure 4. Socket1 results of the tests ran using the metrics L2 Cache Misses (Million), L2 Cache [Misses/Hits] (%), L3 Cache Misses (Million), and L3 Cache [Misses/Hits] (%).')
#         pdf.ln(10)
    
# if info[3] == '1':
#     pdf.image('{}/{}_results_{}_socket1.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#     pdf.write(5, 'Figure 1. Socket1 results of the tests ran using the metrics CPU Utilization (%), Memory Bandwidth (GB/sec), Instructions Per Cycle, Instructions Retired Any (Million).')
#     pdf.ln(10)
#     pdf.image('{}/{}_results_cache_{}_socket1.png'.format(output_dir, pdf_name, info_pcm_basic[4]), w=180)
#     pdf.write(5, 'Figure 2. Socket1 results of the tests ran using the metrics L2 Cache Misses (Million), L2 Cache [Misses/Hits] (%), L3 Cache Misses (Million), and L3 Cache [Misses/Hits] (%).')
#     pdf.ln(10)
    
# #----------------------------------------CONFIGURATIONS---------------------------------------------
# pcm_file, uprof_file, core_utilization_file, reformated_uprof_file, reformated_core_utilization_file, all_file, all_plots_file = make_name_list(input_dir)

# if print_info:
#     pdf.write(5, 'Configurations: \n', 'B')
#     for i in range(len(all_file)):
#         info = break_file_name(all_file[i])
        
#         var_i='ru{}{}{}'.format(info[2], info[4], '0')
#         file_daqconf_i='daqconf-{}-{}-{}-{}'.format(info[4], info[5], info[2], info[3])
#         file_core_i='reformatter_core_utilization-{}-{}-{}-{}-{}'.format(info[1], info[2], info[3], info[4], info[5])
        
#         if os.path.exists('{}/{}'.format(input_dir,file_core_i)):
#             json_info(file_daqconf=file_daqconf_i, file_core=file_core_i, input_directory=daqconfs_cpupins_folder_parent_dir, input_dir=input_dir, var=var_i, pdf=pdf, if_pdf=print_info, repin_threads_file=repin_threads_file[i])
#         else:
#             print('Missing {}'.format(file_core_i))
#             json_info(file_daqconf=file_daqconf_i, file_core='reformatter_core_utilization-all0', input_directory=daqconfs_cpupins_folder_parent_dir, input_dir='{}/tools'.format(daqconfs_cpupins_folder_parent_dir), var=var_i, pdf=pdf, if_pdf=print_info, repin_threads_file=repin_threads_file[i])
#         pdf.cell(0, 10, 'Table {}. CPU core pins information for the "{}" test using dune_daq {}.'.format(i+2, info[5], info[1]))
#         pdf.ln(10)           
        
# pdf.ln(20)
# pdf.set_font('Times', '', 10)
# pdf.write(5, 'The End, made on {}'.format(current_time()))
# pdf.output('{}/{}_report.pdf'.format(output_dir, pdf_name))

# print('The report was create and saved to {}/{}.pdf'.format(output_dir, pdf_name))

