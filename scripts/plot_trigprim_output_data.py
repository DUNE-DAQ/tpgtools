from matplotlib import pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Trigger Primitives plotter")
parser.add_argument('-f' '--file', dest='file', help='Trigger Primitives')
args = parser.parse_args()
file = args.file

data = np.genfromtxt(file, delimiter=' ')
data = data.transpose()
fig = plt.subplot(111)
legend_properties = {'weight': 'bold'}

start_TPs=0
n=100 # Number of TPs to plot
t0 = data[0][:10]
startTime = data[0][start_TPs:n]
#time_shift = startTime[0]
time_shift = t0[0]
channel = data[3][start_TPs:n]
adc_sum = data[4][start_TPs:n]
startTime -= time_shift
startTime *= 16e-9

label="Input TPs - Event Display"
fig.set_xlabel("Relative Time (s)")
fig.set_ylabel("Channel ID")
#fig.set_ylim([120, 128])

#fig.set_title(title, fontweight='bold')
scatter = fig.scatter(startTime, channel, s=5, label=label, c=adc_sum, vmax=np.max(adc_sum)/20)
fig.legend(prop=legend_properties)


# Add color bar
colorbar = plt.colorbar(scatter)
colorbar.set_label('ADC Sum')

plt.show()
plt.savefig("output_trigger_primitives.png")
