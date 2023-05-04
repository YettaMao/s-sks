# This script is to calculate the Observed and Calculated time difference of SKS
# maoyt 20230426 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from math import *
import sys

syn_file = 'epi83/u.xy'
obs_file = '../ncisp6/data/raw/2007_289_21_05_41/obs_r.xy'
sta_file = 'epi83/sta.xy'
thres_syn = 0.003
thres_obs = 0.0015
prem_iaspi = 2.3
vel = 4.48
h = 50

def calculate_arrival_time(data,threshold):
    data= np.array(data)
    diff = np.abs(np.diff(data[:,1]))
    arrival_index = np.argmax(diff>threshold)
    return data[arrival_index+1,0]

#def read_obs(): #read Obs record

def read_data(filename): #read seismic record
    data_arrays = []
    with open(filename) as f:
        for line in f:
            if line.startswith(">>"):
                current_array = []
            else:
                time, amplitude = map(float, line.split())
                if -10 <= time <=10:
                    current_array.append([time,amplitude])
            if not line.strip() or line.startswith(">>"):
                data_arrays.append(current_array)
    return data_arrays

def plot_record(data_arrays):
    fix, ax = plt.subplots()
    for i, data_array in enumerate(data_arrays):
        times = [row[0] for row in data_array]
        amplitudes = [row[1] for row in data_array]
        ax.plot(times, amplitudes, label=f"Series {i+1}")
    ax.set_xlabel("Time")
    ax.set_ylabel("GCARC")
    plt.show()


### read synthetic results and calculate the arrival time
record_syns = read_data(syn_file)
arrival_times_syn = [calculate_arrival_time(record_syn,thres_syn)+prem_iaspi for record_syn in record_syns]
# for i,arrival_time_syn in enumerate(arrival_times_syn): # print results
    # print(f"Arrival time of syn_record {i+1}: {arrival_time_syn: .3f}s")

### read synthetic results and calculate the arrival time
record_obss = read_data(obs_file)
arrival_times_obs = [calculate_arrival_time(record_obs,thres_obs) for record_obs in record_obss]
# for i,arrival_time_obs in enumerate(arrival_times_obs):  # print results
    # print(f"Arrival time of obs_record {i+1}: {arrival_time_obs: .3f}s")

### plot results
a_syn = np.array(arrival_times_syn)
a_obs = np.array(arrival_times_obs)
a_reduce = a_obs - a_syn
red_vel = h/(h + vel* a_reduce)-1
sta = np.array(pd.read_csv(sta_file,sep='\s+',header = None))[:,0]
fig, ax = plt.subplots(figsize=(25,6))
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.scatter(sta,red_vel, s=100, c='#88c999')
plt.axhline(0, linewidth = 3.0, linestyle = "dashed", color = 'darkgray')
plt.xlim((81.935,91.626))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
ax.set_xlabel("GCARC ($\circ)$",fontsize = 20)
ax.set_ylabel("Vs Pertubation ", fontsize = 20)

cmap = plt.cm.get_cmap('Blues')
angle_range = [0, 90]
norm = plt.Normalize(*angle_range)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm,alpha=0.4,label='Orientation', pad=0.03)
cbar.ax.tick_params(labelsize =18)
cbar.set_label('Orientation',size=18,rotation=90, labelpad=10)
plt.axvspan(81.935,84, color = cmap(norm(10)), alpha = 0.4)
plt.axvspan(84,85, color = cmap(norm(20)), alpha = 0.4)
plt.axvspan(85,86, color = cmap(norm(5)), alpha = 0.4)
plt.axvspan(86,86.3, color = cmap(norm(0)), alpha = 0.4)
plt.axvspan(86.3,87, color = cmap(norm(180-110)), alpha = 0.4)
plt.axvspan(87,87.3, color = cmap(norm(180-170)), alpha = 0.4)
plt.axvspan(87.3,87.8, color = cmap(norm(20)), alpha = 0.4)
plt.axvspan(87.8,88.5, color = cmap(norm(20)), alpha = 0.4)
plt.axvspan(88.5,89, color = cmap(norm(180-100)), alpha = 0.4)
plt.axvspan(89,89.9, color = cmap(norm(180-160)), alpha = 0.4)
plt.axvspan(89.9,90.8, color = cmap(norm(5)), alpha = 0.4)
plt.axvspan(90.8,91.626, color = cmap(norm(20)), alpha = 0.4)

plt.scatter(sta,red_vel, s=100, c='#FF8884')

plt.savefig('Vs_pertubation.pdf')
plt.show()


# plot_record(record_syns)
# plot_record(record_obss)
