# This script is to calculate the Observed and Calculated time difference of SKS
# maoyt 20230426 

import numpy as np
import matplotlib.pyplot as plt
from math import *
import sys

syn_file = 'epi83/u.xy'

def calculate_arrival_time(data,threshold):
    diff = np.abs(np.diff(data))
    arrival_index = np.where(diff>threshold)[0][0]
    return arrival_index + 1

#def read_obs(): #read Obs record

def read_syn(filename): #read Syn record
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

def plot_syn(data_arrays):
    fix, ax = plt.subplots()
    for i, data_array in enumerate(data_arrays):
        times = [row[0] for row in data_array]
        amplitudes = [row[1] for row in data_array]
        ax.plot(times, amplitudes, label=f"Series {i+1}")
    ax.set_xlabel("Time")
    ax.set_ylabel("GCARC")
    plt.show()

record_syn = read_syn(syn_file)
plot_syn(record_syn)
# det_t = calculate_arrival_time(wt,0.01)
# print (det_t)