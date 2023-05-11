import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from math import *
import sys


### plot topography
topo = 'reduce_vel/topo_profile.txt'
data = np.loadtxt(topo)
fix,ax = plt.subplots(figsize=(25,2))
plt.plot(data[:,0],data[:,3],color='lightgray')
ax.fill_between(data[:,0],data[:,3],color='lightgray')

plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.xlim(0,10)
ax.set_ylabel("Elevation (m) ", fontsize = 20)
ax.spines['bottom'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)


plt.savefig('topo_profile.pdf')
plt.show()