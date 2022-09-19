#!/bin/python
"""
Usage: python plot_2d_powder.py {config file}

Read data from 2d_pseudo_powder.dat, make plot.

"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

if (len(sys.argv) < 2) :
    print("Usage: plot_2d_powder.py [config file]")
    exit()

config_file = "{0:s}".format(sys.argv[1])
fp = open(config_file, "r")
lines = fp.readlines()
fp.close()

for line in lines:
    if "data_dir" in line:
        data_dir = line.split( "= ")[1].strip()
    if "num_row" in line:
        num_row = int(line.split(" = ")[1])
    if "num_col" in line:
        num_col = int(line.split(" = ")[1])
    if "num_raw_data" in line:
        num_data = int(line.split(" = ")[1])

dir2 = os.path.join(data_dir,"2d_pseudo_powder.dat")
fp = open(dir2, "r")
det = np.array(fp.readlines()).astype(int)
fp.close()

ll = len(det)
num_relevant_frames = det[ll-1]
det = np.resize (det, ll-1)
det = np.reshape(det, (num_row, num_col))
img = np.log10(det+1.0)
fig, ax = plt.subplots()
ax.set_title('Pseudo-powder diffraction pattern')
im = ax.imshow(img, interpolation='nearest', cmap=plt.get_cmap('jet'))
tmp = "Powder pattern from {0:d} images".format(num_relevant_frames)
plt.title(tmp)
plt.savefig('powder_2d.png', bbox_inches='tight')
plt.show()
