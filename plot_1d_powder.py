#!/bin/python
"""
Usage: python plot_1d_powder.py {config file}

Get name of unit cell file from config file, read it.
Read data from 1d_pseudo_powder.dat, make plot.

"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

if (len(sys.argv) < 2) :
    print("Usage: plot_1d_powder.py [config file]")
    exit()

config_file = "{0:s}".format(sys.argv[1])
fp = open(config_file, "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    entry_name = lines[t].split(" = ")[0]
    if ("data_dir" == entry_name):
        data_dir = lines[t].split( "= ")[1].strip()
    if ("unit_cell_file" == entry_name) :
        unit_cell_file = lines[t].split(" = ")[1].strip()

fin = open(unit_cell_file, "r")
lines = fin.readlines()
fin.close()
words = lines[0].split()
a = np.double(words[0])
b = np.double(words[1])
c = np.double(words[2])
alpha = np.double(words[3])
beta = np.double(words[4])
gamma = np.double(words[5])
print("Using unit cell",a,b,c,alpha,beta,gamma)

dir1 = os.path.join(data_dir,"1d_pseudo_powder.dat")
fin = open(dir1, "r")
lines = fin.readlines()
fin.close()

length = len(lines)
qval = np.zeros(length)
radial_hist = np.zeros((2, length))
for i in range(length):
    words = lines[i].split()
    qval[i] = np.double(words[0])
    radial_hist[0][i] = np.double(words[1])
    radial_hist[1][i] = np.double(words[2])

idx_max = 10
qmax = 0.075
epsilon = 1.e-10

# reciprocal basis vectors
conv = np.arctan(1.0)*4.0/180.0
alpharad = alpha*conv
betarad = beta*conv
gammarad = gamma*conv
sina = np.sin(alpharad)
cosa = np.cos(alpharad)
sinb = np.sin(betarad)
cosb = np.cos(betarad)
sing = np.sin(gammarad)
cosg = np.cos(gammarad)
cosas = (cosg*cosb-cosa)/ (sinb*sing)
sinas = np.sqrt(1.0-cosas*cosas)
cosbs = (cosa*cosg-cosb)/ (sina*sing)
sinbs = np.sqrt(1.0-cosbs*cosbs)
cosgs = (cosa*cosb-cosg)/ (sina*sinb)
sings = np.sqrt(1.0-cosgs*cosgs)

sumc = (alpharad+betarad+gammarad)*0.5
vol = 2.0*a*b*c*np.sqrt(np.sin(sumc-alpharad)*np.sin(sumc-betarad)*np.sin(sumc-gammarad)*np.sin(sumc));

recipa = b * c * sina/vol
recipb = c * a * sinb/vol
recipc = a * b * sing/vol
recipalrad = np.arctan2(sinas,cosas)
recipberad = np.arctan2(sinbs,cosbs)
recipgmrad = np.arctan2(sings,cosgs)

recip_mat = np.zeros((3,3))
recip_mat[0][0] = recipa
recip_mat[0][1] = recipb * cosgs
recip_mat[0][2] = recipc * cosbs
recip_mat[1][1] = recipb * sings
recip_mat[1][2] = -recipc * sinbs * cosa
recip_mat[2][2] = recipc * sinbs * sina

vec_a = np.array([recip_mat[0][0], recip_mat[1][0], recip_mat[2][0]])
vec_b = np.array([recip_mat[0][1], recip_mat[1][1], recip_mat[2][1]])
vec_c = np.array([recip_mat[0][2], recip_mat[1][2], recip_mat[2][2]])
rvec = np.zeros(3)

for p in range(2):

    if (p == 0):
        print("1D histogram of inter-peak distances in reciprocal space")
    else:
        print("1D histogram of spatial frequency magnitudes of peaks")

    ymax = np.max(radial_hist[p])*1.1
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.set_title('Azimuthally averaged peaks')
    plt.plot(qval, radial_hist[p], 'b')
    if (p == 0):
        plt.title('Inter-peak distances in recip space, dashed red lines predicted from cell')
    else:
        plt.title('Magnitude of recip space vectors, dashed red lines predicted from cell')
    ystep = int(ymax*0.8)
    if (ystep == 0) :
        ystep = 1
    
    for h in range(-idx_max, idx_max+1):
        for k in range(-idx_max, idx_max+1):
            for l in range(-idx_max, idx_max+1):
                q = 0.
                for i in range(3):
                    rvec[i] = h*vec_a[i] + k*vec_b[i] + l*vec_c[i]
                    q += rvec[i]*rvec[i]
                q = np.sqrt(q)
                if (q > qmax):
                    continue
    
                y = np.arange(0, ymax, ystep)
                x = np.ones_like(y)*q
                plt.plot(x, y, 'r--')
                
    plt.xlim([0, qmax])
    plt.ylim([0, ymax])
    plt.show()
