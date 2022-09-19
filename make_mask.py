#!/bin/python
"""
make detector mask

usage:
make_mask.py [config file] [path to a data frame in cbf or hdf5 format]

need:
config.ini, a data frame

makes:
mask file (usually mask.dat)

"""
from subprocess import *
import numpy as np
import matplotlib.pyplot as plt
import fabio
import h5py
import hdf5plugin
import sys
import os

if (len(sys.argv) < 1) :
    print( "Usage: make_mask.py [config file] [data frame (.cbf file or .h5 master file)]")
    exit()

#plugin_dir = os.path.join(os.path.dirname(bitshuffle.__file__),
#                'plugin')
#os.environ["HDF5_PLUGIN_PATH"] = plugin_dir

short_max = 2**15 - 1
bh_width = 110
n_rect = 0
msk_r_xmin = np.zeros(20, dtype=int)
msk_r_xmax = np.zeros(20, dtype=int)
msk_r_ymin = np.zeros(20, dtype=int)
msk_r_ymax = np.zeros(20, dtype=int)

gui = 0
if (len(sys.argv) > 1) :
    if (sys.argv[1] == "-gui") :
        config_file = "config.ini"
        gui = 1
    else :
        config_file = "{0:s}".format(sys.argv[1])
else :
    config_file = "config.ini"
fp = open(config_file, "r")
lines = fp.readlines()
fp.close()

print("run make_mask.py:\n")
for line in lines:
    words = line.split(" = ")
    if (words[0] == "data_dir"):
        data_dir = words[1].strip()
    if (words[0] == "raw_data_dir"):
        raw_data_dir = words[1].strip()
    if (words[0] == "required_string"):
        required_string = words[1].strip()
    if (words[0] == "num_row"):
        num_row = int(words[1].strip())
    if (words[0] == "num_col"):
        num_col = int(words[1].strip())
    if (words[0] == "cx"):
        c_row = np.double(words[1].strip())
    if (words[0] == "cy"):
        c_col = np.double(words[1].strip())
    if (words[0] == "bh_width"):
        bh_width = int(words[1].strip())
    if (words[0] == "mask_rect"):
        items = words[1].split(" ")
        msk_r_xmin[n_rect] = int(items[0].strip())
        msk_r_xmax[n_rect] = int(items[1].strip())
        msk_r_ymin[n_rect] = int(items[2].strip())
        msk_r_ymax[n_rect] = int(items[3].strip())
        n_rect += 1
    if (words[0] == "mask_file"):
        mask_file = words[1].strip()

if (len(sys.argv) > 2) :
    data_file = sys.argv[2]
else :
    cmd = "ls {0:s} | grep {1:s} | grep -v _master_ > TEMP_FILE".format( \
        raw_data_dir,required_string)
    p = Popen(cmd, shell=True)
    p.wait()
    fp = open("TEMP_FILE", "r")
    data_file = os.path.join (raw_data_dir, fp.readline().strip())
    fp.close()
    cmd = "rm TEMP_FILE"
    p = Popen(cmd, shell=True)
    p.wait()
print("data_file = {0:s}".format(data_file))
if (".cbf" in data_file) :
    print("Data type is CBF\n")
elif (".h5" in data_file) :
    print("Data type is HDF5\n")
else :
    print("Invalid data type\n")
    exit()
print("num_row = {0:d}, num_col = {1:d}, c_row = {2:.1f}, \
c_col = {3:.1f}\n".format(num_row, num_col, c_row, c_col))
print("half-width of beamstop holder = {0:d} pixels".format(bh_width))
if n_rect > 0:
    for k in range(n_rect):
        print("Mask rectangle, x from {0:d} to {1:d}, y from {2:d} to \
{3:d}".format(msk_r_xmin[k], msk_r_xmax[k], msk_r_ymin[k], \
         msk_r_ymax[k]))

if (".cbf" in data_file) :
    imgs = fabio.open(data_file)
    print(imgs)
    img = np.array(imgs.data)
elif (".h5" in data_file) :
    hdf5_file = h5py.File(data_file, 'r')
    root_group = hdf5_file.keys()
    dset = hdf5_file.get('/entry/data/data')
    dset_shape = dset.shape
    dsetp_shape = (2, dset_shape[1], dset_shape[2])
    h5_img = np.zeros(dsetp_shape,dtype=int)
    dim = len(dset_shape)
    dset.read_direct(h5_img, np.s_[1:2,:,:])
    hdf5_file.close()
    if (dim == 3) :
        img_shape = (dset_shape[1], dset_shape[2])
        img = np.zeros(img_shape, dtype=int)
        for i in range(num_row):
            for j in range(num_col):
                img[i][j] = h5_img[0][i][j]
    else :
        print("Unexpected dimensionality of data {0:d}, giving up.".format(dim))
        exit()

mask = np.zeros_like(img, dtype=int)
# mask detector gaps
mask[img < 0] = 1 
# mask hot pixels
mask[img > short_max] = 1

masky = np.array(mask,dtype=np.int64).flatten()
with h5py.File(mask_file,'w') as mask_file:
    #mask_set = mask_file.create_dataset("mask", data=masky)
    block_size = 0
    dataset = mask_file.create_dataset(
        "mask",
        data=masky,
        **hdf5plugin.Bitshuffle(nelems=0, lz4=True),
        dtype='int'
        )

img[mask == 1] = 0.
img = np.log10(img + 1.)
fig, ax = plt.subplots()
ax.set_title("Sample image file, after masking")
im = ax.imshow(img, interpolation='nearest', cmap=plt.get_cmap('jet'))
cb = fig.colorbar(im, ax=ax)
tmp = "File {0:s}".format(data_file)
plt.title(tmp)
if (gui == 1) :
    imfile = os.path.join (data_dir,"masked_det_image.png")
    plt.savefig(imfile, bbox_inches='tight')
else :
    plt.show()
