"""
Set up files and execute make_background

usage:
make_background.py [config file]

creates subdirectories radial_bg_dir and outlier_dir, usually "radial_bg" and
"outlier" under data_dir

makes files cbflist_file, radiallist_file, orienlist_file, outlierlist_file,
peaklist_file, usually "cbflist.dat", "radiallist.dat", "orienlist.dat",
"outlierlist.dat", "peaklist.dat", located in data_dir

After making files, calls mpi_make_background to generate data_info_file,
outlier and peak files in outlier_dir, ave_bg files in radial_bg_dir

"""

from subprocess import *
import sys
import os
import time
from pathlib import Path
import curses
import configparser

    

if (len(sys.argv) < 2) :
    print("Usage: make_background.py [config file]")
    exit()

config_file = Path(sys.argv[1])

config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
config.optionxform = str
config.read(config_file)

home_dir = Path(config['files']['home_dir'])
raw_data_dir = Path(config['files']['raw_data_dir'])
required_string = config['files']['required_string']
data_dir = Path(config['files']['data_dir'])
radial_dir = Path(config['files']['radial_bg_dir'])
outlier_dir = Path(config['files']['outlier_dir'])
cbflist_file = Path(config['files']['cbflist_file'])
radiallist_file = Path(config['files']['radiallist_file'])
orienlist_file = Path(config['files']['orienlist_file'])
peaklist_file = Path(config['files']['peaklist_file'])
outlierlist_file = Path(config['files']['outlierlist_file'])
nproc = int(config['reduce_data']['nproc'])
file_ext = config['make_detector']['file_type']


dir_set = set([])
file_list = []
for file_name in raw_data_dir.rglob(f"*.{file_ext}"):
    file_list.append(file_name)
    dir_set.add(str(file_name.parent))
dir_list = sorted(list(dir_set))
print(dir_list)
file_list = sorted(file_list)

cbf_files = []
radial_files = []
orien_files = []
outlier_files = []
peak_files = []

ctdirs = len(dir_list)
print("number of dirs ",ctdirs)
if (ctdirs < 1):
    sys.exit("No directories specified. Exiting")
ext_len=4
fileid = 0

radial_dir.mkdir(exist_ok=True)
outlier_dir.mkdir(exist_ok=True)

file_id = 1
for f in file_list:
    cbf_files.append(str(f))
    radial_files.append(str(Path(radial_dir,f"ave_bg_{file_id:06d}.bin")))
    orien_files.append(str(Path(radial_dir,f"prob_orien_{file_id:06d}.dat")))
    outlier_files.append(str(Path(outlier_dir,f"outlier_{file_id:06d}.dat")))    
    peak_files.append(str(Path(outlier_dir,f"peak_{file_id:06d}.dat")))
    file_id += 1   

with open(cbflist_file, "w") as fp:
    for line in cbf_files:
        fp.write(f"{line}\n")


with open(radiallist_file, "w") as fp:
    for line in radial_files:
        fp.write(f"{line}\n")


with open(orienlist_file, "w") as fp:
    for line in orien_files:
        fp.write(f"{line}\n")


with open(outlierlist_file, "w") as fp:
    for line in outlier_files:
        fp.write(f"{line}\n")


with open(peaklist_file, "w") as fp:
    for line in peak_files:
        fp.write(f"{line}\n")


num_raw_data = len(cbf_files)
fp = open(config_file, "r")
lines = fp.readlines()
fp.close()

for t in range(len(lines)):
    if ("num_raw_data" in lines[t]):
        lines[t] = "num_raw_data = {0:d}\n".format(num_raw_data)

fp = open(config_file, "w")
for t in range(len(lines)):
    fp.write(lines[t])
fp.close()

sys.exit("Finished")
