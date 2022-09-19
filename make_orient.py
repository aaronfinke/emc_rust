#!/bin/python
from subprocess import *
import sys
import numpy as np
import os
import time
from pathlib import Path
import math
import configparser

if (len(sys.argv) < 2) :
    print("Usage: make_orient.py [config file]")
    exit()

config_file = Path(sys.argv[1])

config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
config.optionxform = str
config.read(config_file)

home_dir = Path(config['files']['home_dir'])
data_dir = Path(config['files']['data_dir'])
unit_cell_file = Path(config['files']['unit_cell_file'])
VN = float(config['orient_peak']['VN'])
gw = float(config['orient_peak']['gw'])
num_raw_data = int(config['make_background']['num_raw_data'])
num_rows = int(config['make_detector']['num_row'])
num_cols =int(config['make_detector']['num_col'])
res_cutoff = float(config['orient_peak']['res_cutoff'])
nproc = int(config['reduce_data']['nproc'])


with open(unit_cell_file,'r') as fp:
    lines = fp.readlines()
# Add sym ops to unit cell file if they are not already there
if (len(lines) < 3) :
    
    cmd = ["../bin/make_sym_op", "../etc/symops.dat", f"{config_file}"]
    print("Running: ",cmd)
    run(cmd,capture_output=True,text=True)

words = lines[0].split()
a = float(words[0])
b = float(words[1])
c = float(words[2])
alpha = float(words[3])
beta = float(words[4])
gamma = float(words[5])
conv = math.atan(1.0)*4.0/180.0
alpharad = alpha*conv
betarad = beta*conv
gammarad = gamma*conv
sina = math.sin(alpharad)
cosa = math.cos(alpharad)
sinb = math.sin(betarad)
cosb = math.cos(betarad)
sing = math.sin(gammarad)
cosg = math.cos(gammarad)
cosas = (cosg*cosb-cosa)/ (sinb*sing)
sinas = math.sqrt(1.0-cosas*cosas)
cosbs = (cosa*cosg-cosb)/ (sina*sing)
sinbs = math.sqrt(1.0-cosbs*cosbs)
cosgs = (cosa*cosb-cosg)/ (sina*sinb)
sings = math.sqrt(1.0-cosgs*cosgs)

sumc = (alpharad+betarad+gammarad)*0.5
vol = 2.0*a*b*c*math.sqrt(math.sin(sumc-alpharad)*math.sin(sumc-betarad)*math.sin(sumc-gammarad)*math.sin(sumc));
recipa = b * c * sina/vol
recipb = c * a * sinb/vol
recipc = a * b * sing/vol
min_rcell = min([recipa, recipb, recipc])

rvol = 1.0 / vol
rsphere = (4.0/3.0) * math.pi / (res_cutoff * res_cutoff * res_cutoff)
est_num_hkl = math.ceil(rsphere / rvol)

delta_theta = 2.0*gw*res_cutoff*min_rcell/VN
num_div = int(0.944/delta_theta)
print("delta_theta, num_div", delta_theta, num_div)

cmd = ["../bin/make_quaternion", config_file, "-bin", str(num_div)]
print("Running... ",cmd)
p = run(cmd, capture_output=True)
print(p.stdout)

filepath = Path(data_dir / f"c_quaternion{num_div:0d}.bin")

with open(filepath,"r") as fp:
    num_rot = np.fromfile(fp, dtype=int)[0]

config['files']['quat_file'] = str(filepath)

with open(config_file,"w") as fp:
    config.write(fp)

total_pix = num_rows * num_cols
print("est_num_hkl ",est_num_hkl)
print("num_raw_data ",num_raw_data)
print("num_rot ",num_rot)
print("total_pix ",total_pix)
mem_req = int((est_num_hkl * 12 * nproc) + \
  (num_raw_data * (536 + 480 + 40000)) + (num_rot * 112 * nproc) + \
  (total_pix * 24 * nproc)) 
# Use numbers below for test version that stores extra information
#mem_req = int((est_num_hkl * 16 * nproc) + \
#  (num_raw_data * (536 + 480 + 880000 + 20)) + (num_rot * 112 * nproc) + \
#  (total_pix * 24 * nproc))
mem_reqM = int((mem_req/1000000) + 1)
mem_reqG = (mem_reqM/1000) + 1
print(f"memory {mem_req} ({mem_reqG}G)")
if (mem_reqG > 256) :
    print("More than 256 GB of memory requested - job aborted.")
    print("To change this test, edit make_orient.py.")
    exit()

sys.exit()
out_log = Path(home_dir / 'make_orient.out')
if out_log.is_file():
    out_log.unlink()
    

script_file = Path(Path.cwd() / 'make_orient_script.sh')
script_header = f"""#!/bin/sh -eu
#SBATCH -t 3-0:0
#SBATCH -n {nproc} 
#SBATCH -N 1
#SBATCH -J emc_mk_orient
#SBATCH --mem=GB
#SBATCH -o {out_log}
module purge
module load GCC/11.2.0 OpenMPI/4.1.1 HDF5/1.12.1
"""
with open(script_file,"w") as fp:
    fp.write(script_header)
    fp.write(f"mpirun -np {nproc} ../bin/mpi_sync_orient_peak {config_file}")

cmd = ['sbatch',script_file]
proc = run(cmd, capture_output=True, text=True)
print(proc.stdout)
jobid = proc.stdout.split()[-1]

time.sleep(5)

cmd = ['squeue','-j',jobid]
progress = run(cmd, capture_output=True, text=True)
while jobid in progress.stdout:
    time.sleep(1)
    progress = run(cmd, capture_output=True, text=True)
sys.exit("Finished")





