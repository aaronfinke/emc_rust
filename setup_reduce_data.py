#!/bin/python
from subprocess import *
import numpy as np
import sys
import os
import time

config_file = "{0:s}".format(sys.argv[1])
fp = open(config_file, "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    entry_name = lines[t].split(" = ")[0]
    if ("res_cutoff" == entry_name):
        res_cutoff = np.double(lines[t].split(" = ")[1])
    if ("high_res_cutoff" == entry_name):
        high_res_cutoff = np.double(lines[t].split(" = ")[1])
    if ("quat_file" == entry_name):
        quat_file = lines[t].split(" = ")[1].strip()
        num_coarse_div = np.int(lines[t].split("quaternion")[1].split(".")[0])
    if ("unit_cell_file" == entry_name):
        unit_cell_file = lines[t].split(" = ")[1].strip()
    if ("home_dir" == entry_name):
        home_dir = lines[t].split(" = ")[1].strip()
    if ("data_dir" == entry_name):
        data_dir = lines[t].split(" = ")[1].strip()
    if ("prob_dir" == entry_name):
        prob_dir = lines[t].split(" = ")[1].strip()
    if ("start_prob_dir" == entry_name):
        base_start_prob_dir = lines[t].split(" = ")[1].strip()
        if ("_A" in base_start_prob_dir):
          base_start_prob_dir = base_start_prob_dir.rstrip("_A")
        elif ("_B" in base_start_prob_dir):
          base_start_prob_dir = base_start_prob_dir.rstrip("_B")
    if ("reduced_cbflist_file" == entry_name):
        reduced_cbflist_file = lines[t].split(" = ")[1].strip()
    if ("nproc" == entry_name):
        nproc = np.int(lines[t].split(" = ")[1].strip())
    if ("num_raw_data" == entry_name):
        num_raw_data = np.int(lines[t].split(" = ")[1].strip())
    if ("qlen" == entry_name):
        qlen = np.int(lines[t].split(" = ")[1].strip())

fp = open (quat_file)
num_rot = np.fromfile(fp, dtype=np.int32)[0]
fp.close()

cmd = "ls {0:s}/radial_bg/prob_orien*.dat | grep -v \' 0 \' > TEMP_LS_FILE".format(data_dir)
p = Popen(cmd, shell=True)
p.wait()
cmd = "wc TEMP_LS_FILE > TEMP_WC_FILE"
p = Popen(cmd, shell=True)
p.wait()
fp = open ("TEMP_WC_FILE","r")
line = fp.readline()
num_data_used = np.int(line.strip()[0])
fp.close()
cmd = "rm TEMP_WC_FILE TEMP_LS_FILE"
p = Popen(cmd, shell=True)
p.wait()

mem_req = np.int64((num_data_used * 40016) + (num_raw_data * 1036) + \
  (nproc * (28 + (4 * num_rot))) + (qlen * 8))
mem_reqM = np.int64((mem_req/1000000) + 1)
mem_reqG = (mem_reqM/1000) + 1
print (f"estimated memory use (reduce_data) {mem_req:d} ({mem_reqG:d}G)")

cmd = "reduce_data {0:s}".format(config_file)
print ("Running:",cmd)
p = Popen(cmd, shell=True)
p.wait()
    
cmd = "write_data {0:s}".format(config_file)
print ("Running:",cmd)
p = Popen(cmd, shell=True)
p.wait()
 
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

rvol = 1.0 / vol
rsphere = (4.0/3.0) * np.pi / (res_cutoff * res_cutoff * res_cutoff)
est_num_hkl = np.ceil(rsphere / rvol)
est_max_num_hkl = np.ceil (rsphere * 1.91 / rvol)
mem_req = np.int64(((nproc * num_rot * 4) / 10000) + (nproc * num_rot * 60) + \
  (nproc * (400000 + (est_num_hkl * 84) + (est_max_num_hkl * 8))))
mem_reqM = np.int64((mem_req/1000000) + 1)
mem_reqG = (mem_reqM/1000) + 1
print (f"estimated memory use (make_Emat) {mem_req:d} ({mem_reqG:d}G)")
  
cur_dir = os.getcwd()
if (config_file.find("/") != 0) :
    config_file = os.path.join (cur_dir, config_file)

out_log = "{0:s}/reduce_data.out".format(home_dir)
if os.path.exists(out_log):
    cmd = "rm -f {0:s}".format(out_log)
    print ("Running ",cmd)
    p = Popen(cmd, shell=True)
    p.wait()

script_file = "{0:s}/reduce_data_script.sh".format(home_dir)
fp = open(script_file, "w")
tmp = "#! /bin/bash\n"
fp.write(tmp)
tmp = "kinit -k -t /nfs/chess/user/macchess/dms35/etc/dms35-keytab dms35@CLASSE.CORNELL.EDU\n"
fp.write(tmp)
tmp = "LOCAL_DIR=${TMPDIR}\n"
fp.write(tmp)
tmp = "SOURCE_DIR={0:s}\n".format(data_dir)
fp.write(tmp)
tmp = "mkdir ${LOCAL_DIR}/Data\n"
fp.write(tmp)
tmp = "rsync -az ${SOURCE_DIR}/ ${LOCAL_DIR}/Data/\n"
fp.write(tmp)
tmp = "convert_config {0:s} {1:s}\n".format(config_file,"${LOCAL_DIR}/Data")
fp.write(tmp)
local_config = "{0:s}".format("${LOCAL_DIR}/Data/config.ini")
tmp = "/usr/lib64/openmpi/bin/mpirun -np {0:d} mpi_make_Emat -low {1:s}\n".format(nproc,local_config)
fp.write(tmp)
tmp = "rsync -az ${LOCAL_DIR}/Data/ ${SOURCE_DIR}\n"
fp.write(tmp)
fp.close()

m_core = nproc
cmd = "qsub -q all.q -S /bin/bash -o {0:s} -pe sge_pe {1:d} -l m_core={2:d} -l mem_free={3:d}G -m e -M dms35@cornell.edu {4:s}".format(out_log,\
  nproc,m_core,mem_reqG,script_file)
print( "Running ",cmd)
p = Popen(cmd, shell=True)
p.wait()

#Wait for job to finish, i.e. for "qstat" to return ""
time.sleep(10)
test_file = "{0:s}/qstat_output".format(home_dir)
while True:
    cmd = "qstat | grep reduce > {0:s}".format(test_file)
    p = Popen(cmd, shell=True)
    p.wait()
    statinfo = os.stat(test_file)
    if (statinfo.st_size < 10):
        break

