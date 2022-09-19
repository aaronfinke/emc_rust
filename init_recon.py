#!/bin/python
"""
Usage: python init_recon.py [config file]

Define project directory by "home_dir" in config file.
Creates subdirectory "Data" under project directory, and "high_prob" under 
this. Other filenames can be specified in config file, but if not will default
to standard locations relative to project directory.

"""
from subprocess import *
import sys
import os
import configparser
from pathlib import Path

def main():
    config_file = f"{sys.argv[1]}"
    try:
        cleanup = sys.argv[2]
    except:
        cleanup = 'False'

    config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
    config.optionxform = str
    config.read(config_file)

    if cleanup.lower() == 'clean':
        config = cleanUpConfig(config)
    else: 
        config = writeConfig(config)

    with open(config_file,'w') as fp:
        config.write(fp)

def cleanUpConfig(config):
        
    config['files']['required_string'] =  ''
    config['files']['data_dir'] =  ''
    config['files']['prob_dir'] =  ''
    config['files']['reduced_data_dir'] =  ''
    config['files']['radial_bg_dir'] =  ''
    config['files']['outlier_dir'] =  ''
    config['files']['start_prob_dir'] =  ''
    config['files']['mask_file'] =  ''
    config['files']['pixmap_file'] =  ''
    config['files']['rec2pix_file'] =  ''
    config['files']['rec_vectors_file'] =  ''
    config['files']['unit_cell_file'] =  ''
    config['files']['quat_file'] =  ''
    config['files']['start_phi_file'] =  ''
    config['files']['start_intens_file'] =  ''
    config['files']['cbflist_file'] =  ''
    config['files']['radiallist_file'] =  ''
    config['files']['orienlist_file'] =  ''
    config['files']['outlierlist_file'] =  ''
    config['files']['peaklist_file'] =  ''
    config['files']['data_info_file'] =  ''
    config['files']['reduced_data_id_file'] =  ''
    config['files']['reduced_peaklist_file'] =  ''
    config['files']['reduced_cbflist_file'] =  ''
    config['files']['reduced_radiallist_file'] =  ''
    config['files']['prob_orien_file'] =  ''
    config['files']['mpi_bgfile'] =  ''
    config['files']['mpi_datafile'] =  ''
    config['files']['r2peak_file'] =  ''
    config['files']['peak2r_file'] =  ''
    config['files']['fine_quat_file'] =  ''
    config['files']['quat_table_file'] =  ''
    config['files']['local_r2peak_file'] =  ''
    config['files']['local_peak2r_file'] =  ''
    return config

def writeConfig(config):
    noneEntries = ['default','Default','none','None']

    home_dir = config['files']['home_dir']
    if not home_dir or any(y in str(home_dir) for y in noneEntries):
        sys.exit("You must specify home_dir in config.ini")
    home_dir = Path(config['files']['home_dir'])
    raw_data_dir = config['files']['raw_data_dir']
    if not raw_data_dir or any(y in str(raw_data_dir) for y in noneEntries):
        sys.exit("You must specify raw_data_dir in config.ini")
    required_string = config['files']['required_string']
    if not required_string or any(y in str(required_string) for y in noneEntries):
        config['files']['required_string'] = 'None'
    data_dir = config['files']['data_dir']
    if not data_dir or any(y in str(data_dir) for y in noneEntries):
        data_dir = Path(home_dir / 'Data')
        config['files']['data_dir'] = str(data_dir)
    if not Path(config['files']['data_dir']).exists():
        Path(config['files']['data_dir']).mkdir()
    data_dir = Path(config['files']['data_dir'])
    prob_dir = config['files']['prob_dir']
    if not prob_dir or any(y in str(prob_dir) for y in noneEntries):
        config['files']['prob_dir'] = str(Path(data_dir / 'high_prob'))
    if not Path(config['files']['prob_dir']).exists():
        Path(config['files']['prob_dir']).mkdir()
    prob_dir = Path(config['files']['prob_dir'])
    reduced_data_dir = config['files']['reduced_data_dir']
    if not reduced_data_dir or any(y in str(reduced_data_dir) for y in noneEntries):
        config['files']['reduced_data_dir'] = str(Path(home_dir / 'reduced_data'))
        if not Path(home_dir / 'reduced_data').exists():
            Path(home_dir / 'reduced_data').mkdir(exist_ok=True)
    radial_bg_dir = config['files']['radial_bg_dir']
    if not radial_bg_dir or any(y in str(radial_bg_dir) for y in noneEntries):
        config['files']['radial_bg_dir'] = str(Path(data_dir / 'radial_bg'))
    outlier_dir = config['files']['outlier_dir']
    if not outlier_dir or any(y in str(outlier_dir) for y in noneEntries):
        config['files']['outlier_dir'] = str(Path(data_dir / 'outlier'))
    start_prob_dir = config['files']['start_prob_dir']
    if not start_prob_dir or any(y in str(start_prob_dir) for y in noneEntries):
        config['files']['start_prob_dir'] = str(Path(data_dir / 'start_prob'))
    mask_file = config['files']['mask_file']
    if not mask_file or any(y in str(mask_file) for y in noneEntries):
        config['files']['mask_file'] = str(Path(home_dir / 'mask.h5'))
    pixmap_file = config['files']['pixmap_file']
    if not pixmap_file or any(y in str(pixmap_file) for y in noneEntries):
        config['files']['pixmap_file'] = str(Path(home_dir / 'pix_map.h5'))
    rec2pix_file = config['files']['rec2pix_file']
    if not rec2pix_file or any(y in str(rec2pix_file) for y in noneEntries):
        config['files']['rec2pix_file'] = str(Path(home_dir / 'rec2pix_map.h5'))
    rec_vectors_file = config['files']['rec_vectors_file']
    if not rec_vectors_file or any(y in str(rec_vectors_file) for y in noneEntries):
        config['files']['rec_vectors_file'] = str(Path(home_dir / 'rec_vectors.h5'))
    unit_cell_file = config['files']['unit_cell_file']
    if not unit_cell_file or any(y in str(unit_cell_file) for y in noneEntries):
        config['files']['unit_cell_file'] = str(Path(home_dir / 'unit_cell.dat'))
    quat_file = config['files']['quat_file']
    if not quat_file or any(y in str(quat_file) for y in noneEntries):
        config['files']['quat_file'] = str(Path(data_dir / 'c_quaternion70.bin'))
    start_phi_file = config['files']['start_phi_file']
    if not start_phi_file or any(y in str(start_phi_file) for y in noneEntries):
        config['files']['start_phi_file'] = str(Path(data_dir / 'start_phi.dat'))
    start_intens_file = config['files']['start_intens_file']
    if not start_intens_file or any(y in str(start_intens_file) for y in noneEntries):
        config['files']['start_intens_file'] = str(Path(data_dir / 'start_intensity.bin'))
    cbflist_file = config['files']['cbflist_file']
    if not cbflist_file or any(y in str(cbflist_file) for y in noneEntries):
        config['files']['cbflist_file'] = str(Path(data_dir / 'cbflist.dat'))
    radiallist_file = config['files']['radiallist_file']
    if not radiallist_file or any(y in str(radiallist_file) for y in noneEntries):
        config['files']['radiallist_file'] = str(Path(data_dir / 'radiallist.dat'))
    orienlist_file = config['files']['orienlist_file']
    if not orienlist_file or any(y in str(orienlist_file) for y in noneEntries):
        config['files']['orienlist_file'] = str(Path(data_dir / 'orienlist.dat'))
    outlierlist_file = config['files']['outlierlist_file']
    if not outlierlist_file or any(y in str(outlierlist_file) for y in noneEntries):
        config['files']['outlierlist_file'] = str(Path(data_dir / 'outlierlist.dat'))
    peaklist_file = config['files']['peaklist_file']
    if not peaklist_file or any(y in str(peaklist_file) for y in noneEntries):
        config['files']['peaklist_file'] = str(Path(data_dir / 'peaklist.dat'))
    data_info_file = config['files']['data_info_file']
    if not data_info_file or any(y in str(data_info_file) for y in noneEntries):
        config['files']['data_info_file'] = str(Path(data_dir / 'data_info.dat'))
    reduced_data_id_file = config['files']['reduced_data_id_file']
    if not reduced_data_id_file or any(y in str(reduced_data_id_file) for y in noneEntries):
        config['files']['reduced_data_id_file'] = str(Path(data_dir / 'reduced_data_id.dat'))
    reduced_peaklist_file = config['files']['reduced_peaklist_file']
    if not reduced_peaklist_file or any(y in str(reduced_peaklist_file) for y in noneEntries):
        config['files']['reduced_peaklist_file'] = str(Path(data_dir / 'reduced_peaklist.dat'))
    reduced_cbflist_file = config['files']['reduced_cbflist_file']
    if not reduced_cbflist_file or any(y in str(reduced_cbflist_file) for y in noneEntries):
        config['files']['reduced_cbflist_file'] = str(Path(data_dir / 'reduced_cbflist.dat'))
    reduced_radiallist_file = config['files']['reduced_radiallist_file']
    if not reduced_radiallist_file or any(y in str(reduced_radiallist_file) for y in noneEntries):
        config['files']['reduced_radiallist_file'] = str(Path(data_dir / 'reduced_radiallist.dat'))
    prob_orien_file = config['files']['prob_orien_file']
    if not prob_orien_file or any(y in str(prob_orien_file) for y in noneEntries):
        config['files']['prob_orien_file'] = str(Path(data_dir / 'prob_orien.bin'))
    mpi_bgfile = config['files']['mpi_bgfile']
    if not mpi_bgfile or any(y in str(mpi_bgfile) for y in noneEntries):
        config['files']['mpi_bgfile'] = str(Path(data_dir / 'mpi_bg_model.bin'))
    mpi_datafile = config['files']['mpi_datafile']
    if not mpi_bgfile or any(y in str(mpi_bgfile) for y in noneEntries):
        config['files']['mpi_datafile'] = str(Path(data_dir / 'mpi_datafile.bin'))
    r2peak_file = config['files']['r2peak_file']
    if not r2peak_file or any(y in str(r2peak_file) for y in noneEntries):
        config['files']['r2peak_file'] = str(Path(data_dir / 'mpi_r2peak_low.bin'))
    peak2r_file = config['files']['peak2r_file']
    if not peak2r_file or any(y in str(peak2r_file) for y in noneEntries):
        config['files']['peak2r_file'] = str(Path(data_dir / 'mpi_peak2r_low.bin'))
    fine_quat_file = config['files']['fine_quat_file']
    if not fine_quat_file or any(y in str(fine_quat_file) for y in noneEntries):
        config['files']['fine_quat_file'] = str(Path(data_dir / 'c_reduced_quaternion100.bin'))
    quat_table_file = config['files']['quat_table_file']
    if not quat_table_file or any(y in str(quat_table_file) for y in noneEntries):
        config['files']['quat_table_file'] = str(Path(data_dir / 'reduced_quaternion_70_100.dat'))
    local_r2peak_file = config['files']['local_r2peak_file']
    if not local_r2peak_file or any(y in str(local_r2peak_file) for y in noneEntries):
        config['files']['local_r2peak_file'] = str(Path(data_dir / 'mpi_r2peak_local.bin'))
    local_peak2r_file = config['files']['local_peak2r_file']
    if not local_peak2r_file or any(y in str(local_peak2r_file) for y in noneEntries):
        config['files']['local_peak2r_file'] = str(Path(data_dir / 'mpi_peak2r_local.bin'))
    return config

if __name__ == '__main__':
    main()