import numpy as np
from pathlib import Path
import sys
import configparser
from math import sqrt, ceil, log, exp
from numba import njit

nbin = 3000
sampling_sz = 2.0

def main(argv):
    if len(argv) < 2:
        print(f"The Config file is missing!")
    
    config_file = Path(argv[1])
    config = configparser.ConfigParser(comment_prefixes='/', allow_no_value=True)
    config.optionxform = str
    config.read(config_file)

    num_row = int(config['make_detector']['num_row'])
    num_col = int(config['make_detector']['num_col'])
    num_data = int(config['make_background']['num_raw_data'])
    min_patch_sz = int(config['make_powder']['min_patch_sz'])
    max_patch_sz = int(config['make_powder']['max_patch_sz'])
    min_num_peak = int(config['make_powder']['min_num_peak'])
    max_num_peak = int(config['make_powder']['max_num_peak'])
    
    detd = float(config['make_detector']['detd'])
    beam_vec = np.array([float(config.get('make_detector','sx')),float(config.get('make_detector','sy')),float(config.get('make_detector','sz'))]) 
    wl = float(config['make_detector']['wl'])
    px = float(config['make_detector']['px'])
    cx = float(config['make_detector']['cx'])
    cy = float(config['make_detector']['cy'])
    home_dir = config['files']['home_dir']
    data_dir = config['files']['data_dir']
    outlier_dir = config['files']['outlier_dir']
    peaklist_file = config['files']['peaklist_file']

    norm = sqrt(beam_vec[0]*beam_vec[0] + beam_vec[1]*beam_vec[1] + beam_vec[2]*beam_vec[2]) 
    D = (detd/px) * norm/abs(beam_vec[2]) 
    wl_D = wl*D 

    beam_vec[0] *= D/norm 
    beam_vec[1] *= D/norm 
    beam_vec[2] *= D/norm 

    total_pix = num_row*num_col
    pix = np.zeros((total_pix,3), dtype = np.double)
    pix = make_pix(total_pix, D, num_col, cx, cy, beam_vec, pix)

    try:
        with open(peaklist_file,'r') as fp:
            peakfiles = fp.read().splitlines()
    except:
        sys.exit('peaklist_file not found')

    count_relevant, radial_hist, frame_peak_count, patch_sz_hist, powder2d = \
        make_powder(total_pix,num_data, nbin, peakfiles, min_patch_sz, max_patch_sz, \
                pix, min_num_peak, max_num_peak, sampling_sz)

    print(f"number of relevant frames = {count_relevant}")
    rescale = 1./wl_D/sampling_sz

    with open( Path(data_dir , '1d_pseudo_powder.dat'), "w" ) as fp:
        for i in range(nbin):
            print(f"{i*rescale:1.5e} {radial_hist[i][0]:1.5e} {radial_hist[i][1]:1.5e}", file=fp)

    with open( Path(data_dir , 'frame_peak_count.dat'), "w" ) as fp:
        for d in range(num_data):
            print(f"{frame_peak_count[d]:d}", file=fp)
    
    with open( Path(data_dir , 'patch_sz_count.dat'), "w" ) as fp:
        for i in range(min_patch_sz, nbin):
            print(f"{i:d} {patch_sz_hist[i]:d}", file= fp)

    with open( Path(data_dir , '2d_pseudo_powder.dat'), "w" ) as fp:
        for i in range(total_pix):
            print(f"{powder2d[i]:d}", file=fp)
        print(f"{count_relevant:d}", file=fp)

@njit
def make_pix(total_pix, D, num_col, cx, cy, beam_vec, pix):
    for t in range(total_pix):
        x = int(t / num_col) - cx
        y = t % num_col - cy
        sx = beam_vec[0] + x
        sy = beam_vec[1] + y
        sz = beam_vec[2]
        norm = sqrt(sx*sx + sy*sy + sz*sz)
        sx *= D/norm
        sy *= D/norm
        sz *= D/norm
        pix[t][0] = sx - beam_vec[0] 
        pix[t][1] = sy - beam_vec[1] 
        pix[t][2] = sz - beam_vec[2] 
    return pix


def make_powder(total_pix,num_data, nbin, peakfiles, min_patch_sz, max_patch_sz, \
                pix, min_num_peak, max_num_peak, sampling_sz):
    peak_rvec = np.zeros((total_pix,3), dtype = float)
    reflection = np.empty(total_pix, dtype=float)
    frame_peak_count = np.empty(num_data, dtype=int)
    patch_sz_hist = np.empty(nbin, dtype=int)
    powder2d = np.empty(total_pix, dtype = int)
    radial_hist = np.zeros((nbin,2), dtype = float)

    count_relevant = 0
    for d, peakfile in enumerate(peakfiles):
        with open(peakfile,"r") as fp:
            num_peak = 0
            for line in fp:
                num_pix = int(line)
                patch_sz_hist[num_pix] += 1
                peakline = fp.readline()
                peakline = [int(i) for i in peakline.split()]
                it = iter(peakline)
                peaks = list(zip(it,it))

                if (num_pix < min_patch_sz or num_pix > max_patch_sz):
                    continue
                else:
                    val_sum = 0.0
                    reflection[num_peak] = 0.0
                    weighted_rvec = np.zeros(3)
                    for peak in peaks:
                        pix_id = peak[0]
                        pix_val = peak[1]
                        if powder2d[pix_id] < pix_val:
                            powder2d[pix_id] = pix_val
                        reflection[num_peak] += pix_val
                        weighted_rvec += pix[pix_id] * pix_val
                        val_sum += pix_val
                    peak_rvec[num_peak] = weighted_rvec/val_sum
                    num_peak += 1

        if num_peak < min_num_peak or num_peak > max_num_peak:
            continue
        frame_peak_count[d] = num_peak
        count_relevant += 1
        for r in range(num_peak):
            for t in range(r):
                dist = 0.0
                diff_rvec = peak_rvec[r] - peak_rvec[t]
                dist = np.sqrt(np.sum(diff_rvec ** 2))*sampling_sz
                idx = int(round(dist))
                if idx < nbin:
                    radial_hist[idx][0] += 1
        for r in range(num_peak):
            dist = 0.0
            for t in range(r):
                diff_rvec = peak_rvec[r] - peak_rvec[t]
                dist = np.sqrt(np.sum(diff_rvec ** 2))*sampling_sz
                idx = int(round(dist))
                if idx < nbin:
                    radial_hist[idx][0] += 1
    
    return count_relevant, radial_hist, frame_peak_count, patch_sz_hist, powder2d


                        









if __name__ == "__main__":
    main(sys.argv)