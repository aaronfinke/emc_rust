use cryiorust::cbf::Cbf;

use cryiorust::frame::Frame;

use make_background::{Config, FrameInfo, PeakLocation};
use rayon::prelude::*;

use std::env;
use std::fs::File;
use std::io::{prelude::*, BufWriter};

use std::error::Error;
use std::process;
use std::time::Instant;
use std::sync::{Arc, Mutex};

#[allow(unused)]
fn main() {
    let now = Instant::now();
    let vmax: f64 = 200.0;
    let cdf_thres = 0.99999;
    let dv = 0.01;
    let num_iter = 5;
    let max_num_peak_on_ring = 5;
    let bg_thres = 1.05;

    let args: Vec<String> = env::args().collect();
    match args.len() {
        2 => {}
        _ => panic!("execution: ./make_background (/path/to/)config.ini"),
    };

    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1);
    });


    let len_val: usize = ((vmax / dv) as f64).round() as usize;
    let mut count_thres: Vec<i32> = {
        let mut val: Vec<f64> = vec![Default::default(); len_val];
        let mut log_val: Vec<f64> = vec![Default::default(); len_val];
        let mut cdf: Vec<f64> = vec![Default::default(); len_val];
        let mut count_thres: Vec<i32> = vec![0; len_val];

        let kmax = (vmax + 10.0 * (vmax).sqrt()) as i32;
        let mut log_fact: Vec<f64> = vec![Default::default(); (kmax + 1).try_into().unwrap()];
        log_fact[0] = 0.0;
        for i in 1..(kmax + 1) as usize {
            log_fact[i] = log_fact[i - 1] + (i as f64).ln();
        }
        for i in 0..len_val {
            val[i] = (i + 1) as f64 * dv;
            log_val[i] = val[i].ln();
            cdf[i] = (-val[i]).exp();
        }
        for i in 0..len_val {
            while cdf[i] < cdf_thres {
                count_thres[i] += 1;
                let logP = -val[i] + count_thres[i] as f64 * log_val[i]
                    - log_fact[count_thres[i] as usize];
                cdf[i] += logP.exp();
            }
        }
        count_thres
    };
    //
    // make_bg section
    //

    let mut data_info = Arc::new(Mutex::new(FrameInfo {
        row_size: 0,
        col_size: 0,
        exposure: 0.0,
        wl: 0.0,
        detd: 0.0,
        c_row: 0.0,
        c_col: 0.0,
    }));

    let dmin = 0;
    let dmax = 1000;
    let bits_per = 32;
    (dmin..dmax).into_par_iter().for_each(|d| {
        let mut peak_id: Vec<i32> = vec![0; config.total_pix];
        let mut peak_val: Vec<i32> = vec![0; config.total_pix];
        let mut peak_label: Vec<i32> = vec![0; config.total_pix];
        let mut peak_map: Vec<i32> = vec![-1; config.total_pix];
        let mut queue: Vec<i32> = vec![0; config.total_pix];
        let mut peak_list: Vec<i32> = vec![0; config.total_pix];
    
        let mut radial_ct: Vec<i32> = vec![0; config.qlen as usize];
        let mut peak_location: Vec<PeakLocation> = vec![PeakLocation(0, 0); config.qlen as usize];
        let mut radial_val: Vec<f64> = vec![0.0; config.qlen as usize];
        let mut radial_weight: Vec<f64> = vec![0.0; config.qlen as usize];
        let mut ave_bg: Vec<f64> = vec![Default::default(); config.qlen as usize];
        let mut ave_bge: Vec<f64> = vec![Default::default(); config.qlen as usize];
        let mut data_frame: Vec<i32> = vec![-1; config.num_pix as usize];
        let det = load_cbf(&config.cbf_files[d], &mut data_info.lock().unwrap()).unwrap_or_else(|err| {
            println!("Problem loading cbf: {err}");
            process::exit(1);
        });
        for (t, pixel) in data_frame.iter_mut().enumerate() {
            let pid: usize = config.rec2pix_map[t] as usize;
            if det[pid] >= 0.0 && det[pid] <= config.hot_pix_thres.into() {
                *pixel = det[pid] as i32
            }
        }
        radial_val.fill(0.0);
        radial_weight.fill(0.0);
        radial_ct.fill(0);

        for t in 0..config.num_pix as usize {
            let qid = config.qid_map[t];
            if qid < 0 {
                continue;
            }
            if data_frame[t] < 0 {
                continue;
            }
            radial_val[qid as usize] += data_frame[t] as f64;
            radial_weight[qid as usize] += config.pix[t].3;
            radial_ct[qid as usize] += 1;
        }
        for t in 0..config.qlen as usize {
            if radial_ct[t] > 0 {
                ave_bg[t] = radial_val[t] / radial_weight[t];
            } else {
                ave_bg[t] = -1.0;
            }
            ave_bge[t] = ave_bg[t];
        }
        for i in 0..num_iter {
            radial_val.fill(0.0);
            radial_weight.fill(0.0);
            radial_ct.fill(0);
            for t in 0..config.num_pix as usize {
                let qid = config.qid_map[t];
                if qid < 0 {
                    continue;
                }
                let photon_count = data_frame[t];
                if photon_count < 0 || ave_bg[qid as usize] < 0. {
                    continue;
                };
                let mut avg_val = ave_bg[qid as usize] * config.pix[t].3;
                if avg_val < dv {
                    avg_val = dv
                };
                let mut idx: i32 = ((avg_val * (1.0 / dv) - 1.).round()) as i32;
                let thres: i32 = {
                    if idx >= len_val as i32 {
                        ad_hoc_thres(avg_val)
                    } else {
                        count_thres[idx as usize]
                    }
                };
                if photon_count > thres {
                    continue;
                }
                radial_val[qid as usize] += photon_count as f64;
                radial_weight[qid as usize] += config.pix[t].3;
                radial_ct[qid as usize] += 1;
            }
            for t in 0..config.qlen as usize {
                if radial_ct[t] > 0 {
                    ave_bg[t] = radial_val[t] / radial_weight[t]
                } else {
                    ave_bg[t] = -1.0
                }
                ave_bge[t] = ave_bg[t];
            }
        }
        let mut outlier_ct = 0;
        for t in 0..config.num_pix as usize {
            let qid = config.qid_map[t];
            if qid < 0 {
                continue;
            }
            let photon_count = data_frame[t];
            if photon_count < 0 || ave_bg[qid as usize] < 0. {
                continue;
            };
            let mut avg_val = ave_bg[qid as usize] * config.pix[t].3;
            if avg_val < dv {
                avg_val = dv
            };
            let mut idx: i32 = ((avg_val * (1.0 / dv) - 1.).round()) as i32;
            let thres: i32 = {
                if idx >= len_val as i32 {
                    ad_hoc_thres(avg_val)
                } else {
                    count_thres[idx as usize]
                }
            };
            if photon_count > thres {
                peak_id[outlier_ct] = config.rec2pix_map[t];
                peak_val[outlier_ct] = photon_count;
                peak_map[peak_id[outlier_ct] as usize] = outlier_ct as i32;
                peak_label[outlier_ct] = -1;
                outlier_ct += 1;
            }
        }
        mask_rings(
            &mut ave_bg,
            &mut ave_bge,
            &mut peak_location,
            config.qlen,
            config.qmin,
            config.dq,
            config.n_qexc as i32,
            &config.qexc,
        );

        //peak_finder section
        {
            let mut nn_pid = vec![0; 4];
            let num_nn = 4;
            let mut weighted_rvec: Vec<f64> = vec![0.0; 3];

            let inv_dq = 1.0 / config.dq;
            radial_ct.fill(0);

            // filter peaks based on number of connected pixels
            let mut num_patch = 0;
            let mut num_peak = 0;
            let mut ct = 0;
            for i in 0..outlier_ct {
                if peak_label[i] > -1 {
                    continue;
                };
                peak_label[i] = num_patch;
                let mut patch_sz = 0;
                queue[patch_sz] = i as i32;
                patch_sz += 1;
                let mut cur_idx = 0;
                while cur_idx < patch_sz {
                    nn_pid[0] = peak_id[queue[cur_idx] as usize] - config.num_col;
                    nn_pid[1] = peak_id[queue[cur_idx] as usize] - 1;
                    nn_pid[2] = peak_id[queue[cur_idx] as usize] + 1;
                    nn_pid[3] = peak_id[queue[cur_idx] as usize] + config.num_col;
                    for t in 0..num_nn {
                        if nn_pid[t] < 0 || nn_pid[t] > (config.total_pix as i32) - 1 {
                            continue;
                        };
                        if peak_map[nn_pid[t] as usize] > -1 {
                            if peak_label[peak_map[nn_pid[t] as usize] as usize] > -1 {
                                continue;
                            };
                            queue[patch_sz] = peak_map[nn_pid[t] as usize];
                            peak_label[queue[patch_sz] as usize] = num_patch;
                            patch_sz += 1;
                        }
                    }
                    cur_idx += 1;
                }
                // at least 2 connected pixels
                if patch_sz > 1 {
                    peak_list[ct] = patch_sz as i32;
                    ct += 1;
                    for t in 0..patch_sz {
                        peak_list[ct + 2 * t] = peak_id[queue[t] as usize];
                        peak_list[ct + 2 * t + 1] = peak_val[queue[t] as usize];
                    }
                    ct += 2 * patch_sz;
                    num_peak += 1;
                }
                num_patch += 1;
            }
            ct = 0;
            for i in 0..num_peak {
                let patch_sz: usize = peak_list[ct] as usize;
                ct += 1;
                let mut val_sum = 0.0;
                weighted_rvec = vec![0.0, 0.0, 0.0];
                for t in 0..patch_sz {
                    let pid: i32 = config.pix_map[peak_list[ct + 2 * t as usize] as usize];
                    let val: f64 = peak_list[ct + (2 * t) as usize + 1] as f64;
                    weighted_rvec[0] += config.pix[pid as usize].0 * val;
                    weighted_rvec[1] += config.pix[pid as usize].1 * val;
                    weighted_rvec[2] += config.pix[pid as usize].2 * val;
                    val_sum += val;
                }

                let mut norm = 0.0;
                for j in 0..3 {
                    weighted_rvec[j] /= val_sum;
                    norm += weighted_rvec[j] * weighted_rvec[j];
                }
                let qid: i32 = ((norm).sqrt() * inv_dq - 0.5).round() as i32;
                if qid >= 0 && qid < config.qlen {
                    radial_ct[qid as usize] += 1
                };
                ct += 2 * patch_sz;
            }

            let mut ring_ct = 0;
            for i in 0..config.qlen {
                if radial_ct[i as usize] > max_num_peak_on_ring {
                    queue[ring_ct] = i;
                    ring_ct += 1;
                }
            }

            //write the peakfile
            {
                ct = 0;
                let mut f = match File::create(&config.peak_files[d]) {
                    Err(why) => panic!("couldn't create {} {}", &config.peak_files[d], why),
                    Ok(file) => file,
                };
                let mut fp = BufWriter::new(f);
                for i in 0..num_peak {
                    let patch_sz = peak_list[ct];
                    ct += 1;
                    let mut val_sum = 0.0;
                    weighted_rvec = vec![0.0, 0.0, 0.0];
                    for t in 0..patch_sz {
                        let pid: i32 = config.pix_map[peak_list[ct + 2 * t as usize] as usize];
                        let val: f64 = peak_list[ct + (2 * t) as usize + 1] as f64;
                        weighted_rvec[0] += config.pix[pid as usize].0 * val;
                        weighted_rvec[1] += config.pix[pid as usize].1 * val;
                        weighted_rvec[2] += config.pix[pid as usize].2 * val;
                        val_sum += val;
                    }
                    let mut norm = 0.0;
                    for j in 0..3 {
                        weighted_rvec[j] /= val_sum;
                        norm += weighted_rvec[j] * weighted_rvec[j];
                    }
                    let qid: i32 = ((norm).sqrt() * inv_dq - 0.5).round() as i32;

                    let mut if_ring = false;
                    for t in 0..ring_ct {
                        if qid == queue[t] {
                            if_ring = true
                        };
                    }
                    if qid > -1 && qid < config.qlen && ave_bge[qid as usize] < 0.0 {
                        if_ring = true
                    }
                    if !if_ring {
                        writeln!(&mut fp, "{}", patch_sz);

                        for t in 0..patch_sz {
                            write!(
                                &mut fp,
                                "{} {} ",
                                peak_list[ct + (2 * t as usize)],
                                peak_list[ct + (2 * t as usize) + 1]
                            );
                        }
                        writeln!(&mut fp, "");
                    } else {
                        if ave_bg[qid as usize] > 0.0 {
                            println!("Reject peak {},{} in excluded ring with qid {}", d, i, qid);
                        }

                        let mut qid_min = qid;
                        let mut qid_max = qid;
                        for t in 0..patch_sz {
                            let pid = config.pix_map
                                [peak_list[((ct as i32) + (2 * t)) as usize] as usize];
                            let mut norm = 0.0;
                            norm += config.pix[pid as usize].0 * config.pix[pid as usize].0;
                            norm += config.pix[pid as usize].1 * config.pix[pid as usize].1;
                            norm += config.pix[pid as usize].2 * config.pix[pid as usize].2;
                            let qid: i32 = ((norm).sqrt() * inv_dq - 0.5).round() as i32;
                            if qid_min > qid {
                                qid_min = qid
                            };
                            if qid_max < qid {
                                qid_max = qid
                            };
                        }
                        qid_max += 1;
                        if qid_min < 0 {
                            qid_min = 0
                        };
                        if qid_max > config.qlen - 1 {
                            qid_max = config.qlen
                        };

                        for t in qid_min..qid_max {
                            if ave_bge[t as usize] > 0.0 {
                                ave_bg[t as usize] = -1.0;
                                ave_bge[t as usize] = -1.0;
                            }
                        }
                    }
                    ct += 2 * patch_sz as usize;
                }
            } //close peakfile
        } //end peak_finder section
          //write to radial_files
        {
            let mut fp = match File::create(&config.radial_files[d]) {
                Err(why) => panic!("couldn't create {} {}", &config.radial_files[d], why),
                Ok(file) => file,
            };
            for i in &ave_bg {
                fp.write(&i.to_le_bytes())
                    .expect("Unable to write to radial_file");
            }
        }
        //write to outlier_files
        {
            let mut fp = match File::create(&config.outlier_files[d]) {
                Err(why) => panic!("couldn't create {} {}", &config.outlier_files[d], why),
                Ok(file) => file,
            };
            for t in 0..outlier_ct {
                let pid = config.pix_map[peak_id[t] as usize];
                let qid = config.qid_map[pid as usize];
                if qid < 0 || ave_bg[qid as usize] < 0.0 {
                    continue;
                };
                writeln!(fp, "{} {}", peak_id[t], peak_val[t]);
            }
        }
    });

    let mut fp = match File::create(&config.data_info_file) {
        Err(why) => panic!("couldn't create {} {}", &config.data_info_file, why),
        Ok(file) => file,
    };
    let d1 = data_info.lock().unwrap();
    let tsrting = format!(
        "{0} {1} {2:.1} {3:.1} {4:1.5e} {5:1.5e} {6:1.5e}",
        d1.row_size,
        d1.col_size,
        d1.c_row,
        d1.c_col,
        d1.exposure,
        d1.wl,
        d1.detd
    );
    write!(fp, "{}", tsrting);

    let elapsed = now.elapsed();
    println!("Time elapsed: {:?}", elapsed);
}

fn mask_rings(
    ave_bg: &mut [f64],
    ave_bge: &mut [f64],
    peak_location: &mut [PeakLocation],
    qlen: i32,
    qmin: i32,
    dq: f64,
    n_qexc: i32,
    qexc: &[Vec<i32>],
) {
    let bg_thres = 1.05;
    let half_width = 3;
    let idx_start = (qmin as f64 / dq).ceil() as i32;
    let idx_end = qlen;
    let mut ct = 0;
    for l in idx_start..idx_end {
        let i: usize = usize::try_from(l).unwrap();
        if ave_bg[i] < 0.0 {
            continue;
        };

        let mut flag = 0;
        let mut cur_val = ave_bg[i];

        for j in 1..half_width {
            let idx = i + j;
            if idx > (idx_end as usize) - 1 || ave_bg[idx] < 0.0 {
                continue;
            };
            if ave_bg[idx] < cur_val {
                cur_val = ave_bg[idx]
            } else {
                flag = 1
            };
        }

        if flag == 1 {
            continue;
        }

        let mut jmin = l;
        loop {
            if jmin == idx_start || ave_bg[(jmin - 1) as usize] > ave_bg[jmin as usize] {
                break;
            };
            jmin -= 1;
        }

        let mut jmax = l;
        loop {
            if jmax == idx_end - 1 || ave_bg[(jmax + 1) as usize] > ave_bg[jmax as usize] {
                break;
            };
            jmax += 1;
        }

        cur_val =
            (l - jmin) as f64 * ave_bg[jmax as usize] + (jmax - l) as f64 * ave_bg[jmin as usize];
        cur_val /= (jmax - jmin) as f64;

        if ave_bg[i] > bg_thres * cur_val {
            peak_location[ct].0 = jmin;
            peak_location[ct].1 = jmax;
            ct += 1;
        }
    }
    for i in 0..ct as usize {
        let jmin = peak_location[i].0;
        let jmax = peak_location[i].1;
        for j in jmin..jmax + 1 {
            ave_bg[j as usize] = -1.0;
            ave_bge[j as usize] = -1.0;
        }
    }
    if n_qexc > 0 {
        for i in 0..n_qexc {
            // Excluded rings - only set flag in
            // ave_bge, to reject peaks but not change radial_bg data.
            for j in qexc[i as usize][0]..qexc[i as usize][1] + 1 {
                ave_bge[j as usize] = -1.0;
            }
        }
    }
}

fn load_cbf(filepath: &str, frame_info: &mut FrameInfo) -> Result<Vec<f64>, Box<dyn Error>> {
    // let cbf: Box<dyn Frame> = frame::open(filepath).expect("Could not open CBF file {filepath}");
    let cbf = Cbf::read_file(filepath).unwrap();
    frame_info.get_frame_info(&cbf)?;
    let framedata = cbf.array().data();
    Ok(framedata.to_vec())
}

fn ad_hoc_thres(ave_val: f64) -> i32 {
    let cdf_thres = 0.99999;
    let mut count_thres: i32 = 0;
    let mut log_fact: f64 = 0.0;
    let log_val: f64 = ave_val.ln();
    let mut cdf: f64 = (-ave_val).exp();
    while cdf < cdf_thres {
        count_thres += 1;
        log_fact += (count_thres as f64).ln();
        let logP = -ave_val + ((count_thres as f64) * log_val) - log_fact;
        cdf += logP.exp();
    }
    count_thres
}
