/*
exhaustive search for the compatible orientations of the
identified peaks using Rayon


needs:
config file, quat file, unit cell file, peaklist file, orienlist file

makes:
orien_files, num_prob_orien.dat (in data_dir)

*/

use libm::*;
use orient_peak::{Config, Results};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::convert::TryInto;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::process;
use std::sync::{Arc, Mutex};
use std::time::Instant;
const MEM_STEP_SIZE: u32 = 100;
const MAX_ORIEN: u32 = 2000;

fn main() -> Result<(), Box<dyn Error>> {
    let t1 = Instant::now();

    let args: Vec<String> = env::args().collect();
    match args.len() {
        2 => {}
        _ => return Err("execution: ./make_background (/path/to/)config.ini".into()),
    };

    let mut config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {err}");
        process::exit(1);
    });
    (config.peak_rvec, config.num_peak) = config.get_num_peak().unwrap();
    let elapsed = t1.elapsed();
    rayon::ThreadPoolBuilder::new().num_threads(config.nproc as usize).build_global().unwrap();
    println!("Setup took {:?}", elapsed);

    let all_results= Arc::new(Mutex::new(vec![]));
    /* End of setup section */

    let _num_load: usize = match config.num_data % MEM_STEP_SIZE as usize {
        0 => config.num_data / MEM_STEP_SIZE as usize,
        _ => config.num_data / (MEM_STEP_SIZE + 1) as usize,
    };

    let now = Instant::now();
    (0.._num_load).into_par_iter().for_each(|i| {
        if let Ok(res) = rotate_rvec(i.try_into().unwrap(), &config) {
        all_results.lock().unwrap().push(res);
    }});
    let elapsed = now.elapsed();
    println!("rotate_rvec took {:?}", elapsed);
    let results: Vec<Results> = all_results.lock().unwrap().clone().into_iter().flatten().collect();
    /* end rotate_rvec section */
    let now = Instant::now();
    let mut d_np_orien: Vec<(i32,u32)> = vec![];
    for mut result in results {
        d_np_orien.push((result.d, result.num_prob_orien));
        let mut zipped = result
            .prob_orien
            .clone()
            .into_iter()
            .zip(result.score.clone().into_iter())
            .collect::<Vec<(i32, f64)>>();
        zipped.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).expect("Can't compare"));
        (
            result.prob_orien,
            result.score,
        ) = zipped.into_iter().unzip();
        {
            let orienfile = PathBuf::from(&config.orient_files[result.d as usize]);
            let mut fp = File::create(orienfile)?;
            let mut k = result.num_prob_orien;
            if k > MAX_ORIEN {
                k = MAX_ORIEN;
            }
            writeln!(&mut fp, "{k}")?;
            for i in 0..k as usize {
                let f: i32 = result.prob_orien[i];
                writeln!(&mut fp, "{f}",)?;
            }
        }
    }
    println!("End mergesort");
    let elapsed = now.elapsed();
    println!("mergesort took {:?}", elapsed);
    let mut ct = 0;
    d_np_orien.sort_by(|a,b| a.0.partial_cmp(&b.0).expect("problem sorting d, num_prob_orien"));
    {
        let fname: PathBuf = [&config.data_dir, "num_prob_orien.dat"].iter().collect();
        let mut fp = File::create(fname)?;
        for d in 0..config.num_data {
            let mut k = d_np_orien[d].1;
            if k > MAX_ORIEN {
                k = MAX_ORIEN;
            }
            writeln!(&mut fp, "{k}")?;
            if d_np_orien[d].1 > 0 {
                ct += 1;
            }
        }
    }
    println!("Number of relevant frames: = {ct}");

    let t2 = t1.elapsed();
    println!("Elapsed time: {:?}", t2);

    Ok(())
}

fn rotate_rvec(load_id: u32, config: &Config) -> Result<Vec<Results>, Box<dyn Error>> {
    /* the rotate_rvec section */
    let mut resulty: Vec<Results> = vec![];
    let klen_llen = config.klen * config.llen;
    let tol2 = config.tolmult * config.tol;
    println!("Tolerance used: {tol2}");

    let dmin:usize = (load_id * MEM_STEP_SIZE) as usize;
        if dmin > config.num_data {
            return Err("it's fine".into());
        }
    let mut dmax:usize  = ((load_id+1) * MEM_STEP_SIZE) as usize;
    if dmax > config.num_data {
        dmax = config.num_data;
    }
    for d in dmin..(dmax) {
        resulty.push(Results::new(d as i32).unwrap());
    }
    println!("dmin = {dmin}, dmax = {dmax}");

    for r in 0..config.num_rot as usize {
        let mut hkl_list: Vec<i32> = vec![Default::default(); config.num_hkl];
        let mut if_visit: Vec<bool> = vec![false; config.num_hkl];
        let mut hdiffs: Vec<f64> = vec![0_f64; (config.max_num_peak * 3) as usize];

        let mut rvec: [f64; 3] = [0_f64; 3];
        let mut weighted_vec: [f64; 3] = [0_f64; 3];
        let mut coefft: [f64; 3] = [0_f64; 3];
        let mut rvec1: [f64; 3] = [0_f64; 3];
        let mut rvec2: [f64; 3] = [0_f64; 3];

        let mut rot: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        let mut t: usize = 0;
        for i in 0..3 {
            for j in 0..3 {
                rot[i][j] = config.rot_mat_table[r][t];
                t += 1;
            }
        }
        for result in &mut resulty {
            let d: usize = result.d as usize;
            if config.num_peak[d] >= config.min_num_peak {

            let mut hkl_ct: u32 = 0;
            // let mut accept: i32 = 0;
            let mut score = 0.0;
            for t in 0..config.num_peak[d] as usize {
                for j in 0..3 {
                    weighted_vec[j] = config.peak_rvec[d][t].vec[j];
                }
                for j in 0..3 {
                    rvec[j] = rot[j][0] * weighted_vec[0];
                    for k in 1..3 {
                        rvec[j] += rot[j][k] * weighted_vec[k];
                    }
                }
                for j in 0..3 {
                    coefft[j] = 0.;
                    for k in 0..3 {
                        coefft[j] += config.proj[j][k] * rvec[k];
                    }
                }

                let hval = round(coefft[0]) as i32;
                if (hval).abs() > config.hmax {
                    continue;
                }
                let kval = round(coefft[1]) as i32;
                if (kval).abs() > config.kmax {
                    continue;
                }
                let lval = round(coefft[2]) as i32;
                if (lval).abs() > config.lmax {
                    continue;
                }

                let dx: f64 = rvec[0]
                    - hval as f64 * config.rlatt_vec[0][0]
                    - kval as f64 * config.rlatt_vec[0][1]
                    - lval as f64 * config.rlatt_vec[0][2];
                let dy: f64 = rvec[1]
                    - hval as f64 * config.rlatt_vec[1][0]
                    - kval as f64 * config.rlatt_vec[1][1]
                    - lval as f64 * config.rlatt_vec[1][2];
                let dz: f64 = rvec[2]
                    - hval as f64 * config.rlatt_vec[2][0]
                    - kval as f64 * config.rlatt_vec[2][1]
                    - lval as f64 * config.rlatt_vec[2][2];
                let dxdif: f64 = dx * dx + dy * dy + dz * dz;
                /*
                 * If images are oscillations, or there is a wide bandwidth, also calculate
                 * deviations with small +/- rotations applied, or with small +/-
                 * increments to wavelength.
                 */
                let (dxdif1, dxdif2) = match config.delphi.partial_cmp(&0.001) {
                    Some(Ordering::Greater) => {
                        /* Get vector assuming rotation of crystal by small angle delphi, about
                         * spindle axis
                         */
                        weighted_vec[1] = config.peak_rvec[d][t].vec[1]
                            + config.peak_rvec[d][t].vec[2] * config.tandp2;
                        weighted_vec[2] = config.peak_rvec[d][t].vec[2]
                            - config.peak_rvec[d][t].vec[1] * config.tandp2;

                        /* Alternatively, get vector assuming small shift in wavelength
                         */
                        for j in 0..3 {
                            rvec1[j] = rot[j][0] * weighted_vec[0];
                            for k in 1..3 {
                                rvec1[j] += rot[j][k] * weighted_vec[k];
                            }
                        }
                        for j in 0..3 {
                            coefft[j] = 0_f64;
                            for k in 0..3 {
                                coefft[j] += config.proj[j][k] * rvec1[k];
                            }
                        }
                        let hval1: i32 = round(coefft[0]) as i32;
                        let kval1: i32 = round(coefft[1]) as i32;
                        let lval1: i32 = round(coefft[2]) as i32;
                        let dx1: f64 = rvec1[0]
                            - hval1 as f64 * config.rlatt_vec[0][0]
                            - kval1 as f64 * config.rlatt_vec[0][1]
                            - lval1 as f64 * config.rlatt_vec[0][2];
                        let dy1: f64 = rvec1[1]
                            - hval1 as f64 * config.rlatt_vec[1][0]
                            - kval1 as f64 * config.rlatt_vec[1][1]
                            - lval1 as f64 * config.rlatt_vec[1][2];
                        let dz1: f64 = rvec1[2]
                            - hval1 as f64 * config.rlatt_vec[2][0]
                            - kval1 as f64 * config.rlatt_vec[2][1]
                            - lval1 as f64 * config.rlatt_vec[2][2];
                        let dxdif1: f64 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
                        /* Get vector for offset in opposite direction */
                        /* For rotation of crystal by small angle delphi, about spindle axis */
                        weighted_vec[1] = config.peak_rvec[d][t].vec[1]
                            - config.peak_rvec[d][t].vec[2] * config.tandp2;
                        weighted_vec[2] = config.peak_rvec[d][t].vec[2]
                            + config.peak_rvec[d][t].vec[1] * config.tandp2;
                        /* For small shift in wavelength */
                        for j in 0..3 {
                            rvec2[j] = rot[j][0] * weighted_vec[0];
                            for k in 1..3 {
                                rvec2[j] += rot[j][k] * weighted_vec[k];
                            }
                        }
                        for j in 0..3 {
                            coefft[j] = 0.;
                            for k in 0..3 {
                                coefft[j] += config.proj[j][k] * rvec2[k];
                            }
                        }
                        let hval2: i32 = round(coefft[0]) as i32;
                        let kval2: i32 = round(coefft[1]) as i32;
                        let lval2: i32 = round(coefft[2]) as i32;
                        let dx2: f64 = rvec2[0]
                            - hval2 as f64 * config.rlatt_vec[0][0]
                            - kval2 as f64 * config.rlatt_vec[0][1]
                            - lval2 as f64 * config.rlatt_vec[0][2];
                        let dy2: f64 = rvec2[1]
                            - hval2 as f64 * config.rlatt_vec[1][0]
                            - kval2 as f64 * config.rlatt_vec[1][1]
                            - lval2 as f64 * config.rlatt_vec[1][2];
                        let dz2: f64 = rvec2[2]
                            - hval2 as f64 * config.rlatt_vec[2][0]
                            - kval2 as f64 * config.rlatt_vec[2][1]
                            - lval2 as f64 * config.rlatt_vec[2][2];
                        let dxdif2: f64 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
                        (dxdif1, dxdif2)
                    }
                    Some(Ordering::Equal) => (0.0, 0.0),
                    Some(Ordering::Less) => (0.0, 0.0),
                    None => (0.0, 0.0),
                };
                /* rvec is too far from closest Bragg peak */
                let mut dxd = dxdif;
                if dxdif > tol2 {
                    if config.delphi < 0.001 {
                        continue;
                    } else {
                        if dxdif1 > tol2 {
                            if dxdif2 > tol2 {
                                continue;
                            }
                        }
                        if dxdif1 < dxdif {
                            dxd = dxdif1;
                        }
                        if dxdif2 < dxd {
                            dxd = dxdif2;
                        }
                    }
                }
                let k: usize = 3 * hkl_ct as usize;
                hdiffs[k] = dxdif;
                if config.delphi > 0.001 {
                    hdiffs[k + 1 as usize] = dxdif1;
                    hdiffs[k + 2] = dxdif2;
                }
                let hkl_id = config.hkl_map[((hval + config.hmax) * klen_llen
                    + (kval + config.kmax) * config.llen
                    + (lval + config.lmax)) as usize];
                if hkl_id < 0 {
                    continue;
                }
                if if_visit[hkl_id as usize] == true {
                    continue;
                } else {
                    if_visit[hkl_id as usize] = true;
                    hkl_list[hkl_ct as usize] = hkl_id;
                    hkl_ct += 1;
                    score += dxd / tol2;
                }
                /* Accept case with at least min_num_peak matches, but don't stop looking
                 * for matches when this number is reached
                 */
                // if hkl_ct == config.min_num_peak as u32 {
                //     // accept = 1;
                // }
            }
            for i in 0..hkl_ct as usize {
                if_visit[hkl_list[i] as usize] = false;
            }
            /* Not enough peaks match the reciprocal lattice in this orientation */

            if hkl_ct >= config.min_num_peak.try_into().unwrap() {
            /*
             * Calculate score for the orientation based on number of matching peaks
             * first, then average discrepancy
             */
            let d1 = hkl_ct as f64;
            score = d1 + 1.0 - (score / d1);
            /* Negate score so an ascending sort will put best values first */
            score = -score;
            /* Store hkl of matched peaks, and rotation matrix, in case we need to dump
             * them for debugging purposes */
            //   for i in 0..hkl_ct {
            //       let k = (results.num_prob_orien[d] * config.max_num_peak) + (3 * i);
            //       let i1 = config.hkl_unmap[hkl_list[i as usize] as usize];
            //       let mut found_n = i1 % config.llen - config.lmax;
            //       results.found_hkl[d][(k + 2) as usize] = found_n as u32;
            //       let i2 = i1 / config.llen;
            //       found_n = i2 % config.klen - config.kmax;
            //       results.found_hkl[d][(k + 1) as usize] = found_n as u32;
            //       found_n = i2 / config.klen - config.hmax;
            //       results.found_hkl[d][k as usize] = found_n as u32;
            //   }
            //   let mut k = results.num_prob_orien[d] * 9;
            //   for i in 0..3 {
            //       for j in 0..3 {
            //         results.rmatx[d][k as usize] = rot[i][j];
            //           k += 1;
            //       }
            //   }
            result.score.push(score);
            result.n_found_hkl.push(hkl_ct);
            result.prob_orien.push(r as i32);
            result.num_prob_orien += 1;
            }}}
    }
    Ok(resulty)
}
