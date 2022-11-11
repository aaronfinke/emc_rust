use crate::config::configure::Config;
use clap::Parser;
use hdf5::File as H5File;
use make_Emat::{Args, Pix};
use ndarray::Array1;
use rayon::iter;
use rayon::prelude::*;
use std::env;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Seek;
use std::io::SeekFrom;
use std::io::Write;
use std::ops::Range;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use std::sync::{Arc,Mutex};

pub mod config;
const MEM_STEP_SIZE: i32 = 10000;
const ROOT: i32 = 0;

#[derive(Clone, Copy, Debug)]
pub struct Quat(f64, f64, f64, f64, f64);

struct _box {
    row_min: i32,
    row_max: i32,
    col_min: i32,
    col_max: i32,
}
fn main() -> Result<(), Box<dyn Error>> {
    let t1 = Instant::now();

    let args = Args::parse();
    let config_file: &str = args.config.to_str().unwrap();
    let if_low: bool = args.low;

    /* Setup section */
    let config = match Config::new(&config_file) {
        Ok(config) => config,
        Err(why) => return Err(format!("Error in config.ini: {:?}", why).into()),
    };

    match if_low {
        true => println!("make Ematrix for low-resolution reconstruction"),
        false => println!("make Ematrix for local update scheme"),
    };

    let mut quat_file_buff = match File::open(config.quat_file) {
        Err(why) => panic!("Couldn't open quatfile : {}", why),
        Ok(file) => BufReader::new(file),
    };
    quat_file_buff.seek(SeekFrom::Start(0))?;
    let mut buf = vec![0; 4];
    quat_file_buff.read_exact(&mut buf)?;
    let num_rot = u32::from_le_bytes(buf.try_into().unwrap());
    println!("num_rot = {num_rot}");
    let mut buf = vec![0; 4];
    /* read four byte buffer, should be zero */
    quat_file_buff.read_exact(&mut buf)?;
    let buffer = u32::from_le_bytes(buf.try_into().unwrap());
    if buffer != 0 {
        return Err("Error in quat file".into());
    }

    let mut quat: Vec<Quat> = vec![Quat(0., 0., 0., 0., 0.); num_rot as usize];
    let mut buf = vec![0_u8; num_rot as usize * 40];
    quat_file_buff.read_exact(&mut buf)?;
    for r in 0..num_rot as usize {
        let t: usize = r * 40;
        quat[r].0 = f64::from_le_bytes((&buf[t..t + 8]).try_into().unwrap());
        quat[r].1 = f64::from_le_bytes((&buf[t + 8..t + 16]).try_into().unwrap());
        quat[r].2 = f64::from_le_bytes((&buf[t + 16..t + 24]).try_into().unwrap());
        quat[r].3 = f64::from_le_bytes((&buf[t + 24..t + 32]).try_into().unwrap());
        quat[r].4 = f64::from_le_bytes((&buf[t + 32..t + 40]).try_into().unwrap());
    }

    let beam_vec = config.beam_vec;
    let norm: f64 = (beam_vec[0].powi(2) + beam_vec[1].powi(2) + beam_vec[2].powi(2)).sqrt();
    let D: f64 = (config.detd / config.px) * norm / beam_vec[2].abs();
    let rescale: f64 = 1_f64 / config.wl / D;
    let qmax: i32 = ((config.wl * D / config.res_cutoff).ceil()) as i32;
    let qmin: i32 = (2. * D * ((config.Rstop / D).atan() / 2.).sin()) as i32;
    let rmax2: f64 = (qmax as f64 * rescale).powi(2);
    let rmin2: f64 = (qmin as f64 * rescale).powi(2);

    let (qmax, qmin, num_pix, pix) = {
        let rvec_file_buff = H5File::open(&config.rec_vectors_file)
            .expect(&format!(
                "Could not open rvec_file: {}",
                config.rec_vectors_file
            ))
            .group("data")
            .unwrap();
        let qmax: i32 = rvec_file_buff
            .attr("qmax")
            .expect("Could not extract qmax from rvec_file")
            .read_scalar()
            .unwrap();
        let qmin: i32 = rvec_file_buff
            .attr("qmin")
            .expect("Could not extract qmin from rvec_file")
            .read_scalar()
            .unwrap();
        let num_pix: i32 = rvec_file_buff
            .attr("Npix")
            .expect("Could not extract Npix from rvec_file")
            .read_scalar()
            .unwrap();
        println!("num_pix= {num_pix}");
        let qx: Array1<f64> = rvec_file_buff.dataset("qx")?.read_1d().unwrap();
        let qy: Array1<f64> = rvec_file_buff.dataset("qy")?.read_1d().unwrap();
        let qz: Array1<f64> = rvec_file_buff.dataset("qz")?.read_1d().unwrap();
        let pix: Vec<Pix> = {
            let mut pix: Vec<Pix> = Vec::new();
            for t in 0..num_pix as usize {
                let mut norm: f64 = 0_f64;
                pix.push(Pix(qx[t] * rescale, qy[t] * rescale, qz[t] * rescale));
                norm += pix[t].0 * pix[t].0;
                norm += pix[t].1 * pix[t].1;
                norm += pix[t].2 * pix[t].2;
                if norm > rmax2 {
                    break;
                }
            }
            pix
        };
        let num_pix = pix.len();
        (qmax, qmin, num_pix, pix)
    };

    println!("num_pix = {num_pix}");

    let rec2pixmap_path = Path::new(&config.rec2pix_file);
    let rec2pix_map: Vec<i32> = match rec2pixmap_path.extension().and_then(OsStr::to_str) {
        Some("dat") => Config::load_pixmap_dat(rec2pixmap_path, num_pix as usize)?,
        Some("h5") => Config::load_pixmap_h5(rec2pixmap_path, "rec2pix_map")?,
        _ => {
            return Err(format!(
                "could not determine pixmap type\nshould be .h5 or .dat (if text file)"
            )
            .into())
        }
    };
    let total_pix: usize = (config.num_row * config.num_col) as usize;
    let pixmap_path = Path::new(&config.pixmap_file);
    let pix_map: Vec<i32> = match pixmap_path.extension().and_then(OsStr::to_str) {
        Some("dat") => Config::load_pixmap_dat(pixmap_path, total_pix)?,
        Some("h5") => Config::load_pixmap_h5(pixmap_path, "pix_map")?,
        _ => {
            return Err(format!(
                "could not determine pixmap type\nshould be .h5 or .dat (if text file)"
            )
            .into())
        }
    };

    let total_pix: i32 = config.num_col * config.num_row;

    let mut unit_cell_buff = match File::open(&config.unit_cell_file) {
        Err(why) => panic!("Couldn't open unit cell file : {}", why),
        Ok(file) => BufReader::new(file),
    };
    let mut first_line = String::new();
    unit_cell_buff.read_line(&mut first_line)?;
    let mut ucline: Vec<&str> = first_line.trim().split_whitespace().collect();
    ucline.pop();
    let mut cell: Vec<f64> = ucline
        .iter()
        .map(|word| word.parse::<f64>().unwrap())
        .collect();
    match cell.len() == 6 {
        True => {}
        False => return Err("Error in unit_cell file".into()),
    }
    let degtorad = 1_f64.atan() * 4.0 / 180.0;
    for i in 3..6 {
        cell[i] *= degtorad;
    }
    let (sina, cosa) = (cell[3].sin(), cell[3].cos());
    let (sinb, cosb) = (cell[4].sin(), cell[4].cos());
    let (sing, cosg) = (cell[5].sin(), cell[5].cos());
    let cosas = (cosg * cosb - cosa) / (sinb * sing);
    let sinas = (1.0 - cosas * cosas).sqrt();
    let cosbs = (cosa * cosg - cosb) / (sina * sing);
    let sinbs = (1.0 - cosbs * cosbs).sqrt();
    let cosgs = (cosa * cosb - cosg) / (sina * sinb);
    let sings = (1.0 - cosgs * cosgs).sqrt();

    let sumc = (cell[3] + cell[4] + cell[5]) * 0.5;
    let vol = 2.0
        * cell[0]
        * cell[1]
        * cell[2]
        * ((sumc - cell[3]).sin() * (sumc - cell[4]).sin() * (sumc - cell[5]).sin() * (sumc).sin())
            .sqrt();

    let rcell: [f64; 6] = [
        cell[1] * cell[2] * sina / vol,
        cell[2] * cell[0] * sinb / vol,
        cell[0] * cell[1] * sing / vol,
        (sinas).atan2(cosas),
        (sinbs).atan2(cosbs),
        (sings).atan2(cosgs),
    ];
    let rlatt_vec: [[f64; 3]; 3] = [
        [rcell[0], 0.0, 0.0],
        [rcell[1] * cosgs, rcell[1] * sings, 0.0],
        [
            rcell[2] * cosbs,
            -rcell[2] * sinbs * cosa,
            rcell[2] * sinbs * sina,
        ],
    ];

    for i in 0..3 {
        println!(
            "Reciprocal lattice vector {}:\n\t{:.3e} {:.3e} {:.3e}",
            i, rlatt_vec[i][0], rlatt_vec[i][1], rlatt_vec[i][2]
        );
    }
    let mut norm = 0_f64;
    let mut tol = 0_f64;

    for j in 0..3 {
        norm = 0_f64;
        for i in 0..3 {
            norm += rlatt_vec[j][i].powi(2);
        }
        if j == 0 {
            tol = norm;
        } else {
            if tol > norm {
                tol = norm;
            }
        }
    }
    let min_rcell_sz = tol.sqrt();
    tol *= (config.gw / config.VN as f64).powi(2);

    /* calculate cofactors of rlatt_vec */
    let mut cofactor: [[f64; 3]; 3] = [[0_f64; 3]; 3];
    let mut mat: [[f64; 3]; 3] = [[0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cofactor[i][j] = match (i + j) % 2 {
                0 => 1_f64,
                _ => -1_f64,
            };
            let mut i1 = 0;
            for i0 in 0..3 {
                if i0 == i {
                    continue;
                };
                let mut j1 = 0;
                for j0 in 0..3 {
                    if j0 == j {
                        continue;
                    };
                    mat[i1][j1] = rlatt_vec[i0][j0];
                    j1 += 1;
                }
                i1 += 1;
            }

            let det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            cofactor[i][j] *= det;
        }
    }

    /* determinant of rlatt_vec */
    let mut det = 0_f64;
    for i in 0..3 {
        det += rlatt_vec[0][i] * cofactor[0][i];
    }
    /* inverse of rlatt_vec */
    let mut proj: [[f64; 3]; 3] = [[0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            proj[i][j] = cofactor[j][i] / det;
        }
    }

    /* determine hmax, kmax, lmax */
    let mut norm = 0_f64;
    for i in 0..3 {
        norm += rlatt_vec[0][i].powi(2);
    }
    let mut hcur: i32 = (rmax2 / norm).sqrt().floor() as i32;

    norm = 0_f64;
    for i in 0..3 {
        norm += rlatt_vec[1][i].powi(2);
    }
    let mut kcur: i32 = (rmax2 / norm).sqrt().floor() as i32;

    norm = 0_f64;
    for i in 0..3 {
        norm += rlatt_vec[2][i].powi(2);
    }
    let mut lcur: i32 = (rmax2 / norm).sqrt().floor() as i32;

    let mut indices: [i32; 3] = [0, 0, 0];
    let mut rvec: [f64; 3] = [0_f64, 0_f64, 0_f64];
    let mut if_enclosed = true;

    while if_enclosed {
        if_enclosed = false;
        for i in 0..2 {
            indices[0] = (2 * i - 1) * hcur;
            for j in -kcur..kcur + 1 {
                indices[1] = j;
                for k in -lcur..lcur + 1 {
                    indices[2] = k;
                    if if_enclosed {
                        break;
                    }
                    let mut norm = 0.0;
                    for i0 in 0..3 {
                        rvec[i0] = 0.0;
                        for j0 in 0..3 {
                            rvec[i0] += rlatt_vec[i0][j0] * indices[j0] as f64;
                        }
                        norm += rvec[i0] * rvec[i0];
                    }
                    if norm < rmax2 {
                        if_enclosed = true;
                    }
                }
            }
        }

        if if_enclosed {
            hcur += 1;
            continue;
        }
        for i in -hcur..hcur + 1 {
            indices[0] = i;
            for j in 0..2 {
                indices[1] = (2 * j - 1) * kcur;
                for k in -lcur..lcur + 1 {
                    indices[2] = k;
                    if if_enclosed {
                        break;
                    }
                    let mut norm = 0_f64;
                    for i0 in 0..3 {
                        rvec[i0] = 0_f64;
                        for j0 in 0..3 {
                            rvec[i0] += rlatt_vec[i0][j0] * indices[j0] as f64;
                        }
                        norm += rvec[i0] * rvec[i0];
                    }
                    if norm < rmax2 {
                        if_enclosed = true;
                    }
                }
            }
        }

        if if_enclosed {
            kcur += 1;
            continue;
        }

        for i in -hcur..hcur + 1 {
            indices[0] = i;
            for j in -kcur..kcur + 1 {
                indices[1] = j;
                for k in 0..2 {
                    indices[2] = (2 * k - 1) * lcur;
                    if if_enclosed {
                        break;
                    }
                    let mut norm = 0_f64;
                    for i0 in 0..3 {
                        rvec[i0] = 0_f64;
                        for j0 in 0..3 {
                            rvec[i0] += rlatt_vec[i0][j0] * indices[j0] as f64;
                        }
                        norm += rvec[i0] * rvec[i0];
                    }
                    if norm < rmax2 {
                        if_enclosed = true;
                    }
                }
            }
        }

        if if_enclosed {
            lcur += 1;
            continue;
        }
    }

    let hmax = hcur;
    let kmax = kcur;
    let lmax = lcur;
    let hlen = 2 * hmax + 1;
    let klen = 2 * kmax + 1;
    let llen = 2 * lmax + 1;

    let max_num_hkl: i32 = hlen * klen * llen;
    let mut is_partial: Vec<i32> = vec![0; max_num_hkl as usize];
    let mut hkl_map: Vec<i32> = vec![-1; max_num_hkl as usize];

    let rescale: f64 = config.VN as f64 / min_rcell_sz;
    let gw_ceil: i32 = config.gw.ceil() as i32;
    let gw2 = config.gw * config.gw;

    let mut num_hkl: usize = 0;
    for i in -hmax..hmax + 1 {
        for j in -kmax..kmax + 1 {
            for k in -lmax..lmax + 1 {
                indices = [i, j, k];
                let idx = (i + hmax) * klen * llen + (j + kmax) * llen + (k + lmax);
                for i0 in 0..3 {
                    rvec[i0] = 0_f64;
                    for j0 in 0..3 {
                        rvec[i0] += rlatt_vec[i0][j0] * (indices[j0] as f64);
                    }
                }
                let vx = rvec[0] * rescale;
                let vy = rvec[1] * rescale;
                let vz = rvec[2] * rescale;

                let tx = vx.round() as i32;
                let ty = vy.round() as i32;
                let tz = vz.round() as i32;

                for x in -(gw_ceil + 1)..(gw_ceil + 2) {
                    for y in -(gw_ceil + 1)..(gw_ceil + 2) {
                        for z in -(gw_ceil + 1)..(gw_ceil + 2) {
                            if is_partial[idx as usize] == 1 {
                                continue;
                            }
                            let dx = (tx + x) as f64 - vx;
                            let dy = (ty + y) as f64 - vy;
                            let dz = (tz + z) as f64 - vz;
                            let dr2 = dx * dx + dy * dy + dz * dz;

                            if dr2 <= gw2 {
                                norm = ((tx + x) * (tx + x)
                                    + (ty + y) * (ty + y)
                                    + (tz + z) * (tz + z))
                                    as f64;
                                norm /= rescale * rescale;
                                if norm > rmax2 || norm < rmin2 {
                                    is_partial[idx as usize] = 1;
                                    continue;
                                }
                            }
                        }
                    }
                }
                if is_partial[idx as usize] == 0 {
                    hkl_map[idx as usize] = num_hkl as i32;
                    num_hkl += 1;
                }
            }
        }
    }

    println!("hmax = {hmax}, kmax = {kmax}, lmax = {lmax}, num_hkl = {num_hkl}");
    let mut if_skip_r: Vec<bool> = vec![true; num_rot as usize];

    if if_low {
        let num_data: i32 = {
            let mut fp = match File::open(&config.reduced_data_id_file) {
                Err(why) => panic!("Couldn't open reduced_data_id_file : {}", why),
                Ok(file) => BufReader::new(file).lines(),
            };
            fp.next()
                .expect("Could not read reduced_data_id_file")
                .unwrap()
                .parse::<i32>()
                .expect("Could not extract num_data from reduced_data_id_file")
        };
        {
            let mut prob_orien_file_buff = match File::open(&config.prob_orien_file) {
                Err(why) => panic!("Couldn't open prob_orien_file : {}", why),
                Ok(file) => BufReader::new(file),
            };
            prob_orien_file_buff.seek(SeekFrom::Start(0))?;
            for d in 0..num_data {
                let mut buf: Vec<u8> = vec![0; 4];
                prob_orien_file_buff.read_exact(&mut buf)?;
                let num: i32 = i32::from_le_bytes(buf.try_into().unwrap());
                for i in 0..num {
                    let mut buf: Vec<u8> = vec![0; 4];
                    prob_orien_file_buff.read_exact(&mut buf)?;
                    let r: i32 = i32::from_le_bytes(buf.try_into().unwrap());
                    if_skip_r[r as usize] = false;
                }
            }
        }
        let mut num = 0;
        for r in 0..num_rot as usize {
            if !if_skip_r[r] {
                num += 1;
            }
        }
        println!("number of rotations = {num_rot}");
        println!("number of relevant rotations = {num}");
    }
    let t2 = t1.elapsed();
    let mut r2peak: Vec<Vec<i32>> = vec![Vec::default(); num_rot as usize];
    let mut cur_r2peak: Vec<i32> = vec![0; num_hkl];
    let mut num_r2peak: Vec<i32> = vec![0; num_rot as usize];
    let mut num_peak2r: Vec<i32> = vec![0; num_rot as usize];
    let mut peak_r: Vec<_box> = Vec::with_capacity(num_hkl);
    let mut if_visit: Vec<bool> = vec![false; num_hkl];

    println!("Setup takes: {:?}", t2);
    /* end setup section */

    /* start make_EMat section */

    Ok(())
}

struct EMatResults {
    r2peak: Vec<Vec<i32>>,
    cur_r2peak: Vec<i32>,
    num_r2peak: Vec<i32>,
    num_peak2r: Vec<i32>,
    peak_r: Vec<_box>,
    if_visit: Vec<bool>,
}

fn make_Emat(
    r: usize,
    config: Config,
    quat: &Vec<Quat>,
    klen: i32,
    llen: i32,
    hmax: i32,
    kmax: i32,
    lmax: i32,
    num_rot: i32,
    num_pix: i32,
    num_col: i32,
    pix: Vec<Pix>,
    proj: [[f64; 3]; 3],
    hkl_map: Vec<i32>,
    rlatt_vec: [[f64; 3]; 3],
    tol: f64,
    rec2pix_map: Vec<i32>,
    if_skip_r: Vec<bool>,
    ) -> Result<(), Box<dyn Error>> {
    let (mut rmin, mut rmax, mut peak_ct): (i32, i32, i32);
    let (klen_llen, mut hkl_id): (i32, usize);
    let (mut idx, mut row_id, mut col_id, mut length): (i32, i32, i32, i32);
    let mut rot_pix: [f64; 3] = [0_f64; 3];
    let mut coefft: [f64; 3] = [0_f64; 3];
    let mut alloc_r2p;

    let load_id: i32 = 1;

    klen_llen = klen * llen;
    rmin = 0;
    rmax = num_rot;

    // rmin = load_id*MEM_STEP_SIZE;
    // if rmin > num_rot-1 {return Ok(())};
    // rmax = match rmin + MEM_STEP_SIZE {
    //     x if x > num_rot => num_rot,
    //     _ => rmin + MEM_STEP_SIZE
    // };

    // load2rank[load_id] = myid;
    // println!("myid= {myid}, load_id={load_id},(rmin,rmax) = ({rmin},{rmax}");
    alloc_r2p = 0;
    // let rmin = 0;
    // let rmax = num_rot;

        let rot = make_rot(r as usize, quat);
        peak_ct = 0;
        if if_skip_r[r] {
            num_r2peak[r] = peak_ct;
            r2peak[r] = vec![0; (1 + 3 * peak_ct) as usize];
            alloc_r2p += 1 + 3 * peak_ct;
            r2peak[r][0] = peak_ct;
            continue;
        }
        for t in 0..num_pix as usize {
            for i in 0..3 {
                rot_pix[i] = rot[i][0] * pix[t].0;
                rot_pix[i] += rot[i][1] * pix[t].1;
                rot_pix[i] += rot[i][2] * pix[t].2;
            }
            for i in 0..3 {
                coefft[i] = proj[i][0] * rot_pix[0];
                for j in 1..3 {
                    coefft[i] += proj[i][j] * rot_pix[j];
                }
            }

            let hval = coefft[0].round() as i32;
            if hval.abs() > hmax {
                continue;
            }
            let kval = coefft[1].round() as i32;
            if kval.abs() > kmax {
                continue;
            }
            let lval = coefft[2].round() as i32;
            if lval.abs() > lmax {
                continue;
            }

            idx = (hval + hmax) * klen_llen + (kval + kmax) * llen + (lval + lmax);
            hkl_id = hkl_map[idx as usize] as usize;
            if hkl_id < 0 {
                continue;
            }

            let dx = rot_pix[0]
                - (hval as f64) * rlatt_vec[0][0]
                - (kval as f64) * rlatt_vec[0][1]
                - (lval as f64) * rlatt_vec[0][2];
            let dy = rot_pix[1]
                - (hval as f64) * rlatt_vec[1][0]
                - (kval as f64) * rlatt_vec[1][1]
                - (lval as f64) * rlatt_vec[1][2];
            let dz = rot_pix[2]
                - (hval as f64) * rlatt_vec[2][0]
                - (kval as f64) * rlatt_vec[2][1]
                - (lval as f64) * rlatt_vec[2][2];

            /* rot_pix is out of the view of its closet Bragg peak */
            if dx * dx + dy * dy + dz * dz > tol {
                continue;
            }

            row_id = rec2pix_map[t] / num_col;
            col_id = rec2pix_map[t] % num_col;

            if if_visit[hkl_id] == false {
                cur_r2peak[peak_ct as usize] = hkl_id as i32;
                peak_ct += 1;
                peak_r[hkl_id].row_min = row_id;
                peak_r[hkl_id].row_max = row_id;
                peak_r[hkl_id].col_min = col_id;
                peak_r[hkl_id].col_max = col_id;
                if_visit[hkl_id] = true;
            } else {
                if (peak_r[hkl_id].row_min > row_id) {
                    peak_r[hkl_id].row_min = row_id;
                }
                if (peak_r[hkl_id].row_max < row_id) {
                    peak_r[hkl_id].row_max = row_id;
                }
                if (peak_r[hkl_id].col_min > col_id) {
                    peak_r[hkl_id].col_min = col_id;
                }
                if (peak_r[hkl_id].col_max < col_id) {
                    peak_r[hkl_id].col_max = col_id;
                }
            }
        }

        for t in 0..peak_ct {
            hkl_id = cur_r2peak[t as usize] as usize;
            num_peak2r[hkl_id] += 1;
            if_visit[hkl_id] = false;
        }

        num_r2peak[r] = peak_ct;
        r2peak[r] = vec![0; (1 + 3 * peak_ct) as usize];
        alloc_r2p += (1 + 3 * peak_ct);
        for t in 0..peak_ct {
            idx = 3 * t + 1;
            hkl_id = cur_r2peak[t as usize] as usize;
            r2peak[r][idx as usize] = hkl_id as i32;
            r2peak[r][(idx + 1) as usize] =
                peak_r[hkl_id].row_min * num_col + peak_r[hkl_id].col_min;
            r2peak[r][(idx + 2) as usize] =
                peak_r[hkl_id].row_max * num_col + peak_r[hkl_id].col_max;
        }
    

    Ok(())
}

fn make_rot(r: usize, quat: &Vec<Quat>) -> [[f64; 3]; 3] {
    let mut rot: [[f64; 3]; 3] = [[0_f64; 3]; 3];

    let q0 = quat[r].0;
    let q1 = quat[r].1;
    let q2 = quat[r].2;
    let q3 = quat[r].3;
    let q01 = q0 * q1;
    let q02 = q0 * q2;
    let q03 = q0 * q3;
    let q11 = q1 * q1;
    let q12 = q1 * q2;
    let q13 = q1 * q3;
    let q22 = q2 * q2;
    let q23 = q2 * q3;
    let q33 = q3 * q3;

    rot[0][0] = 1. - 2. * (q22 + q33);
    rot[0][1] = 2. * (q12 + q03);
    rot[0][2] = 2. * (q13 - q02);
    rot[1][0] = 2. * (q12 - q03);
    rot[1][1] = 1. - 2. * (q11 + q33);
    rot[1][2] = 2. * (q01 + q23);
    rot[2][0] = 2. * (q02 + q13);
    rot[2][1] = 2. * (q23 - q01);
    rot[2][2] = 1. - 2. * (q11 + q22);

    rot
}
