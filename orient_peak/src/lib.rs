use configparser::ini::Ini;
use libm::*;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::SeekFrom;
use std::path::PathBuf;
use std::time::Instant;

const MEM_STEP_SIZE: u32 = 100;
const MAX_ORIEN: u32 = 2000;
#[derive(Default, Clone, Copy)]
pub struct r_vector {
    pub vec: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct Quat(f64, f64, f64, f64, f64);
#[derive(Clone, Copy, Debug)]

pub struct Pix(pub f64, pub f64, pub f64);
pub struct Config {
    pub data_dir: String,
    pub radial_bg_dir: String,
    pub outlier_dir: String,
    pub data_info_file: String,
    pub file_type: String,
    pub home_dir: String,
    pub pixmap_file: String,
    pub rec2pixmap_file: String,
    pub rvec_file: String,
    pub cbflist_file: String,
    pub radiallist_file: String,
    pub outlierlist_file: String,
    pub peaklist_file: String,
    pub orienlist_file: String,
    pub unit_cell_file: String,
    pub num_row: i32,
    pub num_col: i32,
    pub detd: f64,
    pub wl: f64,
    pub px: f64,
    pub cx: f64,
    pub cy: f64,
    pub num_data: usize,
    pub hot_pix_thres: i32,
    pub orient_files: Vec<String>,
    pub peak_files: Vec<String>,
    pub qmax: i32,
    pub qmin: i32,
    pub total_pix: i32,
    pub min_patch_sz: u32,
    pub max_patch_sz: u32,
    pub min_num_peak: u32,
    pub max_num_peak: u32,
    pub res_cutoff: f64,
    pub VN: i32,
    pub gw: f64,
    pub quat: Vec<Quat>,
    pub tolmult: f64,
    pub delphi: f64,
    pub beam_vec: [f64; 3],
    pub rmax2: f64,
    pub rmin2: f64,
    pub pix: Vec<Pix>,
    pub cell: [f64; 6],
    pub rot: [[f64; 3]; 3],
    pub rot_mat_table: Vec<Vec<f64>>,
    pub num_rot: u32,
    pub hmax: i32,
    pub kmax: i32,
    pub lmax: i32,
    pub hlen: i32,
    pub klen: i32,
    pub llen: i32,
    pub hkl_map: Vec<i32>,
    pub hkl_unmap: Vec<i32>,
    pub num_hkl: usize,
    pub rlatt_vec: [[f64; 3]; 3],
    pub proj: [[f64; 3]; 3],
    pub tol: f64,
    pub tandp2: f64,
    pub peak_rvec: Vec<Vec<r_vector>>,
    pub num_peak: Vec<u32>,
    pub nproc: i32,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, Box<dyn Error>> {
        let mut config = Ini::new();
        let config_file = &args[1];
        config.load(config_file)?;

        let data_dir = config.get("files", "data_dir").unwrap();
        let radial_bg_dir = config.get("files", "radial_bg_dir").unwrap();
        let outlier_dir = config.get("files", "outlier_dir").unwrap();
        let data_info_file = config.get("files", "data_info_file").unwrap();
        let file_type = config.get("make_detector", "file_type").unwrap();
        let home_dir = config.get("files", "home_dir").unwrap();
        let pixmap_file = config.get("files", "pixmap_file").unwrap();
        let rec2pixmap_file = config.get("files", "rec2pix_file").unwrap();
        let rvec_file = config.get("files", "rec_vectors_file").unwrap();
        let cbflist_file = config.get("files", "cbflist_file").unwrap();
        let radiallist_file = config.get("files", "radiallist_file").unwrap();
        let outlierlist_file = config.get("files", "outlierlist_file").unwrap();
        let peaklist_file = config.get("files", "peaklist_file").unwrap();
        let orienlist_file = config.get("files", "orienlist_file").unwrap();
        let quatfile = config.get("files", "quat_file").unwrap();
        let unit_cell_file = config.get("files", "unit_cell_file").unwrap();
        let nproc = config.get("reduce_data","nproc").unwrap().parse::<i32>()?;
        let num_row = config
            .get("make_detector", "num_row")
            .unwrap()
            .parse::<i32>()?;
        let num_col = config
            .get("make_detector", "num_col")
            .unwrap()
            .parse::<i32>()?;
        let detd = config
            .get("make_detector", "detd")
            .unwrap()
            .parse::<f64>()?;
        let wl = config.get("make_detector", "wl").unwrap().parse::<f64>()?;
        let px = config.get("make_detector", "px").unwrap().parse::<f64>()?;
        let cx = config.get("make_detector", "cx").unwrap().parse::<f64>()?;
        let cy = config.get("make_detector", "cy").unwrap().parse::<f64>()?;
        let Rstop = config
            .get("make_detector", "Rstop")
            .unwrap()
            .parse::<f64>()?;
        let num_data: usize = config
            .get("make_background", "num_raw_data")
            .unwrap()
            .parse::<usize>()?;
        let hot_pix_thres: i32 = config
            .get("make_background", "hot_pix_thres")
            .unwrap()
            .parse::<f64>()? as i32;
        let min_patch_sz = config
            .get("make_powder", "min_patch_sz")
            .unwrap()
            .parse::<u32>()?;
        let max_patch_sz = config
            .get("make_powder", "max_patch_sz")
            .unwrap()
            .parse::<u32>()?;
        let min_num_peak = config
            .get("make_powder", "min_num_peak")
            .unwrap()
            .parse::<u32>()?;
        let max_num_peak = config
            .get("make_powder", "max_num_peak")
            .unwrap()
            .parse::<u32>()?;

        let res_cutoff = config
            .get("orient_peak", "res_cutoff")
            .unwrap()
            .parse::<f64>()?;
        let VN = config.get("orient_peak", "VN").unwrap().parse::<i32>()?;
        let gw = config.get("orient_peak", "gw").unwrap().parse::<f64>()?;
        let tolmult = config
            .get("orient_peak", "tolerance_mult")
            .unwrap()
            .parse::<f64>()?;
        let delphi = config
            .get("orient_peak", "delta_phi")
            .unwrap()
            .parse::<f64>()?;
        let mut beam_vec: [f64; 3] = [
            config.get("make_detector", "sx").unwrap().parse::<f64>()?,
            config.get("make_detector", "sy").unwrap().parse::<f64>()?,
            config.get("make_detector", "sz").unwrap().parse::<f64>()?,
        ];

        //  read in outlier files
        let mut orient_files: Vec<String> = vec![Default::default(); num_data];
        let orient_file_buff = match File::open(&orienlist_file) {
            Err(why) => panic!("Couldn't open outlierlist_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in orient_file_buff.enumerate() {
            if let Ok(ip) = line {
                orient_files[i] = ip;
            }
        }

        //  read in peaklist files
        let mut peak_files: Vec<String> = vec![Default::default(); num_data];
        let peak_file_buff = match File::open(&peaklist_file) {
            Err(why) => panic!("Couldn't open peaklist_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in peak_file_buff.enumerate() {
            if let Ok(ip) = line {
                peak_files[i] = ip;
            }
        }
        let now = Instant::now();
        //read in quat file
        let quat_file = match File::open(quatfile) {
            Err(why) => panic!("Couldn't open quatfile : {}", why),
            Ok(file) => file,
        };
        let mut quat_file_buff: BufReader<File> = BufReader::new(quat_file);
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
        let elapsed = now.elapsed();
        println!("quat_file read took {:?}", elapsed);

        let norm = (beam_vec[0].powi(2) + beam_vec[1].powi(2) + beam_vec[2].powi(2)).sqrt();
        let D = (detd / px) * norm / (beam_vec[2]).abs();
        let rescale = (1.0 / wl) / D;

        beam_vec[0] *= D / norm;
        beam_vec[1] *= D / norm;
        beam_vec[2] *= D / norm;

        let qmax = (wl * D / res_cutoff).ceil() as i32;
        let qmin = (2.0 * D * sin(atan(Rstop / D) / 2.0)).ceil() as i32;
        let rmax2 = (qmax as f64 * rescale).powi(2);
        let rmin2 = (qmin as f64 * rescale).powi(2);
        let total_pix: i32 = num_row * num_col;

        /* pix[t] has unit 1/A */
        let mut pix: Vec<Pix> = vec![Pix(0., 0., 0.); total_pix as usize];
        for t in 0..total_pix {
            let x = (t / num_col) as f64 - cx;
            let y = (t % num_col) as f64 - cy;
            let mut sx = beam_vec[0] + x;
            let mut sy = beam_vec[1] + y;
            let mut sz = beam_vec[2];
            let norm = sqrt(sx * sx + sy * sy + sz * sz);
            sx *= D / norm;
            sy *= D / norm;
            sz *= D / norm;
            pix[t as usize].0 = (sx - beam_vec[0]) * rescale;
            pix[t as usize].1 = (sy - beam_vec[1]) * rescale;
            pix[t as usize].2 = (sz - beam_vec[2]) * rescale;
        }

        let mut unit_cell_buff = match File::open(&unit_cell_file) {
            Err(why) => panic!("Couldn't open unit cell file : {}", why),
            Ok(file) => BufReader::new(file),
        };
        let mut first_line = String::new();
        unit_cell_buff.read_line(&mut first_line)?;
        let mut ucline: Vec<&str> = first_line.trim().split_whitespace().collect();
        ucline.pop();
        let unit_cell: Vec<f64> = ucline
            .iter()
            .map(|word| word.parse::<f64>().unwrap())
            .collect();
        let cell: [f64; 6] = unit_cell
            .try_into()
            .unwrap_or_else(|_| panic!("Unit Cell wrong format: should be a b c al be gam sg"));

        let mut rot_mat_table: Vec<Vec<f64>> =
            vec![[Default::default(); 9].to_vec(); num_rot as usize];
        let mut rot: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        for r in 0..num_rot as usize {
            Config::make_rot(&quat, r as i32, &mut rot);
            let mut t = 0;
            for i in 0..3 {
                for j in 0..3 {
                    rot_mat_table[r][t] = rot[i][j];
                    t += 1;
                }
            }
        }
        let mut rcell: [f64; 6] = [0_f64; 6];
        let mut rlatt_vec: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        let mut cofactor: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        let mut mat: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        let mut proj: [[f64; 3]; 3] = [[0_f64; 3]; 3];
        let mut rvec: [f64; 3] = [0_f64; 3];
        let mut indices: [i32; 3] = [0; 3];

        let mut cell = cell.clone();
        let degtorad = atan(1.0) * 4.0 / 180.0;

        for i in 3..6 {
            cell[i] *= degtorad;
        }

        let sina = sin(cell[3]);
        let cosa = cos(cell[3]);
        let sinb = sin(cell[4]);
        let cosb = cos(cell[4]);
        let sing = sin(cell[5]);
        let cosg = cos(cell[5]);
        let cosas = (cosg * cosb - cosa) / (sinb * sing);
        let sinas = sqrt(1.0 - cosas * cosas);
        let cosbs = (cosa * cosg - cosb) / (sina * sing);
        let sinbs = sqrt(1.0 - cosbs * cosbs);
        let cosgs = (cosa * cosb - cosg) / (sina * sinb);
        let sings = sqrt(1.0 - cosgs * cosgs);

        let sumc = (cell[3] + cell[4] + cell[5]) * 0.5;
        let vol = 2.0
            * cell[0]
            * cell[1]
            * cell[2]
            * sqrt(sin(sumc - cell[3]) * sin(sumc - cell[4]) * sin(sumc - cell[5]) * sin(sumc));

        rcell[0] = cell[1] * cell[2] * sina / vol;
        rcell[1] = cell[2] * cell[0] * sinb / vol;
        rcell[2] = cell[0] * cell[1] * sing / vol;
        rcell[3] = atan2(sinas, cosas);
        rcell[4] = atan2(sinbs, cosbs);
        rcell[5] = atan2(sings, cosgs);

        rlatt_vec[0][0] = rcell[0];
        rlatt_vec[0][1] = 0.0;
        rlatt_vec[0][2] = 0.0;
        rlatt_vec[1][0] = rcell[1] * cosgs;
        rlatt_vec[1][1] = rcell[1] * sings;
        rlatt_vec[1][2] = 0.0;
        rlatt_vec[2][0] = rcell[2] * cosbs;
        rlatt_vec[2][1] = -rcell[2] * sinbs * cosa;
        rlatt_vec[2][2] = rcell[2] * sinbs * sina;

        for i in 0..3 {
            print!("reciprocal lattice vector {}:\n\t", i);
            for j in 0..3 {
                print!("{:.3e} ", rlatt_vec[i][j]);
            }
            print!("\n");
        }
        let mut tol = 0_f64;
        for j in 0..3 {
            let mut norm = 0_f64;
            for i in 0..3 {
                norm += pow(rlatt_vec[j][i], 2.);
            }
            if j == 0 {
                tol = norm;
            } else {
                if tol > norm {
                    tol = norm;
                }
            }
        }
        tol *= pow(gw / (VN as f64), 2.);

        /* calculate cofactors of rlatt_vec */
        for i in 0..3 {
            for j in 0..3 {
                if (i + j) % 2 == 0 {
                    cofactor[i][j] = 1_f64;
                } else {
                    cofactor[i][j] = -1_f64;
                }

                let mut i1 = 0;
                for i0 in 0..3 {
                    if i0 == i {
                        continue;
                    }
                    let mut j1 = 0;
                    for j0 in 0..3 {
                        if j0 == j {
                            continue;
                        }
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
        let mut det = 0.0;
        for i in 0..3 {
            det += rlatt_vec[0][i] * cofactor[0][i];
        }

        /* inverse of rlatt_vec */
        for i in 0..3 {
            for j in 0..3 {
                proj[i][j] = cofactor[j][i] / det;
            }
        }

        let d1 = delphi * 0.5 * degtorad;
        let tandp2 = tan(d1);
        /*
        if (myid == nproc-1)
          printf ("delphi %f\n",delphi);
        printf ("tandp2 %f\n",tandp2);
        */
        /* determine hmax, kmax, lmax */
        let mut norm = 0.;
        for i in 0..3 {
            norm += rlatt_vec[0][i].powi(2);
        }
        let mut hcur: i32 = floor(sqrt(rmax2 / norm)) as i32;

        norm = 0.;
        for i in 0..3 {
            norm += rlatt_vec[1][i].powi(2);
        }
        let mut kcur: i32 = floor(sqrt(rmax2 / norm)) as i32;

        norm = 0.;
        for i in 0..3 {
            norm += rlatt_vec[2][i].powi(2);
        }
        let mut lcur: i32 = floor(sqrt(rmax2 / norm)) as i32;

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
        let mut hkl_map: Vec<i32> = vec![Default::default(); (hlen * klen * llen) as usize];
        let mut hkl_unmap: Vec<i32> = vec![Default::default(); (hlen * klen * llen) as usize];

        let mut num_hkl: usize = 0;
        for i in -hmax..hmax + 1 {
            for j in -kmax..kmax + 1 {
                for k in -lmax..lmax + 1 {
                    indices[0] = i;
                    indices[1] = j;
                    indices[2] = k;
                    let mut norm = 0_f64;
                    for i0 in 0..3 {
                        rvec[i0] = 0_f64;
                        for j0 in 0..3 {
                            rvec[i0] += rlatt_vec[i0][j0] * indices[j0] as f64;
                        }
                        norm += rvec[i0] * rvec[i0];
                    }
                    let idx = (i + hmax) * klen * llen + (j + kmax) * llen + (k + lmax);
                    if norm > rmax2 || norm < rmin2 {
                        hkl_map[idx as usize] = -1;
                        continue;
                    }
                    hkl_map[idx as usize] = num_hkl as i32;
                    hkl_unmap[num_hkl as usize] = idx;
                    num_hkl += 1;
                }
            }
        }
        println!("hmax = {hmax}, kmax = {kmax}, lmax = {lmax}, num_hkl = {num_hkl}");
        let peak_rvec: Vec<Vec<r_vector>> = vec![vec![]; num_data];
        let num_peak: Vec<u32> = vec![0; num_data];
        Ok(Config {
            data_dir,
            radial_bg_dir,
            outlier_dir,
            data_info_file,
            file_type,
            home_dir,
            pixmap_file,
            rec2pixmap_file,
            rvec_file,
            cbflist_file,
            radiallist_file,
            outlierlist_file,
            peaklist_file,
            orienlist_file,
            unit_cell_file,
            num_row,
            num_col,
            detd,
            wl,
            px,
            cx,
            cy,
            num_data,
            hot_pix_thres,
            orient_files,
            peak_files,
            qmax,
            qmin,
            total_pix,
            min_patch_sz,
            max_patch_sz,
            min_num_peak,
            max_num_peak,
            res_cutoff,
            VN,
            gw,
            quat,
            tolmult,
            delphi,
            beam_vec,
            rmax2,
            rmin2,
            pix,
            cell,
            rot,
            rot_mat_table,
            num_rot,
            hmax,
            kmax,
            lmax,
            hlen,
            klen,
            llen,
            hkl_map,
            hkl_unmap,
            num_hkl,
            rlatt_vec,
            proj,
            tol,
            tandp2,
            peak_rvec,
            num_peak,
            nproc
        })
    }

    fn make_rot(quat: &Vec<Quat>, r: i32, rot: &mut [[f64; 3]; 3]) {
        let q0 = quat[r as usize].0;
        let q1 = quat[r as usize].1;
        let q2 = quat[r as usize].2;
        let q3 = quat[r as usize].3;

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
    }

    pub fn get_num_peak(&self) -> Result<(Vec<Vec<r_vector>>, Vec<u32>), Box<dyn Error>> {
        let dmin = 0;
        let dmax = self.num_data;
        let mut peak_rvec: Vec<Vec<r_vector>> = vec![vec![]; self.num_data];
        let mut num_peak: Vec<u32> = vec![0; self.num_data];
        for d in dmin..dmax {
            let mut weighted_vec: [f64; 3] = [0_f64; 3];
            num_peak[d] = 0;

            let fp_path: PathBuf = PathBuf::from(&self.peak_files[d]);
            let mut fp = match File::open(&fp_path) {
                Err(why) => {
                    return Err(format!("File {:?} can't be opened: {}", &fp_path, why).into())
                }
                Ok(file) => BufReader::new(file),
            };
            let mut buffer = String::new();
            let mut buf = fp.read_line(&mut buffer)?;
            while buf != 0 {
                let patch_sz: i32 = match buffer.trim().parse() {
                    Ok(int) => int,
                    Err(_) => {
                        return Err(format!("Expected digit, got {}: {:?}", buffer, &fp_path).into())
                    }
                };
                let mut val_sum = 0_f64;
                for i in 0..3 {
                    weighted_vec[i] = 0_f64;
                }
                buffer = String::new();
                fp.read_line(&mut buffer)?;
                let moo: Vec<&str> = buffer.trim().split_whitespace().collect();
                let mut spl_line: Vec<i32> =
                    moo.iter().map(|x| -> i32 { x.parse().unwrap() }).collect();
                // let mut spl_line:Vec<i32> = buffer.trim().split_whitespace().map(|x|-> i32 {x.parse().unwrap()}).collect();
                for t in 0..patch_sz {
                    let pix_val: i32 = spl_line.pop().unwrap();
                    let pix_id: i32 = spl_line.pop().unwrap();
                    weighted_vec[0] += self.pix[pix_id as usize].0 * pix_val as f64;
                    weighted_vec[1] += self.pix[pix_id as usize].1 * pix_val as f64;
                    weighted_vec[2] += self.pix[pix_id as usize].2 * pix_val as f64;
                    val_sum += pix_val as f64;
                }
                if patch_sz < self.min_patch_sz.try_into().unwrap()
                    || patch_sz > self.max_patch_sz.try_into().unwrap()
                {
                    buffer = String::new();
                    buf = fp.read_line(&mut buffer)?;
                    continue;
                }
                let mut qval2 = 0_f64;
                let inv_val_sum = 1.0 / val_sum;
                for i in 0..3 {
                    weighted_vec[i] *= inv_val_sum;
                    qval2 += weighted_vec[i] * weighted_vec[i];
                }
                if qval2 > self.rmax2 || qval2 < self.rmin2 {
                    buffer = String::new();
                    buf = fp.read_line(&mut buffer)?;
                    continue;
                }
                peak_rvec[d as usize].push(r_vector {
                    vec: (weighted_vec),
                });
                num_peak[d] += 1;
                buffer = String::new();
                buf = fp.read_line(&mut buffer)?;
            }

            /* not enough peaks in the frame */
            if num_peak[d] < self.min_num_peak {
                continue;
            }
            if num_peak[d] > self.max_num_peak {
                num_peak[d] = 0;
                continue;
            }
        }
        Ok((peak_rvec, num_peak))
    }
}
#[derive(Clone, Debug)]

pub struct Results {
    pub d: i32,
    pub prob_orien: Vec<i32>,
    pub num_prob_orien: u32,
    pub max_num_prob_orien: u32,
    pub found_hkl: Vec<u32>,
    pub n_found_hkl: Vec<u32>,
    pub rmatx: Vec<f64>,
    pub score: Vec<f64>,
}

impl Results {
    pub fn new(d:i32) -> Result<Results, Box<dyn Error>> {
        let d: i32 = d;
        let prob_orien: Vec<i32> = vec![];
        let num_prob_orien: u32 = 0;
        let max_num_prob_orien: u32 = 0;

        let found_hkl: Vec<u32> = vec![];
        let n_found_hkl: Vec<u32> = vec![];
        let rmatx: Vec<f64> = vec![];
        let score: Vec<f64> = vec![];
        Ok(Results {
            d,
            prob_orien,
            num_prob_orien,
            max_num_prob_orien,
            found_hkl,
            n_found_hkl,
            rmatx,
            score,
        })
    }
}
