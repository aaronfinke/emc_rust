use configparser::ini::Ini;
use cryiorust::cbf::Cbf;
use cryiorust::frame::Frame;
use hdf5::File as H5File;
use ndarray::Array1;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ffi::OsStr;
use std::path::Path;

#[derive(Copy, Clone, Debug)]
pub struct FrameInfo {
    pub row_size: i32,
    pub col_size: i32,
    pub exposure: f64,
    pub wl: f64,
    pub detd: f64,
    pub c_row: f64,
    pub c_col: f64,
}
impl FrameInfo {
    pub fn get_frame_info_exact(&mut self, cbf: &dyn Frame) {
        self.col_size = cbf
            .get_header_i64("X-Binary-Size-Fastest-Dimension: ")
            .unwrap() as i32;
        self.row_size = cbf
            .get_header_i64("X-Binary-Size-Second-Dimension: ")
            .unwrap() as i32;
        self.exposure = cbf.get_header_float("# Exposure_time ");
        self.wl = cbf.get_header_float("# Wavelength ");
        self.detd = cbf.get_header_float("# Detector_distance ");
        if let Some(cryiorust::frame::HeaderEntry::Pixels(x)) = cbf.header().get("# Beam_xy ") {
            self.c_col = x[0];
            self.c_row = x[1];
        }
    }
    pub fn get_frame_info(&mut self, cbf: &Cbf) -> Result<(), Box<dyn Error>> {
        self.col_size = cbf
            .get_header_i64("X-Binary-Size-Fastest-Dimension: ")
            .unwrap() as i32;
        self.row_size = cbf
            .get_header_i64("X-Binary-Size-Fastest-Dimension: ")
            .unwrap() as i32;
        self.exposure = cbf.get_header_float(cryiorust::cbf::KEY_EXPOSURE_TIME);
        self.wl = cbf.get_header_float(cryiorust::cbf::KEY_WAVELENGTH);
        self.detd = cbf.get_header_float(cryiorust::cbf::KEY_DETECTOR_DISTANCE);
        if let Some(cryiorust::frame::HeaderEntry::Pixels(x)) =
            cbf.header().get(cryiorust::cbf::KEY_BEAM_XY)
        {
            self.c_col = x[0];
            self.c_row = x[1];
        };
        Ok(())
    }
}

#[derive(Copy, Clone)]
pub struct Mpix(pub f64, pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub struct PeakLocation(pub i32, pub i32);

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
    pub num_row: i32,
    pub num_col: i32,
    pub detd: f64,
    pub wl: f64,
    pub px: f64,
    pub cx: f64,
    pub cy: f64,
    pub qlen: i32,
    pub num_data: usize,
    pub hot_pix_thres: i32,
    pub qexc_d: Vec<Vec<f64>>,
    pub n_qexc: usize,
    pub qexc: Vec<Vec<i32>>,
    pub cbf_files: Vec<String>,
    pub outlier_files: Vec<String>,
    pub radial_files: Vec<String>,
    pub peak_files: Vec<String>,
    pub qx: Array1<f64>,
    pub qy: Array1<f64>,
    pub qz: Array1<f64>,
    pub scale_factor: Array1<f64>,
    pub qmax: i32,
    pub qmin: i32,
    pub num_pix: i32,
    pub total_pix: usize,
    pub pix_map: Vec<i32>,
    pub pix: Vec<Mpix>,
    pub qid_map: Vec<i32>,
    pub rec2pix_map: Vec<i32>,
    pub dq: f64,
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
        let pixmap_path = Path::new(&pixmap_file);
        let rec2pixmap_file = config.get("files", "rec2pix_file").unwrap();
        let rec2pixmap_path = Path::new(&rec2pixmap_file);
        let rvec_file = config.get("files", "rec_vectors_file").unwrap();
        let cbflist_file = config.get("files", "cbflist_file").unwrap();
        let radiallist_file = config.get("files", "radiallist_file").unwrap();
        let outlierlist_file = config.get("files", "outlierlist_file").unwrap();
        let peaklist_file = config.get("files", "peaklist_file").unwrap();
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
        let qlen = config
            .get("make_background", "qlen")
            .unwrap()
            .parse::<i32>()?;
        let num_data: usize = config
            .get("make_background", "num_raw_data")
            .unwrap()
            .parse::<usize>()?;
        let hot_pix_thres: i32 = config
            .get("make_background", "hot_pix_thres")
            .unwrap()
            .parse::<f64>()? as i32;

        let exclude_res_inner_str = match config.get("make_background", "exclude_res_inner") {
            Some(_) => config.get("make_background", "exclude_res_inner").unwrap(),
            None => "".to_string(),
        };
        let exclude_res_outer_str = match config.get("make_background", "exclude_res_outer") {
            Some(_) => config.get("make_background", "exclude_res_outer").unwrap(),
            None => "".to_string(),
        };
        //resolution exclusion zones
        let exclude_res_inner: Vec<&str> = exclude_res_inner_str.split_whitespace().collect();
        let exclude_res_outer: Vec<&str> = exclude_res_outer_str.split_whitespace().collect();
        if exclude_res_inner.len() != exclude_res_outer.len() {
            return Err("exclude_res_inner and exclude_res_outer should have same values".into());
        }
        let n_qexc = exclude_res_inner.len();
        let mut qexc_d: Vec<Vec<f64>> = vec![vec![Default::default(); 2]; n_qexc];
        let mut qexc: Vec<Vec<i32>> = vec![vec![Default::default(); 2]; n_qexc];
        if n_qexc > 0 {
            for i in 0..n_qexc {
                let d1: f64 = exclude_res_inner[i].parse().unwrap();
                let d2: f64 = exclude_res_outer[i].parse().unwrap();
                qexc_d[i][0] = wl * (detd / px) / d1;
                qexc_d[i][1] = wl * (detd / px) / d2;
            }
        }
        //  read in CBF files
        let mut cbf_files: Vec<String> = vec![Default::default(); num_data];
        let cbf_file_buff = match File::open(&cbflist_file) {
            Err(why) => panic!("Couldn't open cbflist_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in cbf_file_buff.enumerate() {
            if let Ok(ip) = line {
                cbf_files[i] = ip;
            }
        }
        //  read in outlier files
        let mut outlier_files: Vec<String> = vec![Default::default(); num_data];
        let outlier_file_buff = match File::open(&outlierlist_file) {
            Err(why) => panic!("Couldn't open outlierlist_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in outlier_file_buff.enumerate() {
            if let Ok(ip) = line {
                outlier_files[i] = ip;
            }
        }
        //  read in radial files
        let mut radial_files: Vec<String> = vec![Default::default(); num_data];
        let radial_file_buff = match File::open(&radiallist_file) {
            Err(why) => panic!("Couldn't open radiallist_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in radial_file_buff.enumerate() {
            if let Ok(ip) = line {
                radial_files[i] = ip;
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
        // read in rvec file
        let rvec_file_buff = H5File::open(&rvec_file)?.group("data").unwrap();
        let qx: Array1<f64> = rvec_file_buff.dataset("qx")?.read_1d().unwrap();
        let qy: Array1<f64> = rvec_file_buff.dataset("qy")?.read_1d().unwrap();
        let qz: Array1<f64> = rvec_file_buff.dataset("qz")?.read_1d().unwrap();
        let scale_factor: Array1<f64> = rvec_file_buff.dataset("scale_factor")?.read_1d().unwrap();
        let qmax: i32 = rvec_file_buff.attr("qmax")?.read_scalar().unwrap();
        let qmin: i32 = rvec_file_buff.attr("qmin")?.read_scalar().unwrap();
        let num_pix: i32 = rvec_file_buff.attr("Npix")?.read_scalar().unwrap();
        let total_pix: usize = (num_row * num_col).try_into().unwrap();

        //read in pixmap file
        
        let pix_map: Vec<i32> = match pixmap_path.extension().and_then(OsStr::to_str)
                    {
                        Some("dat") => Config::load_pixmap_dat(pixmap_path, total_pix)?,
                        Some("h5") => Config::load_pixmap_h5(pixmap_path, "pix_map", total_pix)?,
                        _ => return Err("could not determine pixmap type\nshould be .h5 or .dat (if text file)".into())
                    };
        let qval_max: f64 = qmax.try_into().unwrap();
        let dq: f64 = qval_max / (qlen as f64);
        if n_qexc > 0 {
            for i in 0..n_qexc {
                let mut idx: i32 = ((qexc_d[i][0] / dq - 0.5).round()) as i32;
                if idx > (qlen - 1) {
                    idx = qlen - 1;
                }
                qexc[i][0] = idx;
                idx = ((qexc_d[i][1] / dq - 0.5).round()) as i32;
                if idx < 0 {
                    idx = 0
                };
                if idx > (qlen - 1) {
                    idx = qlen - 1
                };
                qexc[i][1] = idx;
                println!("Excluding q bins {} to {}", qexc[i][0], qexc[i][1]);
            }
        }
        let mut pix: Vec<Mpix> = vec![Mpix(0.0, 0.0, 0.0, 0.0); num_pix as usize];
        let mut qid_map: Vec<i32> = vec![0; num_pix as usize];
        for t in 0..num_pix as usize {
            pix[t].0 = qx[t];
            pix[t].1 = qy[t];
            pix[t].2 = qz[t];
            pix[t].3 = scale_factor[t];
            let qval = (pix[t].0 * pix[t].0 + pix[t].1 * pix[t].1 + pix[t].2 * pix[t].2).sqrt();
            let idx: i32 = ((qval / dq - 0.5).round()) as i32;
            if idx < 0 || idx > qlen - 1 {
                qid_map[t] = -1
            } else {
                qid_map[t] = idx
            };
        }

        let rec2pix_map: Vec<i32> = match rec2pixmap_path.extension().and_then(OsStr::to_str)
        {
            Some("dat") => Config::load_pixmap_dat(rec2pixmap_path, num_pix as usize)?,
            Some("h5") => Config::load_pixmap_h5(rec2pixmap_path, "rec2pix_map", num_pix as usize)?,
            _ => return Err("could not determine pixmap type\nshould be .h5 or .dat (if text file)".into())
        };

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
            num_row,
            num_col,
            detd,
            wl,
            px,
            cx,
            cy,
            qlen,
            num_data,
            hot_pix_thres,
            qexc_d,
            n_qexc,
            qexc,
            cbf_files,
            outlier_files,
            radial_files,
            peak_files,
            qx,
            qy,
            qz,
            scale_factor,
            qmax,
            qmin,
            num_pix,
            total_pix,
            pix_map,
            pix,
            qid_map,
            rec2pix_map,
            dq,
        })
    }
    fn load_pixmap_dat(pixmap_path: &Path, total_pix: usize) -> Result<Vec<i32>, Box<dyn Error>> {
        let mut pix_map: Vec<i32> = vec![Default::default();total_pix];
        let pix_map_buff = match File::open(&pixmap_path) {
        Err(why) => return Err(format!("Couldn't open pixmap_file: {}",why).into()),
        Ok(file) => BufReader::new(file).lines(),
        };
        for (i, line) in pix_map_buff.enumerate() {
        if let Ok(ip) = line {
            pix_map[i] = ip.parse::<i32>().unwrap();
        }
    }
    println!("Pixmap {:?} read in", pixmap_path);
    Ok(pix_map)
    }

    fn load_pixmap_h5(pixmap_path: &Path, dataset_name:&str, total_pix: usize) -> Result<Vec<i32>,Box<dyn Error>> {
        // let mut pix_map: Vec<i32> = vec![Default::default();self.total_pix];
        let pix_map_buff = match H5File::open(&pixmap_path) {
            Err(why) => return Err(format!("Couldn't open pixmap_file: {}",why).into()),
            Ok(file) => file,
            };
        let data = pix_map_buff.group("data")?;
        let pix_map_data = data.dataset(dataset_name)?.read_1d()?;
        let pix_map: Vec<i32> = pix_map_data.to_vec();
        assert_eq!(total_pix, pix_map.len(), "pixmap does not equal total pix");
        println!("Pixmap {:?} read in", pixmap_path);
        Ok(pix_map)
    }
}
