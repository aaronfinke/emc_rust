use configparser::ini::Ini;
use hdf5::File as H5File;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::io::{Seek, SeekFrom};
use std::path::{Path, PathBuf};

pub struct Config {
    pub home_dir: String,
    pub raw_data_dir: String,
    pub required_string: String,
    pub reduced_data_dir: String,
    pub data_dir: String,
    pub prob_dir: String,
    pub radial_bg_dir: String,
    pub outlier_dir: String,
    pub start_prob_dir: String,
    pub mask_file: String,
    pub pixmap_file: String,
    pub rec2pix_file: String,
    pub rec_vectors_file: String,
    pub unit_cell_file: String,
    pub quat_file: String,
    pub start_phi_file: String,
    pub start_intens_file: String,
    pub cbflist_file: String,
    pub radiallist_file: String,
    pub orienlist_file: String,
    pub outlierlist_file: String,
    pub peaklist_file: String,
    pub data_info_file: String,
    pub reduced_data_id_file: String,
    pub reduced_cbflist_file: String,
    pub reduced_radiallist_file: String,
    pub reduced_peaklist_file: String,
    pub r2peak_file: String,
    pub peak2r_file: String,
    pub mpi_bgfile: String,
    pub mpi_datafile: String,
    pub prob_orien_file: String,
    pub num_prob_orien_file: String,
    pub fine_quat_file: String,
    pub quat_table_file: String,
    pub local_r2peak_file: String,
    pub local_peak2r_file: String,
    pub file_type: String,
    pub num_row: i32,
    pub num_col: i32,
    pub detd: f64,
    pub wl: f64,
    pub px: f64,
    pub cx: f64,
    pub cy: f64,
    pub bh_width: i32,
    pub Rstop: f64,
    pub res_max: f64,
    pub num_raw_data: usize,
    pub hot_pix_thres: i32,
    pub qlen: i32,
    pub min_patch_sz: i32,
    pub max_patch_sz: i32,
    pub min_num_peak: i32,
    pub max_num_peak: i32,
    pub res_cutoff: f64,
    pub VN: i32,
    pub gw: f64,
    pub tolmult: f64,
    pub delphi: f64,
    pub nproc: i32,
    pub iter_data_block: i32,
    pub high_res_cutoff: f64,
    pub beam_vec: [f64; 3],
}

impl Config {
    pub fn new(args: &str) -> Result<Config, Box<dyn Error>> {
        let mut config = Ini::new();
        let config_file = args;
        config.load(config_file)?;
        //files and directories
        let home_dir = match config.get("files", "home_dir") {
            Some(name) => name,
            None => return Err("can't find home_dir".into()),
        };
        let raw_data_dir = match config.get("files", "raw_data_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing raw_data_dir").into()),
        };
        let required_string = match config.get("files", "required_string") {
            Some(name) => name,
            None => return Err(format!("Error parsing required_string").into()),
        };
        let reduced_data_dir = match config.get("files", "reduced_data_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing reduced_data_dir").into()),
        };
        let data_dir = match config.get("files", "data_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing data_dir").into()),
        };
        let prob_dir = match config.get("files", "prob_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing prob_dir").into()),
        };
        let radial_bg_dir = match config.get("files", "radial_bg_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing radial_bg_dir").into()),
        };
        let outlier_dir = match config.get("files", "outlier_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing outlier_dir").into()),
        };
        let start_prob_dir = match config.get("files", "start_prob_dir") {
            Some(name) => name,
            None => return Err(format!("Error parsing start_prob_dir").into()),
        };
        let mask_file = match config.get("files", "mask_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing mask_file").into()),
        };
        let pixmap_file = match config.get("files", "pixmap_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing pixmap_file").into()),
        };
        let rec2pix_file = match config.get("files", "rec2pix_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing rec2pix_file").into()),
        };
        let rec_vectors_file = match config.get("files", "rec_vectors_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing rec_vectors_file").into()),
        };
        let unit_cell_file = match config.get("files", "unit_cell_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing unit_cell_file").into()),
        };
        let quat_file = match config.get("files", "quat_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing quat_file").into()),
        };
        let start_phi_file = match config.get("files", "start_phi_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing start_phi_file").into()),
        };
        let start_intens_file = match config.get("files", "start_intens_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing start_intens_file").into()),
        };
        let cbflist_file = match config.get("files", "cbflist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing cbflist_file").into()),
        };
        let radiallist_file = match config.get("files", "radiallist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing radiallist_file").into()),
        };
        let orienlist_file = match config.get("files", "orienlist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing orienlist_file").into()),
        };
        let outlierlist_file = match config.get("files", "outlierlist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing outlierlist_file").into()),
        };
        let peaklist_file = match config.get("files", "peaklist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing peaklist_file").into()),
        };
        let data_info_file = match config.get("files", "data_info_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing data_info_file").into()),
        };
        let reduced_data_id_file = match config.get("files", "reduced_data_id_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing reduced_data_id_file").into()),
        };
        let reduced_cbflist_file = match config.get("files", "reduced_cbflist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing reduced_cbflist_file").into()),
        };
        let reduced_radiallist_file = match config.get("files", "reduced_radiallist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing reduced_radiallist_file").into()),
        };
        let reduced_peaklist_file = match config.get("files", "reduced_peaklist_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing reduced_peaklist_file").into()),
        };
        let r2peak_file = match config.get("files", "r2peak_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing r2peak_file").into()),
        };
        let peak2r_file = match config.get("files", "peak2r_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing peak2r_file").into()),
        };
        let mpi_bgfile = match config.get("files", "mpi_bgfile") {
            Some(name) => name,
            None => return Err(format!("Error parsing mpi_bgfile").into()),
        };
        let mpi_datafile = match config.get("files", "mpi_datafile") {
            Some(name) => name,
            None => return Err(format!("Error parsing mpi_datafile").into()),
        };
        let prob_orien_file = match config.get("files", "prob_orien_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing prob_orien_file").into()),
        };
        let num_prob_orien_file = match config.get("files", "num_prob_orien_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing num_prob_orien_file").into()),
        };
        let fine_quat_file = match config.get("files", "fine_quat_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing fine_quat_file").into()),
        };
        let quat_table_file = match config.get("files", "quat_table_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing quat_table_file").into()),
        };
        let local_r2peak_file = match config.get("files", "local_r2peak_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing local_r2peak_file").into()),
        };
        let local_peak2r_file = match config.get("files", "local_peak2r_file") {
            Some(name) => name,
            None => return Err(format!("Error parsing local_peak2r_file").into()),
        };
        // make_detector
        let file_type = match config.get("make_detector", "file_type") {
            Some(name) => name,
            None => return Err(format!("Error parsing file_type").into()),
        };

        let num_row: i32 = match config.get("make_detector", "num_row") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing num_row: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let num_col = match config.get("make_detector", "num_col") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing num_col: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let detd = match config.get("make_detector", "detd") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing detd: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let wl = match config.get("make_detector", "wl") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing wl: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let px = match config.get("make_detector", "px") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing px: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let cx = match config.get("make_detector", "cx") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing cx: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let cy = match config.get("make_detector", "cy") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing cy: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let bh_width = match config.get("make_detector", "bh_width") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing bh_width: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let Rstop = match config.get("make_detector", "Rstop") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing Rstop: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let res_max = match config.get("make_detector", "res_max") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing res_max: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };

        let sx = match config.get("make_detector", "sx") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing sx: {}", why).into()),
            },
            None => return Err(format!("Error parsing sx").into()),
        };
        let sy = match config.get("make_detector", "sy") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing sy: {}", why).into()),
            },
            None => return Err(format!("Error parsing sy").into()),
        };
        let sz = match config.get("make_detector", "sz") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing sz: {}", why).into()),
            },
            None => return Err(format!("Error parsing sz").into()),
        };

        let beam_vec: [f64; 3] = [sx, sy, sz];
        // make_background
        let num_raw_data: usize = match config.get("make_background", "num_raw_data") {
            Some(num) => match num.parse::<usize>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing num_raw_data: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_raw_data").into()),
        };

        let hot_pix_thres: i32 = match config.get("make_background", "hot_pix_thres") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing hot_pix_thres: {}", why).into()),
            },
            None => return Err(format!("Error parsing hot_pix_thres").into()),
        };
        let qlen: i32 = match config.get("make_background", "qlen") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing qlen: {}", why).into()),
            },
            None => return Err(format!("Error parsing qlen").into()),
        };
        // make_powder
        let min_patch_sz = match config.get("make_powder", "min_patch_sz") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing min_patch_sz: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let max_patch_sz = match config.get("make_powder", "max_patch_sz") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing max_patch_sz: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let min_num_peak = match config.get("make_powder", "min_num_peak") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing min_num_peak: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let max_num_peak = match config.get("make_powder", "max_num_peak") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing max_num_peak: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        // orient_peak
        let res_cutoff = match config.get("orient_peak", "res_cutoff") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing res_cutoff: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let VN = match config.get("orient_peak", "VN") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing VN: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let gw = match config.get("orient_peak", "gw") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing gw: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let tolmult = match config.get("orient_peak", "tolerance_mult") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing tolmult: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        let delphi = match config.get("orient_peak", "delta_phi") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing delphi: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        //reduce_data
        let nproc = match config.get("reduce_data", "nproc") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing nproc: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        //low_res_emc
        let iter_data_block = match config.get("low_res_emc", "iter_data_block") {
            Some(num) => match num.parse::<i32>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing iter_data_block: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };
        // high_res_emc
        let high_res_cutoff = match config.get("high_res_emc", "high_res_cutoff") {
            Some(num) => match num.parse::<f64>() {
                Ok(num) => num,
                Err(why) => return Err(format!("Error parsing high_res_cutoff: {}", why).into()),
            },
            None => return Err(format!("Error parsing num_row").into()),
        };

        //  make quat

        Ok(Config {
            home_dir,
            raw_data_dir,
            required_string,
            reduced_data_dir,
            data_dir,
            prob_dir,
            radial_bg_dir,
            outlier_dir,
            start_prob_dir,
            mask_file,
            pixmap_file,
            rec2pix_file,
            rec_vectors_file,
            unit_cell_file,
            quat_file,
            start_phi_file,
            start_intens_file,
            cbflist_file,
            radiallist_file,
            orienlist_file,
            outlierlist_file,
            peaklist_file,
            data_info_file,
            reduced_data_id_file,
            reduced_cbflist_file,
            reduced_radiallist_file,
            reduced_peaklist_file,
            r2peak_file,
            peak2r_file,
            mpi_bgfile,
            mpi_datafile,
            prob_orien_file,
            num_prob_orien_file,
            fine_quat_file,
            quat_table_file,
            local_r2peak_file,
            local_peak2r_file,
            file_type,
            num_row,
            num_col,
            detd,
            wl,
            px,
            cx,
            cy,
            bh_width,
            Rstop,
            res_max,
            num_raw_data,
            hot_pix_thres,
            qlen,
            min_patch_sz,
            max_patch_sz,
            min_num_peak,
            max_num_peak,
            res_cutoff,
            VN,
            gw,
            tolmult,
            delphi,
            nproc,
            iter_data_block,
            high_res_cutoff,
            beam_vec,
        })
    }
    pub fn load_pixmap_dat(
        pixmap_path: &Path,
        total_pix: usize,
    ) -> Result<Vec<i32>, Box<dyn Error>> {
        let mut pix_map: Vec<i32> = vec![Default::default(); total_pix];
        let pix_map_buff = match std::fs::File::open(&pixmap_path) {
            Err(why) => return Err(format!("Couldn't open pixmap_file: {}", why).into()),
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

    pub fn load_pixmap_h5(
        pixmap_path: &Path,
        dataset_name: &str,
    ) -> Result<Vec<i32>, Box<dyn Error>> {
        // let mut pix_map: Vec<i32> = vec![Default::default();self.total_pix];
        let pix_map_buff = match H5File::open(&pixmap_path) {
            Err(why) => return Err(format!("Couldn't open pixmap_file: {}", why).into()),
            Ok(file) => file,
        };
        let data = pix_map_buff.group("data")?;
        let pix_map_data = data.dataset(dataset_name)?.read_1d()?;
        let pix_map: Vec<i32> = pix_map_data.to_vec();
        println!("Pixmap {:?} read in", pixmap_path);
        Ok(pix_map)
    }
}
