use cryiorust::cbf::Cbf;
use cryiorust::frame::Frame;
use hdf5::File as H5File;
use ndarray::Array1;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

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

pub fn load_pixmap_dat(pixmap_path: &Path, total_pix: usize) -> Result<Vec<i32>, Box<dyn Error>> {
    let mut pix_map: Vec<i32> = vec![Default::default(); total_pix];
    let pix_map_buff = match File::open(&pixmap_path) {
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
    total_pix: usize,
) -> Result<Vec<i32>, Box<dyn Error>> {
    // let mut pix_map: Vec<i32> = vec![Default::default();self.total_pix];
    let pix_map_buff = match H5File::open(&pixmap_path) {
        Err(why) => return Err(format!("Couldn't open pixmap_file: {}", why).into()),
        Ok(file) => file,
    };
    let data = pix_map_buff.group("data")?;
    let pix_map_data = data.dataset(dataset_name)?.read_1d()?;
    let pix_map: Vec<i32> = pix_map_data.to_vec();
    assert_eq!(total_pix, pix_map.len(), "pixmap does not equal total pix");
    println!("Pixmap {:?} read in", pixmap_path);
    Ok(pix_map)
}
