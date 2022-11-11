use clap::Parser;
use hdf5::File as H5File;
use std::error::Error;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy)]

pub struct Quat(pub f64, pub f64, pub f64, pub f64, pub f64);
#[derive(Clone, Copy)]
pub struct Pix(pub f64, pub f64, pub f64);

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    /// make Ematrix for low-resolution reconstruction
    #[clap(long, parse(from_flag))]
    pub low: bool,

    /// Path to config_file
    #[clap(value_parser, value_name = "FILE")]
    pub config: PathBuf,
}

pub fn load_pixmap_dat(pixmap_path: &Path, total_pix: usize) -> Result<Vec<i32>, Box<dyn Error>> {
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
