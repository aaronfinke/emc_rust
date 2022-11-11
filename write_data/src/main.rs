use crate::config::configure::Config;
use cryiorust::cbf::Cbf;
use cryiorust::frame::Frame;
use hdf5::File as H5File;
use ndarray::Array1;
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
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

use write_data::{load_pixmap_dat, load_pixmap_h5, FrameInfo};

pub mod config;
/*
read raw data, rewrite data in short integers in HDF5 format

usage:
write_data [config file]

needs:
config_file,
prob_orien_file, reduced_cbflist_file, reduced_peaklist_file,
reduced_data_id file
pixmap_file, rec_vectors_file, rec2pix_file (specified in config file)
num_prob_orien.dat (in data_dir),
raw data files (in raw_data_dir)

makes:
mpi_datafile  (location specified in config file)

*/
const SHORT_MAX: i32 = 32767;

fn main() -> Result<(), Box<dyn Error>> {
    let t1 = Instant::now();

    let args: Vec<String> = env::args().collect();
    match args.len() {
        2 => {}
        _ => return Err("execution: ./make_background (/path/to/)config.ini".into()),
    };

    /* Setup section */
    let config = match Config::new(&args) {
        Ok(config) => config,
        Err(why) => return Err(format!("Error in config.ini: {:?}", why).into()),
    };

    //read in reduced_data_id file
    let mut reduced_data_id_buff = match File::open(&config.reduced_data_id_file) {
        Err(why) => panic!("Couldn't open reduced_data_id_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    //num_data is the first line in reduced_data_id_file
    let num_data: usize = reduced_data_id_buff
        .next()
        .unwrap()?
        .parse::<usize>()
        .expect(&format!(
            "Can't extract num_data from {:?}",
            config.reduced_data_id_file
        ));
    let mut data_id: Vec<i32> = Vec::with_capacity(num_data);
    for line in reduced_data_id_buff {
        if let Ok(ip) = line {
            data_id.push(ip.parse::<i32>().expect(&format!(
                "Can't extract data_id from {:?}",
                config.reduced_data_id_file
            )));
        }
    }
    assert_eq!(
        num_data,
        data_id.len(),
        "num_data and data_ids do not match!"
    );

    let data_info: Vec<FrameInfo> = Vec::with_capacity(num_data);

    //read in reduced_cbf file
    let mut cbf_files: Vec<String> = Vec::with_capacity(num_data);
    let cbf_file_buff = match File::open(&config.reduced_cbflist_file) {
        Err(why) => panic!("Couldn't open reduced_data_id_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for line in cbf_file_buff {
        if let Ok(ip) = line {
            cbf_files.push(ip);
        }
    }

    //read in reduced_peaklist_file file
    let mut peakfiles: Vec<String> = Vec::with_capacity(num_data);
    let peaklist_file_buff = match File::open(&config.reduced_peaklist_file) {
        Err(why) => panic!("Couldn't open reduced_data_id_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for line in peaklist_file_buff {
        if let Ok(ip) = line {
            peakfiles.push(ip);
        }
    }

    let (qmax, qmin, num_pix) = {
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
        let qx: Array1<f64> = rvec_file_buff.dataset("qx")?.read_1d().unwrap();
        let qy: Array1<f64> = rvec_file_buff.dataset("qy")?.read_1d().unwrap();
        let qz: Array1<f64> = rvec_file_buff.dataset("qz")?.read_1d().unwrap();
    
        (qmax, qmin, num_pix)
    };
    let rec2pixmap_path = Path::new(&config.rec2pix_file);
    let rec2pix_map: Vec<i32> = match rec2pixmap_path.extension().and_then(OsStr::to_str) {
        Some("dat") => load_pixmap_dat(rec2pixmap_path, num_pix as usize)?,
        Some("h5") => load_pixmap_h5(rec2pixmap_path, "rec2pix_map", num_pix as usize)?,
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
        Some("dat") => load_pixmap_dat(pixmap_path, total_pix)?,
        Some("h5") => load_pixmap_h5(pixmap_path, "pix_map", total_pix)?,
        _ => {
            return Err(format!(
                "could not determine pixmap type\nshould be .h5 or .dat (if text file)"
            )
            .into())
        }
    };
    assert_eq!(total_pix, pix_map.len());
    /* end setup section */

    /* check_reduce section */

    let count_d2r: Vec<i32> = {
        let mut count_d2r: Vec<i32> = Vec::with_capacity(num_data);
        let mut prob_orien_buff: File = match File::open(&config.prob_orien_file) {
            Err(why) => panic!("Couldn't open prob_orien_file : {}", why),
            Ok(file) => file,
        };
        for d in 0..num_data {
            let mut buff = vec![0_u8; 4];
            prob_orien_buff
                .read_exact(&mut buff)
                .expect("buff not read in");
            count_d2r.push(i32::from_le_bytes(buff.try_into().unwrap()));
            prob_orien_buff
                .seek(SeekFrom::Current((count_d2r[d] * 4) as i64))
                .expect("seek failed");
        }
        count_d2r
    };
    let num_orien: Vec<i32> = {
        let mut num_orien: Vec<i32> = Vec::with_capacity(config.num_raw_data);
        let num_prob_orien_buff = match File::open(&config.num_prob_orien_file) {
            Err(why) => panic!("Couldn't open num_prob_orien_file : {}", why),
            Ok(file) => BufReader::new(file).lines(),
        };
        for line in num_prob_orien_buff {
            if let Ok(ip) = line {
                num_orien.push(ip.parse::<i32>().expect(&format!(
                    "Can't extract num_orien from {:?}",
                    config.num_prob_orien_file
                )))
            };
        }
        num_orien
    };
    for d in 0..num_data {
        let idx: usize = data_id[d] as usize;
        if count_d2r[d] != num_orien[idx] {
            return Err("reduce_data fails!!".into());
        }
    }

    /* end check_reduce section */

    /* write_data section */

    let total_pix: usize = (config.num_row * config.num_col) as usize;
    let det: Vec<i32> = Vec::with_capacity(total_pix);
    // register_bitshuffle_plugin();

    let out_fp = match H5File::create(&config.mpi_datafile) {
        Err(why) => panic!("couldn't open {}: {}", &config.mpi_datafile, why),
        Ok(file) => file,
    };
    let data = out_fp.create_group("data").expect("Can't create group");
    println!("num_data = {num_data}");
    for d in 0..num_data {
        // let t2 = Instant::now();
        let mut data_frame: Vec<i16> = Vec::with_capacity(num_pix as usize);
        let det: Vec<i32> = match &config.file_type as &str  {
            "cbf" => load_cbf(&cbf_files[d]).expect("Error loading cbf_file"),
            "h5" => return Err("no support for h5 yet".into()),
            _ => return Err("no support for h5 yet".into()),
        };
        if config.file_type == String::from("cbf") {
            let det = load_cbf(&cbf_files[d]).expect("Error loading cbf_file");
            } else if config.file_type == String::from("h5") {
            return Err("no support for h5 yet".into());
        } else {
            return Err("check file type, either cbf or hdf5".into());
        }
        let idx = data_id[d];
        for t in 0..num_pix as usize {
            let pid = rec2pix_map[t];
            let photon_count = match det[pid as usize] {
                x if x < 0 => -1,
                x if x > config.hot_pix_thres => -1,
                x if x > SHORT_MAX => -1,
                _ => det[pid as usize],
            };
            data_frame.push(photon_count as i16);
        }
        let fp_path: PathBuf = PathBuf::from(&peakfiles[d]);
        let mut fp = match File::open(&fp_path) {
            Err(why) => return Err(format!("File {:?} can't be opened: {}", &fp_path, why).into()),
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
            if patch_sz > config.max_patch_sz {
                buffer = String::new();
                fp.read_line(&mut buffer)?;
                let moo: Vec<&str> = buffer.trim().split_whitespace().collect();
                let mut spl_line: Vec<i32> =
                    moo.iter().map(|x| -> i32 { x.parse().unwrap() }).collect();
                for t in 0..patch_sz {
                    let pix_val: i32 = spl_line.pop().unwrap();
                    let pix_id: i32 = spl_line.pop().unwrap();
                    let pid = pix_map[pix_id as usize];
                    if pid < 0 {
                        buffer = String::new();
                        buf = fp.read_line(&mut buffer)?;
                        continue;
                    } else {
                        data_frame[pid as usize] = -1;
                    }
                }
            } else {
                buffer = String::new();
                fp.read_line(&mut buffer)?;
                // let moo: Vec<&str> = buffer.trim().split_whitespace().collect();
                // let mut spl_line: Vec<i32> =
                //     moo.iter().map(|x| -> i32 { x.parse().unwrap() }).collect();

                // for t in 0..patch_sz {
                //     let pix_val: i32 = spl_line.pop().unwrap();
                //     let pix_id: i32 = spl_line.pop().unwrap();
                // }
            }

            buffer = String::new();
            buf = fp.read_line(&mut buffer)?;
        }
        let ds_name = format!("data_{:06}", d);
        let ds: &str = &ds_name;
        data.new_dataset_builder()
            .with_data(&data_frame)
            .create(ds)
            .expect("Dataset write failed");
        // let t3 = t2.elapsed();
        // println!("Time elapsed: {:?}", t3);
    }

    let t4 = t1.elapsed();
    println!("Time elapsed: {:?}", t4);
    Ok(())
}

fn load_cbf(filepath: &str) -> Result<Vec<i32>, Box<dyn Error>> {
    // let cbf: Box<dyn Frame> = frame::open(filepath).expect("Could not open CBF file {filepath}");
    let cbf = match Cbf::read_file(filepath) {
        Ok(cbf) => cbf,
        Err(why) => return Err(format!("Could not open file {}: {}",filepath, why).into()),
    };
    let framedata = cbf.array().data();
    let x = framedata.iter().map(|z| *z as i32).collect();
    Ok(x)
}
