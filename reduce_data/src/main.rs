use crate::config::configure::Config;
use std::error::Error;
use std::time::Instant;
use std::env;
use std::io::BufReader;
use std::io::SeekFrom;
use std::fs::File;
use std::io::Seek;
use std::io::Read;
use std::io::BufRead;

pub mod config;
fn main() -> Result<(), Box<dyn Error>> {
    let t1 = Instant::now();

    let args: Vec<String> = env::args().collect();
    match args.len() {
        2 => {}
        _ => return Err("execution: ./make_background (/path/to/)config.ini".into()),
    };

    let config = match Config::new(&args){
        Ok(config) => config,
        Err(why) => return Err(format!("Error in config.ini: {:?}",why).into() )
    };

    let num_rot: i32 = {
        let quat_file = match File::open(&config.quat_file) {
            Err(why) => panic!("Couldn't open quatfile : {}", why),
            Ok(file) => file,
        };    
        let mut quat_file_buff: BufReader<File> = BufReader::new(quat_file);
        quat_file_buff.seek(SeekFrom::Start(0))?;
        let mut buf = vec![0; 4];
        quat_file_buff.read_exact(&mut buf)?;
        i32::from_le_bytes(buf.try_into().unwrap())
    };
    let num_div: i32 = config.quat_file.split("quaternion")
                        .collect::<Vec<&str>>()[1].split(".")
                        .collect::<Vec<&str>>()[0]
                        .parse::<i32>().unwrap();
    //  read in outlier files
    let mut orient_files: Vec<String> = vec![Default::default(); config.num_raw_data];
    let orient_file_buff = match File::open(&config.orienlist_file) {
        Err(why) => panic!("Couldn't open outlierlist_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for (i, line) in orient_file_buff.enumerate() {
        if let Ok(ip) = line {
            orient_files[i] = ip;
        }
    }
    //  read in peaklist files
    let mut peak_files: Vec<String> = vec![Default::default(); config.num_raw_data];
    let peak_file_buff = match File::open(&config.peaklist_file) {
        Err(why) => panic!("Couldn't open peaklist_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for (i, line) in peak_file_buff.enumerate() {
        if let Ok(ip) = line {
            peak_files[i] = ip;
        }
    }
    //read in radial files
    let mut radialfiles: Vec<String> = vec![Default::default(); config.num_raw_data];
    let radial_file_buff = match File::open(&config.radiallist_file) {
        Err(why) => panic!("Couldn't open radial_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for (i, line) in radial_file_buff.enumerate() {
        if let Ok(ip) = line {
            radialfiles[i] = ip;
        }
    }
    //read in cbf files
    let mut cbf_files: Vec<String> = vec![Default::default(); config.num_raw_data];
    let cbf_file_buff = match File::open(&config.cbflist_file) {
        Err(why) => panic!("Couldn't open cbf_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for (i, line) in cbf_file_buff.enumerate() {
        if let Ok(ip) = line {
            cbf_files[i] = ip;
        }
    }

    let mut num_orien: Vec<i32> = vec![0;config.num_raw_data];
    let mut prob_orien: Vec<Vec<i32>> = vec![vec![];config.num_raw_data];

    let num_prob_orien_file_buff = match File::open(&config.num_prob_orien_file) {
        Err(why) => panic!("Couldn't open num_prob_orien_file : {}", why),
        Ok(file) => BufReader::new(file).lines(),
    };
    for (i, line) in num_prob_orien_file_buff.enumerate() {
        if let Ok(ip) = line {
            num_orien[i] = ip.parse::<i32>()?;
        }
    }

    for d in 0..config.num_raw_data {
        if num_orien[d] == 0 {continue;}
        let fp = match File::open(&orient_files[d]) {
            Err(why) => panic!("Couldn't open orien_file {} : {}", orient_files[d], why),
            Ok(file) => BufReader::new(file).lines(),
        };
        let mut probs : Vec<i32> = fp.into_iter().map(|x| x.expect("not a number").parse::<i32>().unwrap()).collect();
        let num:i32 = probs.remove(0);
        if num != num_orien[d] as i32 {println!("num = {}, num_orien[{d}] = {}", num, num_orien[d]);}
        prob_orien[d] = probs;
    }
    let num_data = num_orien.iter().filter(|&i| *i > 0).count();
    println!("num_raw_data = {}, num_data = {num_data}", config.num_raw_data);
    println!("num_rot = {num_rot}, num_div = {num_div}");


    Ok(())
}
