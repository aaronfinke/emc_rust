use crate::config::configure::Config;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Seek;
use std::io::SeekFrom;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;
pub mod config;

/*
reduce data and balance the work loads between processors 

usage:
reduce_data [config file]

needs:
config_file, quat_file, cbflist_file, radiallist_file, orienlist_file,
peaklist_file, 
num_prob_orien.dat (in data_dir),
individual ave_bg, peak, orien files (in radial_bg and outlier dirs),
raw data files (in raw_data dir)

makes:
prob_orien_file, mpi_bgfile, reduced_cbflist_file, reduced_radiallist_file,
reduced_peaklist_file, reduced_data_id_file (usually in data_dir)

*/

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
    let num_div: i32 = config.quat_file.split("quaternion").collect::<Vec<&str>>()[1]
        .split(".")
        .collect::<Vec<&str>>()[0]
        .parse::<i32>()
        .unwrap();
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

    let mut num_orien: Vec<i32> = vec![0; config.num_raw_data];
    let mut prob_orien: Vec<Vec<i32>> = vec![vec![]; config.num_raw_data];

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
        if num_orien[d] == 0 {
            continue;
        }
        let fp = match File::open(&orient_files[d]) {
            Err(why) => panic!("Couldn't open orien_file {} : {}", orient_files[d], why),
            Ok(file) => BufReader::new(file).lines(),
        };
        let mut probs: Vec<i32> = fp
            .into_iter()
            .map(|x| x.expect("not a number").parse::<i32>().unwrap())
            .collect();
        let num: i32 = probs.remove(0);
        if num != num_orien[d] as i32 {
            println!("num = {}, num_orien[{d}] = {}", num, num_orien[d]);
        }
        prob_orien[d] = probs;
    }
    let num_prov_data = num_orien.iter().filter(|&i| *i > 0).count();
    println!(
        "num_raw_data = {}, num_data = {num_prov_data}",
        config.num_raw_data
    );
    println!("num_rot = {num_rot}, num_div = {num_div}");

    /* End setup section */

    /* Start load_balance section */

    let mut count_d2r: Vec<i32> = vec![];
    let mut data_id: Vec<i32> = vec![];

    let mut num_data = 0;
    for d in 0..config.num_raw_data {
        if num_orien[d] > 0 {
            count_d2r.push(num_orien[d]);
            data_id.push(d as i32);
            num_data += 1;
        }
    }

    (count_d2r, data_id) = merge_sort(count_d2r, data_id);

    let s_num_data: i32 = match num_data % config.nproc {
        0 => num_data / config.nproc,
        _ => num_data / config.nproc + 1,
    };
    let mut frame_count: Vec<i32> = vec![s_num_data - 1; config.nproc as usize];
    let mut proc_num_orien: Vec<i32> = vec![0; config.nproc as usize];
    let mut proc_occupancy: Vec<Vec<i32>> = vec![vec![]; config.nproc as usize];
    let mut proc_data_id: Vec<Vec<i32>> = vec![vec![]; config.nproc as usize];
    for i in 0..config.nproc as usize {
        proc_occupancy[i] = vec![0; num_rot as usize];
        proc_data_id[i] = vec![0; s_num_data as usize];
    }
    let mut proc_id: Vec<i32> = (0..config.nproc).collect();
    for d in 0..s_num_data - 1 {
        let idx_offset = d * config.nproc;
        (proc_num_orien, proc_id) = merge_sort(proc_num_orien, proc_id);
        for i in 0..config.nproc {
            let idx: usize = (idx_offset + i) as usize;
            let count_d2r_d = count_d2r[idx];
            let prob_orien_d = &prob_orien[data_id[idx] as usize];
            let j: usize = (config.nproc - 1 - i) as usize;
            for r in 0..count_d2r_d as usize {
                let rid: usize = prob_orien_d[r] as usize;
                if proc_occupancy[proc_id[j] as usize][rid] == 0 {
                    proc_occupancy[proc_id[j] as usize][rid] = 1;
                    proc_num_orien[j] += 1;
                }
            }
            proc_data_id[proc_id[j] as usize][d as usize] = data_id[idx];
        }
    }

    (proc_num_orien, proc_id) = merge_sort(proc_num_orien, proc_id);
    let idx_offset = (s_num_data - 1) * config.nproc;
    for d in idx_offset..num_data {
        let j: usize = (config.nproc - 1 - d % config.nproc) as usize;
        let count_d2r_d = count_d2r[d as usize];
        let prob_orien_d = &prob_orien[data_id[d as usize] as usize];
        for r in 0..count_d2r_d as usize {
            let rid: usize = prob_orien_d[r] as usize;
            if proc_occupancy[proc_id[j] as usize][rid] == 0 {
                proc_occupancy[proc_id[j] as usize][rid] = 1;
                proc_num_orien[j] += 1;
            }
        }
        proc_data_id[proc_id[j] as usize][(s_num_data - 1) as usize] = data_id[d as usize];
        frame_count[j] += 1;
    }

    let mut load_balanced_data_id: Vec<i32> = vec![0; num_data as usize];
    num_data = 0;
    for i in 0..config.nproc {
        let j: usize = (config.nproc - 1 - i) as usize;
        for d in 0..frame_count[j] as usize {
            let idx = proc_data_id[proc_id[j] as usize][d];
            load_balanced_data_id[num_data as usize] = idx;
            num_data += 1;
        }
        println!("myid = {i}, frame_count = {}", frame_count[j]);
    }

    /* end load_balance section */
    /* start write_files section */

    {
        let prob_orien_file = PathBuf::from(&config.prob_orien_file);
        let mut fp = File::create(prob_orien_file)?;
        for d in 0..num_data {
            let idx: usize = load_balanced_data_id[d as usize] as usize;
            fp.write(&num_orien[idx].to_le_bytes())?;
            prob_orien[idx].iter().for_each(|x| {
                fp.write(&x.to_le_bytes())
                    .expect("Could not write to prob_orien_file");
            });
        }
    }
    {
        let reduced_cbflist_file = PathBuf::from(&config.reduced_cbflist_file);
        let mut fp = File::create(reduced_cbflist_file)?;
        for d in 0..num_data as usize {
            let idx = load_balanced_data_id[d];
            writeln!(fp, "{}", cbf_files[idx as usize])?;
        }
    }
    {
        let reduced_radiallist_file = PathBuf::from(&config.reduced_radiallist_file);
        let mut fp = File::create(reduced_radiallist_file)?;
        for d in 0..num_data as usize {
            let idx = load_balanced_data_id[d];
            writeln!(fp, "{}", radialfiles[idx as usize])?;
        }
    }
    {
        let reduced_peaklist_file = PathBuf::from(&config.reduced_peaklist_file);
        let mut fp = File::create(reduced_peaklist_file)?;
        for d in 0..num_data as usize {
            let idx = load_balanced_data_id[d];
            writeln!(fp, "{}", peak_files[idx as usize])?;
        }
    }
    {
        let mpi_bgfile = PathBuf::from(&config.mpi_bgfile);
        let mut fp = File::create(mpi_bgfile)?;
        for d in 0..num_data as usize {
            let idx = load_balanced_data_id[d] as usize;
            let radialfile = PathBuf::from(&radialfiles[idx]);
            let fp2 = File::open(radialfile).expect("could not open radialfile");
            let mut radial_file_buff: BufReader<File> = BufReader::new(fp2);
            for _ in 0..config.qlen {
                let mut buff = vec![0_u8; 8];
                radial_file_buff.read_exact(&mut buff)?;
                let qval = f64::from_le_bytes(buff.try_into().unwrap());
                fp.write(&qval.to_le_bytes())?;
            }
        }
    }
    {
        let reduced_data_id_file = PathBuf::from(&config.reduced_data_id_file);
        let mut fp = File::create(reduced_data_id_file)?;
        writeln!(fp, "{num_data}")?;
        for d in 0..num_data as usize {
            writeln!(fp, "{}", load_balanced_data_id[d])?;
        }
    }


    let t2 = t1.elapsed();
    println!("Time elapsed: {:?}", t2);
    Ok(())
}

/// sort vectors <i>A</i> and <i>B</i> by vector <i>A</i>
fn merge_sort(a: Vec<i32>, b: Vec<i32>) -> (Vec<i32>, Vec<i32>) {
    let mut zipped = a
        .clone()
        .into_iter()
        .zip(b.clone().into_iter())
        .collect::<Vec<(i32, i32)>>();
    zipped.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).expect("Can't compare"));
    let (a_sorted, b_sorted) = zipped.into_iter().unzip();
    return (a_sorted, b_sorted);
}
