use configparser::ini::Ini;
use hdf5_bitshuffle::register_bitshuffle_plugin;
use hdf5::{File as H5File};
use ndarray::{Array1};
use std::{time::Instant};
use std::path::Path;
use std::env;

fn main() {
    let now = Instant::now();
    register_bitshuffle_plugin();
    let mut config = Ini::new();
    let args: Vec<String> = env::args().collect();
    match args.len() {
        2 => {},
        _ => panic!("execution: ./make_detector (/path/to/)config.ini")
    };
    let config_path = &args[1];
    let _map = match config.load(config_path) {
        Err(_) => panic!("Config file not found"),
        Ok(cfg) => cfg
    };
    let num_row = config.get("make_detector","num_row").unwrap()
                    .parse::<i32>().unwrap();
    let num_col = config.get("make_detector","num_col").unwrap()
                    .parse::<i32>().unwrap();
    let detd = config.get("make_detector","detd").unwrap()
                    .parse::<f64>().unwrap();
    let wl = config.get("make_detector","wl").unwrap()
                    .parse::<f64>().unwrap();
    let px = config.get("make_detector","px").unwrap()
                    .parse::<f64>().unwrap();
    let cx = config.get("make_detector","cx").unwrap()
                    .parse::<f64>().unwrap();
    let cy = config.get("make_detector","cy").unwrap()
                    .parse::<f64>().unwrap();
    let sx = config.get("make_detector","sx").unwrap()
                    .parse::<f64>().unwrap();
    let sy = config.get("make_detector","sy").unwrap()
                    .parse::<f64>().unwrap();
    let sz = config.get("make_detector","sz").unwrap()
                    .parse::<f64>().unwrap();
    let res_max = config.get("make_detector","res_max").unwrap().parse::<f64>().unwrap();
    let rstop = config.get("make_detector","rstop").unwrap().parse::<f64>().unwrap();
    let mask_file = config.get("files","mask_file").unwrap();
    let pixmap_file = config.get("files","pixmap_file").unwrap();
    let rec2pix_file = config.get("files","rec2pix_file").unwrap();
    let rec_vectors_file = config.get("files","rec_vectors_file").unwrap();

    let total_pixel = num_row * num_col;
    let mask = read_hdf5mask(mask_file);
    let mut beam_vec = vec![sx,sy,sz];
    let mut pix_map:Vec<i32> = std::iter::repeat(-1).take(total_pixel.try_into().unwrap()).collect::<Vec<_>>();
    let mut norm = 0.0;
    for i in &beam_vec {
        norm += i * i;
    }
    norm = (norm).sqrt();
    let d_space = (detd / px) * norm/(beam_vec[2]).abs();
    for i in &mut beam_vec {
        *i *= d_space/norm;
    }
    let qmax : i32 = (wl*d_space/res_max).ceil() as i32;
    let theta = (rstop/d_space).atan();
    let qmin: i32 = (2.0*d_space*(theta/2.0).sin()).ceil() as i32;

    let mut npix:i32 = 0;
    let mut max_scale_factor = 0.0;
    let mut qval_arr: Vec<f64> =   std::iter::repeat(0.0).take(total_pixel.try_into().unwrap()).collect::<Vec<_>>();
    let mut rec2pix_map:Vec<i32> = std::iter::repeat(-1).take(total_pixel.try_into().unwrap()).collect::<Vec<_>>();

    for i in 0..total_pixel {
        let x:f64 =(i/num_col) as f64 - cx;
        let y:f64 = (i % num_col) as f64 - cy;
        let mut sx = beam_vec[0] + (x as f64); 
        let mut sy = beam_vec[1] + (y as f64);
        let mut sz = beam_vec[2];
        let mut norm = (sx*sx + sy*sy + sz*sz).sqrt();
        let polarization = 1.0 - (sy/norm).powi(2);
        let solid_angle = sz.abs()/(norm).powi(3);
        let scale_factor = polarization * solid_angle;

        sx *= d_space / norm;
        sy *= d_space / norm;
        sz *= d_space / norm;
        let qx = sx - beam_vec[0]; 
        let qy = sy - beam_vec[1];
        let qz = sz - beam_vec[2];
       norm = (qx*qx + qy*qy + qz*qz).sqrt();
       qval_arr[i as usize] = norm;
        if mask[i as usize] == 0{
            if (qmin as f64) <= norm && norm <= (qmax as f64) {
                npix += 1;
                if scale_factor > max_scale_factor {
                    max_scale_factor = scale_factor;
                }
                rec2pix_map[i as usize] = i;
            }
        }
    }
    merge_sort(&mut qval_arr, &mut rec2pix_map);

    let path1 = Path::new(&rec_vectors_file);
    let display1 = path1.display();
    let path2 = Path::new(&rec2pix_file);
    let display2 = path2.display();
    let path3 = Path::new(&pixmap_file);
    let display3 = path3.display();

    let out_fp = match H5File::create(&path1) {
        Err(why) => panic!("couldn't open {}: {}", display1, why),
        Ok(file) => file,
    };
    let out_fp2 = match H5File::create(&path2) {
        Err(why) => panic!("couldn't open {}: {}", display2, why),
        Ok(file) => file,
    };
    let out_fp3 = match H5File::create(&path3) {
        Err(why) => panic!("couldn't open {}: {}", display3, why),
        Ok(file) => file,
    };
    let mut count = 0;
    let mut qx_v:Vec<f64> = vec![];
    let mut qy_v:Vec<f64> = vec![];
    let mut qz_v:Vec<f64> = vec![];
    let mut scale_factor_v:Vec<f64> = vec![];
    let mut rec2pix_map_f:Vec<i32> = vec![];
    

    for t in 0..total_pixel as usize {
        if rec2pix_map[t] < 0 {continue};
        let x = (((rec2pix_map[t]/num_col) as i32) as f64) - cx;
        let y: f64 = (rec2pix_map[t] % num_col) as f64 - cy;
        let mut sx = beam_vec[0] + (x as f64); 
        let mut sy = beam_vec[1] + (y as f64);
        let mut sz = beam_vec[2];

        let norm = (sz*sx + sy*sy + sz*sz).sqrt();
        let polarization = 1.0 - (sy/norm).powi(2);
        let solid_angle = sz.abs()/(norm).powi(3);
        let scale_factor = polarization*solid_angle/max_scale_factor;
        sx *= d_space/norm ;
        sy *= d_space/norm ;
        sz *= d_space/norm ;
        let qx = sx - beam_vec[0] ;
        let qy = sy - beam_vec[1] ;
        let qz = sz - beam_vec[2] ;
        qx_v.push(qx);
        qy_v.push(qy);
        qz_v.push(qz);
        scale_factor_v.push(scale_factor);
        pix_map[rec2pix_map[t] as usize] = count;
        count += 1;
        rec2pix_map_f.push(rec2pix_map[t]);
        }
    let group = out_fp.create_group("data").expect("Can't create group");
    group.new_dataset_builder()
        .with_data(&qx_v)
        .create("qx")
        .expect("Dataset write failed");
    group.new_dataset_builder()
        .with_data(&qy_v)
        .create("qy")
        .expect("Dataset write failed");
    group.new_dataset_builder()
        .with_data(&qz_v)
        .create("qz")
        .expect("Dataset write failed");
    group.new_dataset_builder()
        .with_data(&scale_factor_v)
        .create("scale_factor")
        .expect("Dataset write failed");
    group.new_attr::<i32>().create("qmin")
        .expect("Dataset write failed")
        .write_scalar(&qmin)
        .expect("Dataset write failed");
    group.new_attr::<i32>().create("qmax")
        .expect("Dataset write failed")
        .write_scalar(&qmax)
        .expect("Dataset write failed");
    group.new_attr::<i32>().create("Npix")
        .expect("Dataset write failed")
        .write_scalar(&npix)
        .expect("Dataset write failed");
    
    let group2 = out_fp2.create_group("data").expect("Group write failed");
    group2.new_dataset_builder()
        .with_data(&rec2pix_map_f)
        .create("rec2pix_map")
        .expect("Dataset write failed");
    
    let group3 = out_fp3.create_group("data").expect("Failed to create group");
    group3.new_dataset_builder()
        .with_data(&pix_map)
        .create("pix_map")
        .expect("Dataset write failed");


    println!("cx = {:.1}, cy = {:.1}, Rstop = {:.1}",cx, cy, rstop);
    println!("qmax = {qmax}, qmin = {qmin}, num_pix = {npix}");
    let elapsed = now.elapsed();
    println!("Time elapsed: {:?}",elapsed);

}



fn read_hdf5mask(mask_file: String) -> Array1<i32> {
    let file = H5File::open(mask_file).unwrap(); 
    let mask = file.dataset("/mask").unwrap();
    let mask1 = mask.read_1d::<i32>().unwrap();

    mask1
}

fn merge_sort<T: PartialOrd + Copy,L:PartialOrd + Copy>(input: &mut [T], input2: &mut [L]) {
    if input.len() < 2 {return;}
    
    let len = input.len();
    let mid = len / 2;
    merge_sort(&mut input[..mid], &mut input2[..mid]);
    merge_sort(&mut input[mid..], &mut input2[mid..]);

    let mut tmp = Vec::with_capacity(len);
    let mut tmp2:Vec<L> = Vec::with_capacity(len);
    
    let mut i = 0;
    let mut j = mid;

    while i < mid && j < len {
        if input[i] < input[j] {
            tmp.push(input[i]);
            tmp2.push(input2[i]);
            i += 1;
        } else {
            tmp.push(input[j]);
            tmp2.push(input2[j]);
            j += 1;
        }
    }
    if i < mid {
        tmp.extend_from_slice(&input[i..mid]);
        tmp2.extend_from_slice(&input2[i..mid]);
    } else if j < len {
        tmp.extend_from_slice(&input[j..len]);
        tmp2.extend_from_slice(&input2[j..len]);
    }

    input.copy_from_slice(&tmp[..]);
    input2.copy_from_slice(&tmp2[..]);
}
