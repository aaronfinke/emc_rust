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
    pub outlier_files: Vec<String>,
    pub peak_files: Vec<String>,
    pub qmax: i32,
    pub qmin: i32,
    pub total_pix: i32,
    pub min_patch_sz: i32, 
    pub max_patch_sz: i32, 
    pub min_num_peak: i32, 
    pub max_num_peak: i32, 
    pub res_cutoff: f64, 
    pub VN: i32, 
    pub gw: f64, 
    pub quat: Vec<Quat>,
    pub tolmult: f64,
    pub delphi: f64,
    pub beam_vec: [f64;3],
    pub rmax2: f64,
    pub rmin2: f64,
    pub pix: Vec<Pix>,
    pub cell: [f64;6]
    
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
        let Rstop = config.get("make_detector", "Rstop").unwrap().parse::<f64>()?;
        let num_data: usize = config.get("make_background", "num_raw_data")
            .unwrap()
            .parse::<usize>()?;
        let hot_pix_thres: i32 = config.get("make_background", "hot_pix_thres")
            .unwrap()
            .parse::<f64>()? as i32;
        let min_patch_sz = config.get("make_powder", "min_patch_sz").unwrap().parse::<i32>()?;
        let max_patch_sz = config.get("make_powder", "max_patch_sz").unwrap().parse::<i32>()?;
        let min_num_peak = config.get("make_powder", "min_num_peak").unwrap().parse::<i32>()?;
        let max_num_peak = config.get("make_powder", "max_num_peak").unwrap().parse::<i32>()?;

        let res_cutoff = config.get("orient_peak", "res_cutoff").unwrap().parse::<f64>()?;
        let VN = config.get("orient_peak", "VN").unwrap().parse::<i32>()?;
        let gw = config.get("orient_peak", "gw").unwrap().parse::<f64>()?;
        let tolmult = config.get("orient_peak", "tolerance_mult").unwrap().parse::<f64>()?;
        let delphi = config.get("orient_peak", "delta_phi").unwrap().parse::<f64>()?;
        let mut beam_vec: [f64;3] = [
            config.get("make_detector", "sx").unwrap().parse::<f64>()?,
            config.get("make_detector", "sy").unwrap().parse::<f64>()?,
            config.get("make_detector", "sz").unwrap().parse::<f64>()?
        ];


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
        let quat_file = match File::open(quatfile) {
            Err(why) => panic!("Couldn't open quatfile : {}", why),
            Ok(file) => file,
        };
        let mut quat_file_buff: BufReader<File> = BufReader::new(quat_file);
        quat_file_buff.seek(SeekFrom::Start(0)).expect("moo");
        let mut buf = vec![0; 4];
        quat_file_buff.read_exact(&mut buf).expect("moo");
        let num_rot = u32::from_le_bytes(buf.try_into().unwrap());
        println!("num_rot = {num_rot}");
        let mut buf = vec![0; 4];
        quat_file_buff.read_exact(&mut buf).expect("moo");
        let buffer = u32::from_le_bytes(buf.try_into().unwrap());
        
        let mut quat : Vec<Quat> = vec![Quat(0.,0.,0.,0.,0.); num_rot as usize];    
        for r in 0..num_rot as usize {
            let mut buf = vec![0; 8];
            quat_file_buff.read_exact(&mut buf).expect("moo");
            quat[r].0 = f64::from_le_bytes(buf.try_into().unwrap());
            buf = vec![0; 8];
            quat_file_buff.read_exact(&mut buf).expect("moo");
            quat[r].1 = f64::from_le_bytes(buf.try_into().unwrap());
            buf = vec![0; 8];
            quat_file_buff.read_exact(&mut buf).expect("moo");
            quat[r].2 = f64::from_le_bytes(buf.try_into().unwrap());
            buf = vec![0; 8];
            quat_file_buff.read_exact(&mut buf).expect("moo");
            quat[r].3 = f64::from_le_bytes(buf.try_into().unwrap());
            buf = vec![0; 8];
            quat_file_buff.read_exact(&mut buf).expect("moo");
            quat[r].4 = f64::from_le_bytes(buf.try_into().unwrap());
        }
        let elapsed = now.elapsed();
        println!("Elapsed time: {:?}",elapsed);
        let norm = (beam_vec[0].powi(2) + beam_vec[1].powi(2) + beam_vec[2].powi(2)).sqrt();
        let D = (detd/px) * norm/(beam_vec[2]).abs();
        let rescale = (1.0/wl)/D;

        beam_vec[0] *= D/norm ;
        beam_vec[1] *= D/norm ;
        beam_vec[2] *= D/norm ;

        let qmax = (wl*D/res_cutoff).ceil() as i32;
        let qmin = (2.0*D*sin(atan(Rstop/D)/2.0)).ceil() as i32;
        let rmax2 = (qmax as f64*rescale).powi(2);
        let rmin2 = (qmin as f64*rescale).powi(2);
        let total_pix: i32 = num_row*num_col;

        /* pix[t] has unit 1/A */
        let mut pix: Vec<Pix> = vec![Pix(0.,0.,0.); total_pix as usize];
        for t in 0..total_pix {
            let x = (t/num_col) as f64 - cx;
            let y = (t % num_col) as f64 - cy;
            let mut sx = beam_vec[0] + x ;
            let mut sy = beam_vec[1] + y ;
            let mut sz = beam_vec[2] ;
            let norm = sqrt(sx*sx + sy*sy + sz*sz);
            sx *= D/norm ;
            sy *= D/norm ;
            sz *= D/norm ;
            pix[t as usize].0 = (sx - beam_vec[0])*rescale ;
            pix[t as usize].1 = (sy - beam_vec[1])*rescale ;
            pix[t as usize].2 = (sz - beam_vec[2])*rescale ;
        
        }

        let mut unit_cell_buff = match File::open(&unit_cell_file) {
            Err(why) => panic!("Couldn't open unit cell file : {}", why),
            Ok(file) => BufReader::new(file),
        };
        let mut first_line = String::new();
        unit_cell_buff.read_line(&mut first_line)?;
        let mut ucline: Vec<&str> = first_line.trim().split_whitespace().collect();
        ucline.pop();
        let unit_cell: Vec<f64> = ucline.iter().map(|word| word.parse::<f64>().unwrap()).collect();
        let cell: [f64;6] = unit_cell.try_into().unwrap_or_else(
            |_| panic!("Unit Cell wrong format: should be a b c al be gam sg"));

        


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
            outlier_files,
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
            cell
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
