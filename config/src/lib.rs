use configparser::ini::Ini;
use std::error::Error;


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
    pub num_data: usize,
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
    pub high_res_cutoff: i32,
    pub beam_vec: [f64;3],

}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, Box<dyn Error>> {
        let mut config = Ini::new();
        let config_file = &args[1];
        config.load(config_file)?;
        //files and directories
        let home_dir = config.get("files","home_dir").unwrap();
        let raw_data_dir = config.get("files","raw_data_dir").unwrap();
        let required_string = config.get("files","required_string").unwrap();
        let reduced_data_dir = config.get("files","reduced_data_dir").unwrap();
        let data_dir = config.get("files","data_dir").unwrap();
        let prob_dir = config.get("files","prob_dir").unwrap();
        let radial_bg_dir = config.get("files","radial_bg_dir").unwrap();
        let outlier_dir = config.get("files","outlier_dir").unwrap();
        let start_prob_dir = config.get("files","start_prob_dir").unwrap();
        let mask_file = config.get("files","mask_file").unwrap();
        let pixmap_file = config.get("files","pixmap_file").unwrap();
        let rec2pix_file = config.get("files","rec2pix_file").unwrap();
        let rec_vectors_file = config.get("files","rec_vectors_file").unwrap();
        let unit_cell_file = config.get("files","unit_cell_file").unwrap();
        let quat_file = config.get("files","quat_file").unwrap();
        let start_phi_file = config.get("files","start_phi_file").unwrap();
        let start_intens_file = config.get("files","start_intens_file").unwrap();
        let cbflist_file = config.get("files","cbflist_file").unwrap();
        let radiallist_file = config.get("files","radiallist_file").unwrap();
        let orienlist_file = config.get("files","orienlist_file").unwrap();
        let outlierlist_file = config.get("files","outlierlist_file").unwrap();
        let peaklist_file = config.get("files","peaklist_file").unwrap();
        let data_info_file = config.get("files","data_info_file").unwrap();
        let reduced_data_id_file = config.get("files","reduced_data_id_file").unwrap();
        let reduced_cbflist_file = config.get("files","reduced_cbflist_file").unwrap();
        let reduced_radiallist_file = config.get("files","reduced_radiallist_file").unwrap();
        let reduced_peaklist_file = config.get("files","reduced_peaklist_file").unwrap();
        let r2peak_file = config.get("files","r2peak_file").unwrap();
        let peak2r_file = config.get("files","peak2r_file").unwrap();
        let mpi_bgfile = config.get("files","mpi_bgfile").unwrap();
        let mpi_datafile = config.get("files","mpi_datafile").unwrap();
        let prob_orien_file = config.get("files","prob_orien_file").unwrap();
        let num_prob_orien_file = config.get("files","num_prob_orien_file").unwrap();
        let fine_quat_file = config.get("files","fine_quat_file").unwrap();
        let quat_table_file = config.get("files","quat_table_file").unwrap();
        let local_r2peak_file = config.get("files","local_r2peak_file").unwrap();
        let local_peak2r_file = config.get("files","local_peak2r_file").unwrap();
        // make_detector
        let file_type = config.get("make_detector", "file_type").unwrap();
        let num_row = config.get("make_detector", "num_row").unwrap().parse::<i32>()?;
        let num_col = config.get("make_detector", "num_col").unwrap().parse::<i32>()?;
        let detd = config.get("make_detector", "detd").unwrap().parse::<f64>()?;
        let wl = config.get("make_detector", "wl").unwrap().parse::<f64>()?;
        let px = config.get("make_detector", "px").unwrap().parse::<f64>()?;
        let cx = config.get("make_detector", "cx").unwrap().parse::<f64>()?;
        let cy = config.get("make_detector", "cy").unwrap().parse::<f64>()?;
        let bh_width = config.get("make_detector", "bh_width").unwrap().parse::<i32>()?;
        let Rstop = config.get("make_detector", "Rstop").unwrap().parse::<f64>()?;
        let res_max = config.get("make_detector", "res_max").unwrap().parse::<f64>()?;
        let beam_vec: [f64; 3] = [
            config.get("make_detector", "sx").unwrap().parse::<f64>()?,
            config.get("make_detector", "sy").unwrap().parse::<f64>()?,
            config.get("make_detector", "sz").unwrap().parse::<f64>()?,
        ];
        // make_background
        let num_data: usize = config.get("make_background", "num_raw_data").unwrap().parse::<usize>()?;
        let hot_pix_thres: i32 = config.get("make_background", "hot_pix_thres").unwrap().parse::<f64>()? as i32;
        let qlen: i32 = config.get("make_background", "qlen").unwrap().parse::<i32>()?;
        // make_powder
        let min_patch_sz = config.get("make_powder", "min_patch_sz").unwrap().parse::<i32>()?;
        let max_patch_sz = config.get("make_powder", "max_patch_sz").unwrap().parse::<i32>()?;
        let min_num_peak = config.get("make_powder", "min_num_peak").unwrap().parse::<i32>()?;
        let max_num_peak = config.get("make_powder", "max_num_peak").unwrap().parse::<i32>()?;
        // orient_peak
        let res_cutoff = config.get("orient_peak", "res_cutoff").unwrap().parse::<f64>()?;
        let VN = config.get("orient_peak", "VN").unwrap().parse::<i32>()?;
        let gw = config.get("orient_peak", "gw").unwrap().parse::<f64>()?;
        let tolmult = config.get("orient_peak", "tolerance_mult").unwrap().parse::<f64>()?;
        let delphi = config.get("orient_peak", "delta_phi").unwrap().parse::<f64>()?;
        //reduce_data
        let nproc = config.get("reduce_data", "nproc").unwrap().parse::<i32>()?;
        //low_res_emc
        let iter_data_block = config.get("low_res_emc", "iter_data_block").unwrap().parse::<i32>()?;
        // high_res_emc
        let high_res_cutoff = config.get("high_res_emc", "high_res_cutoff").unwrap().parse::<i32>()?;


        //  read in outlier files
        


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
            num_data,
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
}
