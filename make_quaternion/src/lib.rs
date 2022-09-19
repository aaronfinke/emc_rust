use clap::Parser;
use std::path::PathBuf;
#[derive(Default, Clone, Copy)]
pub struct Qpoint {
    pub vec: [[i32; 2]; 4],
    pub weight: f64,
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    /// Do icosahedral quaternions
    #[clap(long, parse(from_flag))]
    pub ico: bool,

    /// Save as binary file
    #[clap(long, parse(from_flag))]
    pub bin: bool,

    /// Save as HDF5 file
    #[clap(long, parse(from_flag))]
    pub hdf5: bool,

    /// Path to config_file
    #[clap(value_parser, value_name = "FILE")]
    pub config: PathBuf,

    /// Number of divisions
    #[clap(value_parser)]
    pub num_divs: i32,
}
