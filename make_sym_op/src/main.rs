use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::error;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::process;
use std::io::Write;
#[derive(Serialize, Deserialize, Debug)]
struct Spacegroup {
    sg_name: String,
    sg_num: i32,
    n_symops: i32,
    symops: Vec<Vec<f32>>,
    is_centro: bool,
}

fn main() -> Result<(), Box<dyn error::Error>> {
    let args: Vec<String> = env::args().collect();
    let mut unit_cell_file = "";
    let mut symops_file = env::current_dir()?;
    match args.len() {
        2 => {
            if &args[1] == "-h" {
                println!(
                "Usage: make_sym_op [symops file] [unit cell file]");
                println!(
                "unit_cell file should have one line:"
                );
                println!(
                "a b c alpha beta gamma sgname"
                );
                println!("sgname should be all lowercase no spaces");
                process::exit(0);
            }
            unit_cell_file = &args[1];
            symops_file.push("symops.json");
        }
        3 => {
            unit_cell_file = &args[2];
            symops_file = PathBuf::from(&args[1]);
        }
        _ => return Err("Usage: make_sym_op [symops file] [unit cell file]".into()),
    };
    let sfile = File::open(&symops_file)?;
    let sreader = io::BufReader::new(sfile);
    let spacegroups: HashMap<String, Spacegroup> = serde_json::from_reader(sreader).unwrap();

    let first_line: String = {
        let unit_cell = File::open(unit_cell_file)?;
        let mut ureader = io::BufReader::new(unit_cell);
        let mut first_line = String::new();
        ureader.read_line(&mut first_line)?;
        first_line
    };
    let unit_cell: Vec<&str> = first_line.trim().split(" ").collect();
    let sg = unit_cell[unit_cell.len()-1];
    let curr_sg: &Spacegroup = match spacegroups.get(sg){
        Some(_) => spacegroups.get(sg).unwrap(),
        None => return Err(format!("Space group not found: {:?}", sg).into())
    };
    // get symops
    let symops = &curr_sg.symops;
    let mut symops_vec: Vec<String> = vec![];
    for sym in symops {
        let mut t: String = String::new();
        for x in sym {
            let y = format!("{:.1} ",x);
            t.push_str(&y.to_string())
        }
        symops_vec.push(t);
    }
    //generate inverse symops for non-centric SGs
    if !curr_sg.is_centro {
        for sym in symops {
            let mut t: String = String::new();
            let mut f = sym.clone();
            f.iter_mut().for_each(|x| {
                if *x != 0.0 {*x *= -1.0}});
            for x in f {
                let y = format!("{:.1} ",x);
                t.push_str(&y.to_string())
            }
            symops_vec.push(t);
        }
    }

    let num_symops: usize = symops_vec.len();
    let mut fp = match File::create(&unit_cell_file) {
        Err(why) => return Err(format!("couldn't create {} {}", &unit_cell_file, why).into()),
        Ok(file) => file,
    };
    writeln!(&mut fp, "{}",first_line.trim())?;
    writeln!(&mut fp, "{}",num_symops)?;
    for symop in symops_vec {
        writeln!(&mut fp, "{}",symop)?;
    }

    println!("{}", first_line.trim());
    println!("symops written to unit_cell.dat");
    Ok(())
}
