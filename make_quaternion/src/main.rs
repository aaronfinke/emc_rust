use clap::Parser;
use configparser::ini::Ini;
use make_quaternion::Args;
use make_quaternion::Qpoint;
use std::error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn error::Error>> {
    const num_vert: usize = 120;
    const num_edge: usize = 720;
    const num_face: usize = 1200;
    const num_cell: usize = 600;
    const nnn: usize = 12;
    const tau: f64 = 1.618033988749895;
    const pi: f64 = 3.141592653589793;
    const f0: f64 = 5.0 / 6.0;
    const f1: f64 = 35.0 / 36.0;

    #[derive(Clone)]
    struct Quat(f64, f64, f64, f64, f64);

    let mut vertices: [[f64; 4]; num_vert] = [[0.0; 4]; num_vert];
    let mut edges: [[i32; 2]; num_edge] = [[0; 2]; num_edge];
    let mut faces: [[i32; 3]; num_face] = [[0; 3]; num_face];
    let mut cells: [[i32; 4]; num_cell] = [[0; 4]; num_cell];
    let mut nn_list: [[i32; nnn]; num_vert] = [[0; nnn]; num_vert];
    let mut edge2cell: [[i32; 4]; num_edge] = [[0; 4]; num_edge];
    let mut face2cell: [[i32; 4]; num_face] = [[0; 4]; num_face];
    let mut vec_vertices: [[[i32; 2]; 4]; num_vert] = [[[0; 2]; 4]; num_vert];

    let args = Args::parse();

    let data_dir = {
        let config_file: PathBuf = args.config;
        let mut config = Ini::new();
        config.load(config_file)?;
        config.get("files", "data_dir").unwrap()
    };
    let num_div: i32 = args.num_divs;
    let is_ico: bool = args.ico;
    let is_bin: bool = args.bin;

    match (is_ico, is_bin) {
        (true, true) => println!("ico option is on, output format: binary"),
        (true, false) => println!("ico option is on!"),
        (false, true) => println!("output format: binary"),
        (false, false) => {}
    }

    // make_vertex section
    let mut idx = 0;
    // 16 vertices
    for h in 0..2 {
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    vertices[idx][0] = h as f64 - 0.5;
                    vertices[idx][1] = i as f64 - 0.5;
                    vertices[idx][2] = j as f64 - 0.5;
                    vertices[idx][3] = k as f64 - 0.5;

                    vec_vertices[idx][0][0] = (2 * h - 1) * num_div;
                    vec_vertices[idx][1][0] = (2 * i - 1) * num_div;
                    vec_vertices[idx][2][0] = (2 * j - 1) * num_div;
                    vec_vertices[idx][3][0] = (2 * k - 1) * num_div;

                    vec_vertices[idx][0][1] = 0;
                    vec_vertices[idx][1][1] = 0;
                    vec_vertices[idx][2][1] = 0;
                    vec_vertices[idx][3][1] = 0;
                    idx += 1;
                }
            }
        }
    }

    for h in 0..2 {
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    vertices[idx][j] = 2. * (h as f64) - 1.0;
                    vec_vertices[idx][j][0] = (2 * h - 1) * 2 * num_div;
                    vec_vertices[idx][j][1] = 0;
                } else {
                    vertices[idx][j] = 0_f64;
                    vec_vertices[idx][j][0] = 0;
                    vec_vertices[idx][j][1] = 0;
                }
            }
            idx += 1;
        }
    }
    // the remaining 96 vertices
    let perm_idx: [[i32; 4]; 12] = [
        [0, 1, 2, 3],
        [0, 2, 3, 1],
        [0, 3, 1, 2],
        [1, 2, 0, 3],
        [1, 0, 3, 2],
        [1, 3, 2, 0],
        [2, 0, 1, 3],
        [2, 3, 0, 1],
        [2, 1, 3, 0],
        [3, 1, 0, 2],
        [3, 0, 2, 1],
        [3, 2, 1, 0],
    ];
    let mut vert: [f64; 4] = [0.0; 4];
    let mut vec_vert: [[i32; 2]; 4] = [[0; 2]; 4];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                vert[0] = (2 * i - 1) as f64 * tau / 2.;
                vert[1] = (2 * j - 1) as f64 * 0.5;
                vert[2] = (2 * k - 1) as f64 / (2_f64 * tau);
                vert[3] = 0_f64;

                vec_vert[0][0] = 0;
                vec_vert[0][1] = (2 * i - 1) * num_div;
                vec_vert[1][0] = (2 * j - 1) * num_div;
                vec_vert[1][1] = 0;
                vec_vert[2][0] = -(2 * k - 1) * num_div;
                vec_vert[2][1] = (2 * k - 1) * num_div;
                vec_vert[3][0] = 0;
                vec_vert[3][1] = 0;

                for m in 0..12 {
                    for n in 0..4 {
                        vertices[idx][n] = vert[perm_idx[m][n] as usize];
                        vec_vertices[idx][n][0] = vec_vert[perm_idx[m][n] as usize][0];
                        vec_vertices[idx][n][1] = vec_vert[perm_idx[m][n] as usize][1];
                    }
                    idx += 1;
                }
            }
        }
    }

    // end make_vertex section

    // start make_edge section
    let epsilon: f64 = 1.0e-6;
    let mut edge_count: i32 = 0;
    let mut nn_count: [i32; num_vert] = [0; num_vert];

    let mut min_dist2: f64 = 0.0;
    for i in 0..4 {
        min_dist2 += (vertices[0][i] - vertices[1][i]).powi(2);
    }
    let mut tmp = 0.0;
    for i in 2..num_vert {
        tmp = 0.0;
        for j in 0..4 {
            tmp += (vertices[0][j] - vertices[i][j]).powi(2);
        }
        if tmp < min_dist2 {
            min_dist2 = tmp;
        }
    }
    //offset by a small number to avoid the round-off error
    min_dist2 += epsilon;

    for i in 0..num_vert {
        for j in i + 1..num_vert {
            tmp = 0.;
            for k in 0..4 {
                tmp += (vertices[i][k] - vertices[j][k]).powi(2);
            }
            if tmp < min_dist2 {
                edges[edge_count as usize][0] = i as i32;
                edges[edge_count as usize][1] = j as i32;
                nn_list[i][nn_count[i] as usize] = j as i32;
                nn_list[j][nn_count[j] as usize] = i as i32;
                nn_count[i] += 1;
                nn_count[j] += 1;
                edge_count += 1;
            }
        }
    }

    //end make_edge section
    //start make_face section

    let mut face_count: i32 = 0;

    for i in 0..num_edge {
        for j in 0..nnn {
            if nn_list[edges[i][0] as usize][j] <= edges[i][1] {
                continue;
            }

            idx = nn_list[edges[i][0] as usize][j] as usize;
            tmp = 0_f64;
            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[edges[i][1] as usize][k]).powi(2);
            }
            if tmp < min_dist2 {
                faces[face_count as usize][0] = edges[i][0];
                faces[face_count as usize][1] = edges[i][1];
                faces[face_count as usize][2] = idx as i32;
                face_count += 1;
            }
        }
    }

    //end make_face section
    //start make_cell section

    let mut cell_count = 0;

    for i in 0..num_face {
        for j in 0..nnn {
            if nn_list[faces[i][0] as usize][j] <= faces[i][2] {
                continue;
            }

            idx = nn_list[faces[i][0] as usize][j] as usize;
            tmp = 0_f64;
            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[faces[i][1] as usize][k]).powi(2);
            }
            if tmp > min_dist2 {
                continue;
            }

            tmp = 0_f64;

            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[faces[i][2] as usize][k]).powi(2);
            }

            if tmp > min_dist2 {
                continue;
            }

            cells[cell_count][0] = faces[i][0];
            cells[cell_count][1] = faces[i][1];
            cells[cell_count][2] = faces[i][2];
            cells[cell_count][3] = idx as i32;
            cell_count += 1;
        }
    }

    //end make_cell section
    //start make_map section

    // face2cell
    for i in 0..num_face {
        for j in 0..nnn {
            idx = nn_list[faces[i][0] as usize][j] as usize;
            if idx as i32 == faces[i][1] || idx as i32 == faces[i][2] {
                continue;
            }

            tmp = 0_f64;
            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[faces[i][1] as usize][k]).powi(2);
            }

            if tmp > min_dist2 {
                continue;
            };

            tmp = 0_f64;
            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[faces[i][2] as usize][k]).powi(2);
            }

            if tmp > min_dist2 {
                continue;
            }
            face2cell[i][0] = faces[i][0];
            face2cell[i][1] = faces[i][1];
            face2cell[i][2] = faces[i][2];
            face2cell[i][3] = idx as i32;
            break;
        }
    }

    // edge2cell
    for i in 0..num_edge {
        for j in 0..nnn {
            idx = nn_list[edges[i][0] as usize][j] as usize;
            if idx as i32 == edges[i][1] {
                continue;
            }
            tmp = 0_f64;

            for k in 0..4 {
                tmp += (vertices[idx][k] - vertices[edges[i][1] as usize][k]).powi(2);
            }

            if tmp > min_dist2 {
                continue;
            }

            edge2cell[i][0] = edges[i][0];
            edge2cell[i][1] = edges[i][1];
            edge2cell[i][2] = idx as i32;

            for k in j + 1..nnn {
                idx = nn_list[edges[i][0] as usize][k] as usize;
                if idx as i32 == edge2cell[i][1] {
                    continue;
                }

                tmp = 0_f64;
                for m in 0..4 {
                    tmp += (vertices[idx][m] - vertices[edge2cell[i][1] as usize][m]).powi(2);
                }
                if tmp > min_dist2 {
                    continue;
                };

                tmp = 0_f64;
                for m in 0..4 {
                    tmp += (vertices[idx][m] - vertices[edge2cell[i][2] as usize][m]).powi(2);
                }

                if tmp > min_dist2 {
                    continue;
                };

                edge2cell[i][3] = idx as i32;
                break;
            }
            break;
        }
    }

    // end edge2cell section

    //start quat_setup section
    let vertex_points: Vec<Qpoint> = {
        let mut visited_vert: [bool; num_vert] = [false; num_vert];
        let mut v_q: [f64; 4] = [0.0; 4];
        let mut v_c: [f64; 4] = [0.0; 4];

        let num_rot: i64 = match is_ico {
            true => {
                (10 * (5 * num_div * num_div * num_div + num_div) / (num_vert as i32 / 2)).into()
            }
            false => (10 * (5 * num_div * num_div * num_div + num_div)).into(),
        };

        let mut quat: Vec<Quat> = vec![Quat(0.0, 0.0, 0.0, 0.0, 0.0); num_rot.try_into().unwrap()];
        let mut vertex_points: Vec<Qpoint> = vec![
            Qpoint {
                ..Default::default()
            };
            num_vert
        ];

        for i in 0..num_cell {
            for j in 0..4 {
                if visited_vert[cells[i][j] as usize] == true {
                    continue;
                };
                visited_vert[cells[i][j] as usize] = true;

                for k in 0..4 {
                    for m in 0..2 {
                        vertex_points[cells[i][j] as usize].vec[k][m] =
                            vec_vertices[cells[i][j] as usize][k][m];
                    }
                }

                for k in 0..4 {
                    v_c[k] = 0_f64;
                    for m in 0..4 {
                        v_c[k] += vertices[cells[i][m] as usize][k];
                    }
                    v_q[k] = vertices[cells[i][j] as usize][k];
                }
                let w = f0 * weight(v_q, v_c);
                vertex_points[cells[i][j] as usize].weight = w;
            }
        }
        vertex_points
    };
    //end quat_setup section

    //start refine_edge section
    let (edge_points, num_edge_point): (Vec<Qpoint>, usize) = {
        let mut v_q: [f64; 4] = [0.0; 4];
        let mut v_c: [f64; 4] = [0.0; 4];
        let mut vec_d_v: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut edge_point_count = 0;

        let num_edge_point: usize = num_edge * (num_div as usize - 1);
        let mut edge_points: Vec<Qpoint> = vec![
            Qpoint {
                ..Default::default()
            };
            num_edge_point
        ];

        for i in 0..num_edge {
            for j in 0..4 {
                vec_d_v[j][0] = (vec_vertices[edges[i][1] as usize][j][0]
                    - vec_vertices[edges[i][0] as usize][j][0])
                    / num_div;
                vec_d_v[j][1] = (vec_vertices[edges[i][1] as usize][j][1]
                    - vec_vertices[edges[i][0] as usize][j][1])
                    / num_div;
            }
            for j in 0..4 {
                v_c[j] = 0_f64;
                for k in 0..4 {
                    v_c[j] += vertices[edge2cell[i][k] as usize][j];
                }
            }
            for j in 1..num_div {
                for k in 0..4 {
                    edge_points[edge_point_count].vec[k][0] =
                        vec_vertices[edges[i][0] as usize][k][0] + j * vec_d_v[k][0];
                    edge_points[edge_point_count].vec[k][1] =
                        vec_vertices[edges[i][0] as usize][k][1] + j * vec_d_v[k][1];
                    v_q[k] = (edge_points[edge_point_count].vec[k][0] as f64
                        + tau * edge_points[edge_point_count].vec[k][1] as f64)
                        / (2.0 * num_div as f64);
                }

                let w = f1 * weight(v_q, v_c);
                edge_points[edge_point_count].weight = w;
                edge_point_count += 1;
            }
        }
        (edge_points, num_edge_point)
    };
    //end refine_edge section

    //start refine_face section
    let (face_points, num_face_point): (Vec<Qpoint>, usize) = {
        let mut v_q: [f64; 4] = [0.0; 4];
        let mut v_c: [f64; 4] = [0.0; 4];
        let mut vec_d_v1: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut vec_d_v2: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut face_point_count = 0;

        let num_face_point = num_face * (num_div as usize - 2) * (num_div as usize - 1) / 2;
        let mut face_points: Vec<Qpoint> = vec![
            Qpoint {
                ..Default::default()
            };
            num_face_point
        ];

        for i in 0..num_face {
            for j in 0..4 {
                vec_d_v1[j][0] = (vec_vertices[faces[i][1] as usize][j][0]
                    - vec_vertices[faces[i][0] as usize][j][0])
                    / num_div;
                vec_d_v1[j][1] = (vec_vertices[faces[i][1] as usize][j][1]
                    - vec_vertices[faces[i][0] as usize][j][1])
                    / num_div;
                vec_d_v2[j][0] = (vec_vertices[faces[i][2] as usize][j][0]
                    - vec_vertices[faces[i][0] as usize][j][0])
                    / num_div;
                vec_d_v2[j][1] = (vec_vertices[faces[i][2] as usize][j][1]
                    - vec_vertices[faces[i][0] as usize][j][1])
                    / num_div;
            }
            for j in 0..4 {
                v_c[j] = 0_f64;
                for k in 0..4 {
                    v_c[j] += vertices[face2cell[i][k] as usize][j];
                }
            }
            for j in 1..num_div - 1 {
                for k in 1..num_div - j {
                    for m in 0..4 {
                        face_points[face_point_count].vec[m][0] = vec_vertices
                            [faces[i][0] as usize][m][0]
                            + j * vec_d_v1[m][0]
                            + k * vec_d_v2[m][0];
                        face_points[face_point_count].vec[m][1] = vec_vertices
                            [faces[i][0] as usize][m][1]
                            + j * vec_d_v1[m][1]
                            + k * vec_d_v2[m][1];
                        v_q[m] = (face_points[face_point_count].vec[m][0] as f64
                            + tau * face_points[face_point_count].vec[m][1] as f64)
                            / (2.0 * num_div as f64);
                    }
                    let w: f64 = weight(v_q, v_c);
                    face_points[face_point_count].weight = w;
                    face_point_count += 1;
                }
            }
        }
        (face_points, num_face_point)
    };
    //end refine_face section

    //start refine_cell section
    let (cell_points, num_cell_point): (Vec<Qpoint>, usize) = {
        let mut v_q: [f64; 4] = [0.0; 4];
        let mut v_c: [f64; 4] = [0.0; 4];
        let mut vec_d_v1: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut vec_d_v2: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut vec_d_v3: [[i32; 2]; 4] = [[0; 2]; 4];
        let mut cell_point_count = 0;

        let num_cell_point: usize =
            ((num_cell as i32) / 6 * (num_div - 3) * (num_div - 2) * (num_div - 1)) as usize;
        let mut cell_points: Vec<Qpoint> = vec![
            Qpoint {
                ..Default::default()
            };
            num_cell_point
        ];
        for i in 0..num_cell {
            for j in 0..4 {
                vec_d_v1[j][0] = (vec_vertices[cells[i][1] as usize][j][0]
                    - vec_vertices[cells[i][0] as usize][j][0])
                    / num_div;
                vec_d_v1[j][1] = (vec_vertices[cells[i][1] as usize][j][1]
                    - vec_vertices[cells[i][0] as usize][j][1])
                    / num_div;
                vec_d_v2[j][0] = (vec_vertices[cells[i][2] as usize][j][0]
                    - vec_vertices[cells[i][0] as usize][j][0])
                    / num_div;
                vec_d_v2[j][1] = (vec_vertices[cells[i][2] as usize][j][1]
                    - vec_vertices[cells[i][0] as usize][j][1])
                    / num_div;
                vec_d_v3[j][0] = (vec_vertices[cells[i][3] as usize][j][0]
                    - vec_vertices[cells[i][0] as usize][j][0])
                    / num_div;
                vec_d_v3[j][1] = (vec_vertices[cells[i][3] as usize][j][1]
                    - vec_vertices[cells[i][0] as usize][j][1])
                    / num_div;
            }
            for j in 0..4 {
                v_c[j] = 0_f64;
                for k in 0..4 {
                    v_c[j] += vertices[cells[i][k] as usize][j];
                }
            }
            for j in 1..num_div - 2 {
                for k in 1..num_div - 1 - j {
                    for m in 1..num_div - j - k {
                        for n in 0..4 {
                            cell_points[cell_point_count].vec[n][0] = vec_vertices
                                [cells[i][0] as usize][n][0]
                                + j * vec_d_v1[n][0]
                                + k * vec_d_v2[n][0]
                                + m * vec_d_v3[n][0];
                            cell_points[cell_point_count].vec[n][1] = vec_vertices
                                [cells[i][0] as usize][n][1]
                                + j * vec_d_v1[n][1]
                                + k * vec_d_v2[n][1]
                                + m * vec_d_v3[n][1];
                            v_q[n] = (cell_points[cell_point_count].vec[n][0] as f64
                                + tau * cell_points[cell_point_count].vec[n][1] as f64)
                                / (2.0 * num_div as f64);
                        }
                        let w: f64 = weight(v_q, v_c);
                        cell_points[cell_point_count].weight = w;
                        cell_point_count += 1;
                    }
                }
            }
        }
        (cell_points, num_cell_point)
    };
    // end refine_cell section

    // start print_ico_quat section
    if is_ico {
        // continue
    }
    // start print_full_quat section
    else {
        let mut q_v: [f64; 4] = [0.0; 4];
        let num_rot: usize = (num_vert + num_edge_point + num_face_point + num_cell_point) / 2;

        if num_rot != 10 * (5 * num_div * num_div * num_div + num_div) as usize {
            return Err("wrong num_rot!".into());
        };
        println!("num_rot = {num_rot}");
        let mut fp: BufWriter<File> = match is_bin {
            true => {
                let filename = format!("c_quaternion{num_div}.bin");
                let path: PathBuf = [data_dir, filename].iter().collect();
                let f = File::create(path)?;
                BufWriter::new(f)
            }
            false => {
                let filename = format!("c_quaternion{num_div}.dat");
                let path: PathBuf = [data_dir, filename].iter().collect();
                let f = File::create(path)?;
                BufWriter::new(f)
            }
        };

        if is_bin {
            fp.write(&num_rot.to_le_bytes())?;
        }

        // select half of quaternions on vertices
        let mut ct: i32 = 0;
        for r in 0..num_vert {
            let flag = select_quat(vertex_points[r]);
            if flag != 1 {
                continue;
            }

            let mut q_norm = 0_f64;
            for i in 0..4 {
                q_v[i] = (vertex_points[r].vec[i][0] as f64
                    + tau * vertex_points[r].vec[i][1] as f64)
                    / (2.0 * num_div as f64);
                q_norm += q_v[i].powi(2);
            }
            q_norm = q_norm.sqrt();
            for i in 0..4 {
                q_v[i] /= q_norm;
            }

            if is_bin {
                fp.write(&q_v[0].to_le_bytes())?;
                fp.write(&q_v[1].to_le_bytes())?;
                fp.write(&q_v[2].to_le_bytes())?;
                fp.write(&q_v[3].to_le_bytes())?;
                fp.write(&vertex_points[r].weight.to_le_bytes())?;
            } else {
                let x = format!(
                    "{:.12} {:.12} {:.12} {:.12} {:.12}",
                    q_v[0], q_v[1], q_v[2], q_v[3], vertex_points[r].weight
                );
                writeln!(&mut fp, "{}", x)?;
            }
            ct += 1;
        }
        if ct != num_vert as i32 / 2 {
            return Err("Wrong number of quaternions on vertices!".into());
        }

        /* select half of the quaternions on edges */
        ct = 0;
        for r in 0..num_edge_point {
            let flag = select_quat(edge_points[r]);
            if flag != 1 {
                continue;
            }

            let mut q_norm: f64 = 0.0;
            for i in 0..4 {
                q_v[i] = (edge_points[r].vec[i][0] as f64 + tau * edge_points[r].vec[i][1] as f64)
                    / (2.0 * num_div as f64);
                q_norm += q_v[i].powi(2);
            }

            q_norm = q_norm.sqrt();
            for i in 0..4 {
                q_v[i] /= q_norm;
            }
            if is_bin {
                fp.write(&q_v[0].to_le_bytes())?;
                fp.write(&q_v[1].to_le_bytes())?;
                fp.write(&q_v[2].to_le_bytes())?;
                fp.write(&q_v[3].to_le_bytes())?;
                fp.write(&edge_points[r].weight.to_le_bytes())?;
            } else {
                let x = format!(
                    "{:.12} {:.12} {:.12} {:.12} {:.12}",
                    q_v[0], q_v[1], q_v[2], q_v[3], edge_points[r].weight
                );
                writeln!(&mut fp, "{}", x)?;
            }
            ct += 1;
        }
        if ct != num_edge_point as i32 / 2 {
            return Err("Wrong number of quaternions on edges!".into());
        }

        /* select half of the quaternions on faces */
        ct = 0;
        for r in 0..num_face_point {
            let flag = select_quat(face_points[r]);
            if flag != 1 {
                continue;
            }

            let mut q_norm: f64 = 0_f64;
            for i in 0..4 {
                q_v[i] = (face_points[r].vec[i][0] as f64 + tau * face_points[r].vec[i][1] as f64)
                    / (2.0 * num_div as f64);
                q_norm += q_v[i].powi(2);
            }
            q_norm = q_norm.sqrt();
            for i in 0..4 {
                q_v[i] /= q_norm
            }
            if is_bin {
                fp.write(&q_v[0].to_le_bytes())?;
                fp.write(&q_v[1].to_le_bytes())?;
                fp.write(&q_v[2].to_le_bytes())?;
                fp.write(&q_v[3].to_le_bytes())?;
                fp.write(&face_points[r].weight.to_le_bytes())?;
            } else {
                let x = format!(
                    "{:.12} {:.12} {:.12} {:.12} {:.12}",
                    q_v[0], q_v[1], q_v[2], q_v[3], face_points[r].weight
                );
                writeln!(&mut fp, "{}", x)?;
            }
            ct += 1;
        }

        if ct != num_face_point as i32 / 2 {
            return Err("Wrong number of quaternions on faces!".into());
        }
        ct = 0;
        for r in 0..num_cell_point {
            let flag = select_quat(cell_points[r]);
            if flag != 1 {
                continue;
            }

            let mut q_norm: f64 = 0_f64;
            for i in 0..4 {
                q_v[i] = (cell_points[r].vec[i][0] as f64 + tau * cell_points[r].vec[i][1] as f64)
                    / (2.0 * num_div as f64);
                q_norm += q_v[i].powi(2);
            }
            q_norm = q_norm.sqrt();
            for i in 0..4 {
                q_v[i] /= q_norm
            }
            if is_bin {
                fp.write(&q_v[0].to_le_bytes())?;
                fp.write(&q_v[1].to_le_bytes())?;
                fp.write(&q_v[2].to_le_bytes())?;
                fp.write(&q_v[3].to_le_bytes())?;
                fp.write(&cell_points[r].weight.to_le_bytes())?;
            } else {
                let x = format!(
                    "{:.12} {:.12} {:.12} {:.12} {:.12}",
                    q_v[0], q_v[1], q_v[2], q_v[3], cell_points[r].weight
                );
                writeln!(&mut fp, "{}", x)?;
            }
            ct += 1;
        }

        if ct != num_cell_point as i32 / 2 {
            return Err("Wrong number of quaternions on cells!".into());
        }
    }

    Ok(())
}

fn weight(v_q: [f64; 4], v_c: [f64; 4]) -> f64 {
    let mut w = 0.0;
    let mut norm_q = 0.0;
    let mut norm_c = 0.0;

    for i in 0..4 {
        norm_q += v_q[i].powi(2);
        norm_c += v_c[i].powi(2);
    }

    norm_q = (norm_q).sqrt();
    norm_c = (norm_c).sqrt();
    for i in 0..4 {
        w += v_q[i] * v_c[i];
    }

    w /= norm_q.powi(4) * norm_c;

    w
}

fn select_quat(q_p: Qpoint) -> i32 {
    let mut flag = 0;
    for i in 0..4 {
        for j in 0..2 {
            if q_p.vec[i][j] > 0 {
                flag = 1;
                break;
            } else if q_p.vec[i][j] < 0 {
                flag = -1;
                break;
            }
        }
        if flag != 0 {
            break;
        }
    }
    return flag;
}

fn compare_ico_quat(r_p: Qpoint, q_p: Qpoint) -> i32 {
    //return 1 if r_p precedes q_p, -1 if r_p follows q_p, 0 if r_p and q_p are the same
    let mut flag = 0;

    for i in 0..4 {
        for j in 0..2 {
            if (r_p.vec[i][j] < q_p.vec[i][j]) {
                flag = 1;
                break;
            } else if (r_p.vec[i][j] > q_p.vec[i][j]) {
                flag = -1;
                break;
            }
        }
        if flag != 0 {
            break;
        }
    }

    return flag;
}
