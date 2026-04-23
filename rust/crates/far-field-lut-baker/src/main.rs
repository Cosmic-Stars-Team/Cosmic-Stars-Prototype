use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use far_field_lut_baker::{
    BakeConfig, DEFAULT_BOUNDARY_RADIUS_RS, DEFAULT_CLUSTER_STRENGTH, DEFAULT_EPSILON_SCALE,
    DEFAULT_HEIGHT, DEFAULT_IAS15_EPSILON, DEFAULT_IAS15_INITIAL_DT_SCALE, DEFAULT_RS,
    DEFAULT_WIDTH, generate_lut, write_exr,
};

#[derive(Debug, Parser)]
#[command(about = "Bake a 2D EXR LUT for far-field relativistic ray deflection.")]
struct Args {
    #[arg(long, default_value = "data/gen/far_field_ray_deflection_lut_4096.exr")]
    output: PathBuf,

    #[arg(long, default_value_t = DEFAULT_WIDTH)]
    width: usize,

    #[arg(long, default_value_t = DEFAULT_HEIGHT)]
    height: usize,

    #[arg(long, default_value_t = DEFAULT_RS)]
    rs: f64,

    #[arg(long = "boundary-radius-rs", default_value_t = DEFAULT_BOUNDARY_RADIUS_RS)]
    boundary_radius_rs: f64,

    #[arg(long = "b-max-rs", default_value_t = DEFAULT_BOUNDARY_RADIUS_RS)]
    b_max_rs: f64,

    #[arg(long = "cluster-strength", default_value_t = DEFAULT_CLUSTER_STRENGTH)]
    cluster_strength: f64,

    #[arg(long = "epsilon-scale", default_value_t = DEFAULT_EPSILON_SCALE)]
    epsilon_scale: f64,

    #[arg(long)]
    rgb: bool,

    #[arg(long = "ias15-epsilon", default_value_t = DEFAULT_IAS15_EPSILON)]
    ias15_epsilon: f64,

    #[arg(
        long = "ias15-initial-dt-scale",
        default_value_t = DEFAULT_IAS15_INITIAL_DT_SCALE
    )]
    ias15_initial_dt_scale: f64,

    #[arg(long, default_value_t = 0)]
    workers: isize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let config = BakeConfig {
        width: args.width,
        height: args.height,
        rs: args.rs,
        boundary_radius_rs: args.boundary_radius_rs,
        b_max: Some(args.b_max_rs * args.rs),
        cluster_strength: args.cluster_strength,
        epsilon_scale: args.epsilon_scale,
        ias15_epsilon: args.ias15_epsilon,
        ias15_initial_dt_scale: args.ias15_initial_dt_scale,
        workers: args.workers,
    };

    let lut = generate_lut(&config)?;
    let output = write_exr(&args.output, &lut, args.rgb)?;
    println!("Wrote LUT EXR to {}", output.display());

    Ok(())
}
