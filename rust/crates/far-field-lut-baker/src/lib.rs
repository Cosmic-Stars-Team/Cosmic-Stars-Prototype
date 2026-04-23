use std::{
    f64::consts::PI,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail, ensure};
use exr::prelude::{Image, SpecificChannels, Vec2, WritableImage, write_rgb_file};
use rayon::prelude::*;
use rebound::{
    create_particle,
    simulation::{
        Integrator, Simulation, SimulationCallbacksWrite, SimulationIntegratorWrite,
        SimulationParticlesRead, SimulationParticlesWrite, SimulationSettingsWrite,
        SimulationWrite,
    },
};

pub const DEFAULT_WIDTH: usize = 4096;
pub const DEFAULT_HEIGHT: usize = 4096;
pub const DEFAULT_RS: f64 = 1.0;
pub const DEFAULT_BOUNDARY_RADIUS_RS: f64 = 15.0;
pub const DEFAULT_CLUSTER_STRENGTH: f64 = 8.1;
pub const DEFAULT_EPSILON_SCALE: f64 = 1.0e-6;
pub const DEFAULT_IAS15_EPSILON: f64 = 1.0e-10;
pub const DEFAULT_IAS15_INITIAL_DT_SCALE: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct BakeConfig {
    pub width: usize,
    pub height: usize,
    pub rs: f64,
    pub boundary_radius_rs: f64,
    pub b_max: Option<f64>,
    pub cluster_strength: f64,
    pub epsilon_scale: f64,
    pub ias15_epsilon: f64,
    pub ias15_initial_dt_scale: f64,
    pub workers: isize,
}

impl Default for BakeConfig {
    fn default() -> Self {
        Self {
            width: DEFAULT_WIDTH,
            height: DEFAULT_HEIGHT,
            rs: DEFAULT_RS,
            boundary_radius_rs: DEFAULT_BOUNDARY_RADIUS_RS,
            b_max: None,
            cluster_strength: DEFAULT_CLUSTER_STRENGTH,
            epsilon_scale: DEFAULT_EPSILON_SCALE,
            ias15_epsilon: DEFAULT_IAS15_EPSILON,
            ias15_initial_dt_scale: DEFAULT_IAS15_INITIAL_DT_SCALE,
            workers: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LutBuffer {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f32>,
}

impl LutBuffer {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            data: vec![0.0; width * height],
        }
    }

    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.data[y * self.width + x]
    }

    pub fn set(&mut self, x: usize, y: usize, value: f32) {
        self.data[y * self.width + x] = value;
    }
}

pub fn critical_impact_parameter(rs: f64) -> f64 {
    1.5 * 3.0_f64.sqrt() * rs
}

pub fn map_pixel_y_to_u(y: usize, height: usize, rs: f64, boundary_radius_rs: f64) -> Result<f64> {
    ensure!(height >= 2, "height must be at least 2");
    ensure!(rs > 0.0, "rs must be > 0");
    ensure!(y < height, "y={y} is out of range for height={height}");

    let v = y as f64 / (height - 1) as f64;
    Ok(v / (boundary_radius_rs * rs))
}

pub fn map_pixel_x_to_b(
    x: usize,
    width: usize,
    b_crit: f64,
    b_max: f64,
    epsilon: f64,
    cluster_strength: f64,
) -> Result<f64> {
    ensure!(width >= 2, "width must be at least 2");
    ensure!(cluster_strength > 0.0, "cluster_strength must be > 0");
    ensure!(epsilon > 0.0, "epsilon must be > 0");
    ensure!(
        b_max > b_crit + epsilon,
        "b_max must be larger than b_crit + epsilon"
    );
    ensure!(x < width, "x={x} is out of range for width={width}");

    let t = x as f64 / (width - 1) as f64;
    let b_min = b_crit + epsilon;
    let span = b_max - b_min;
    let scaled = (cluster_strength * t).exp_m1() / cluster_strength.exp_m1();
    Ok(b_min + span * scaled)
}

pub fn map_b_to_uv_x(
    b: f64,
    b_crit: f64,
    b_max: f64,
    epsilon: f64,
    cluster_strength: f64,
) -> Result<f64> {
    ensure!(cluster_strength > 0.0, "cluster_strength must be > 0");
    ensure!(epsilon > 0.0, "epsilon must be > 0");
    let b_min = b_crit + epsilon;
    ensure!(b_max > b_min, "b_max must be larger than b_crit + epsilon");
    if !(b_min..=b_max).contains(&b) {
        bail!("b={b} is outside [{b_min}, {b_max}]");
    }

    let normalized = (b - b_min) / (b_max - b_min);
    Ok((normalized * cluster_strength.exp_m1()).ln_1p() / cluster_strength)
}

pub fn calculate_deflection(b: f64, u: f64, config: &BakeConfig) -> Result<f64> {
    ensure!(u >= 0.0, "u must be >= 0");
    ensure!(config.rs >= 0.0, "rs must be >= 0");
    if config.rs == 0.0 {
        return Ok(0.0);
    }

    let boundary_u = 1.0 / (config.boundary_radius_rs * config.rs);
    if u >= boundary_u {
        return Ok(0.0);
    }

    let values = trace_column_deflections(
        &[u, boundary_u],
        b,
        config.rs,
        config.boundary_radius_rs,
        config.ias15_epsilon,
        config.ias15_initial_dt_scale,
    )?;
    Ok(values[0])
}

pub fn generate_lut(config: &BakeConfig) -> Result<LutBuffer> {
    ensure!(
        config.width >= 2 && config.height >= 2,
        "width and height must be at least 2"
    );
    ensure!(config.rs > 0.0, "rs must be > 0");
    ensure!(config.epsilon_scale > 0.0, "epsilon_scale must be > 0");
    ensure!(config.ias15_epsilon > 0.0, "ias15_epsilon must be > 0");
    ensure!(
        config.ias15_initial_dt_scale > 0.0,
        "ias15_initial_dt_scale must be > 0"
    );

    let resolved_workers = resolve_worker_count(config.workers);
    ensure!(
        resolved_workers > 0,
        "workers must resolve to a positive integer"
    );

    let b_crit = critical_impact_parameter(config.rs);
    let epsilon = config.epsilon_scale * config.rs;
    let resolved_b_max = resolve_b_max(config);
    ensure!(
        resolved_b_max > b_crit + epsilon,
        "resolved b_max must exceed b_crit + epsilon"
    );

    let u_values = (0..config.height)
        .map(|y| map_pixel_y_to_u(y, config.height, config.rs, config.boundary_radius_rs))
        .collect::<Result<Vec<_>>>()?;

    let mut lut = LutBuffer::new(config.width, config.height);

    if resolved_workers == 1 {
        let (_, chunk) = generate_column_range(&ColumnRangeParams {
            x_start: 0,
            x_stop: config.width,
            width: config.width,
            u_values: &u_values,
            rs: config.rs,
            boundary_radius_rs: config.boundary_radius_rs,
            b_crit,
            b_max: resolved_b_max,
            epsilon,
            cluster_strength: config.cluster_strength,
            ias15_epsilon: config.ias15_epsilon,
            ias15_initial_dt_scale: config.ias15_initial_dt_scale,
        })?;
        lut.data.copy_from_slice(&chunk);
        return Ok(lut);
    }

    let chunk_ranges = build_column_chunks(config.width, resolved_workers);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(resolved_workers)
        .build()
        .context("failed to build rayon thread pool")?;

    let results = pool.install(|| {
        chunk_ranges
            .into_par_iter()
            .map(|(x_start, x_stop)| {
                generate_column_range(&ColumnRangeParams {
                    x_start,
                    x_stop,
                    width: config.width,
                    u_values: &u_values,
                    rs: config.rs,
                    boundary_radius_rs: config.boundary_radius_rs,
                    b_crit,
                    b_max: resolved_b_max,
                    epsilon,
                    cluster_strength: config.cluster_strength,
                    ias15_epsilon: config.ias15_epsilon,
                    ias15_initial_dt_scale: config.ias15_initial_dt_scale,
                })
            })
            .collect::<Vec<_>>()
    });

    for result in results {
        let (x_start, chunk) = result?;
        let chunk_width = chunk_width(config.width, x_start, &chunk, config.height)?;
        for y in 0..config.height {
            let src_start = y * chunk_width;
            let src_end = src_start + chunk_width;
            let dst_start = y * config.width + x_start;
            let dst_end = dst_start + chunk_width;
            lut.data[dst_start..dst_end].copy_from_slice(&chunk[src_start..src_end]);
        }
    }

    Ok(lut)
}

pub fn write_exr(output_path: impl AsRef<Path>, lut: &LutBuffer, rgb: bool) -> Result<PathBuf> {
    let output = output_path.as_ref().to_path_buf();
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    if rgb {
        write_rgb_file(&output, lut.width, lut.height, |x, y| {
            (lut.get(x, y), 0.0_f32, 0.0_f32)
        })
        .with_context(|| format!("failed to write RGB EXR {}", output.display()))?;
        return Ok(output);
    }

    let channels = SpecificChannels::build()
        .with_channel::<f32>("Y")
        .with_pixels(|Vec2(x, y)| (lut.get(x, y),));
    let image = Image::from_channels((lut.width, lut.height), channels);
    image
        .write()
        .to_file(&output)
        .with_context(|| format!("failed to write grayscale EXR {}", output.display()))?;

    Ok(output)
}

fn unwrap_angle(previous_angle: f64, mut current_angle: f64) -> f64 {
    let mut delta = current_angle - previous_angle;
    while delta <= -PI {
        current_angle += 2.0 * PI;
        delta = current_angle - previous_angle;
    }
    while delta > PI {
        current_angle -= 2.0 * PI;
        delta = current_angle - previous_angle;
    }
    current_angle
}

fn ray_state_from_b_u(b: f64, u: f64, rs: f64) -> Result<((f64, f64), (f64, f64))> {
    ensure!(b > 0.0, "b must be > 0");
    ensure!(u > 0.0, "u must be > 0");
    ensure!(rs >= 0.0, "rs must be >= 0");

    let radius = 1.0 / u;
    let v_tangent = b * u;
    let v_radial_sq = 1.0 - (b * b * u * u) + (rs * b * b * u * u * u);
    ensure!(
        v_radial_sq >= 0.0,
        "non-escaping configuration for the supplied (b, u, rs)"
    );

    let v_radial = -v_radial_sq.sqrt();
    Ok(((radius, 0.0), (v_radial, v_tangent)))
}

fn direction_angle(vx: f64, vy: f64) -> f64 {
    vy.atan2(vx)
}

fn make_ray_simulation(
    b: f64,
    start_u: f64,
    rs: f64,
    boundary_radius: f64,
    ias15_epsilon: f64,
    ias15_initial_dt_scale: f64,
) -> Result<Simulation> {
    let (position, velocity) = ray_state_from_b_u(b, start_u, rs)?;
    let l2 = b * b;

    let mut sim = Simulation::new();
    sim.set_integrator(Integrator::Ias15);
    sim.ri_ias15().set_epsilon(ias15_epsilon).set_min_dt(0.0);
    sim.set_dt(f64::max(
        1.0e-3 * rs.max(1.0),
        ias15_initial_dt_scale * (position.0 - boundary_radius),
    ))?;
    sim.add_particle(create_particle! {
        mass: 0.0,
        x: position.0,
        y: position.1,
        vx: velocity.0,
        vy: velocity.1,
    })?;

    sim.set_additional_forces(move |mut sim_ref| unsafe {
        let raw = sim_ref.raw_mut();
        let particle = &mut *(*raw).particles.add(0);
        let radius_sq = particle.x * particle.x + particle.y * particle.y + particle.z * particle.z;
        let radius = radius_sq.sqrt();
        let scale = -1.5 * rs * l2 / (radius_sq * radius_sq * radius);
        particle.ax += scale * particle.x;
        particle.ay += scale * particle.y;
        particle.az += scale * particle.z;
    });

    Ok(sim)
}

fn trace_column_deflections(
    u_targets: &[f64],
    b: f64,
    rs: f64,
    boundary_radius_rs: f64,
    ias15_epsilon: f64,
    ias15_initial_dt_scale: f64,
) -> Result<Vec<f64>> {
    if u_targets.is_empty() {
        return Ok(Vec::new());
    }
    ensure!(
        u_targets.iter().all(|u| *u >= 0.0),
        "u_targets must be >= 0"
    );
    ensure!(
        u_targets.windows(2).all(|pair| pair[0] <= pair[1]),
        "u_targets must be sorted in ascending order"
    );
    if rs == 0.0 {
        return Ok(vec![0.0; u_targets.len()]);
    }
    ensure!(ias15_epsilon > 0.0, "ias15_epsilon must be > 0");
    ensure!(
        ias15_initial_dt_scale > 0.0,
        "ias15_initial_dt_scale must be > 0"
    );

    let boundary_radius = boundary_radius_rs * rs;
    let boundary_u = 1.0 / boundary_radius;
    ensure!(
        u_targets.iter().all(|u| *u <= boundary_u),
        "u_targets must stay in the far-field domain outside the boundary sphere"
    );

    let mut out = vec![0.0; u_targets.len()];
    let positive_indices = u_targets
        .iter()
        .enumerate()
        .filter_map(|(index, u)| (*u > 0.0).then_some(index))
        .collect::<Vec<_>>();
    if positive_indices.is_empty() {
        return Ok(out);
    }

    let finite_u_targets = positive_indices
        .iter()
        .map(|index| u_targets[*index])
        .collect::<Vec<_>>();
    let start_u = finite_u_targets[0];
    let start_radius = 1.0 / start_u;
    if start_radius <= boundary_radius {
        return Ok(out);
    }

    let mut sim = make_ray_simulation(
        b,
        start_u,
        rs,
        boundary_radius,
        ias15_epsilon,
        ias15_initial_dt_scale,
    )?;

    let (mut previous_radius, mut previous_angle) = current_radius_and_angle(&sim)?;
    let mut finite_angles = vec![0.0; finite_u_targets.len()];
    finite_angles[0] = previous_angle;
    let mut next_target = 1;
    let boundary_angle = loop {
        sim.step();
        let (current_radius, raw_current_angle) = current_radius_and_angle(&sim)?;
        let current_angle = unwrap_angle(previous_angle, raw_current_angle);

        while next_target < finite_u_targets.len() {
            let target_radius = 1.0 / finite_u_targets[next_target];
            if !(current_radius <= target_radius && target_radius <= previous_radius) {
                break;
            }

            let blend = if previous_radius == current_radius {
                1.0
            } else {
                (previous_radius - target_radius) / (previous_radius - current_radius)
            };
            finite_angles[next_target] = previous_angle + blend * (current_angle - previous_angle);
            next_target += 1;
        }

        if current_radius <= boundary_radius {
            let blend = if previous_radius == current_radius {
                1.0
            } else {
                (previous_radius - boundary_radius) / (previous_radius - current_radius)
            };
            break previous_angle + blend * (current_angle - previous_angle);
        }

        previous_radius = current_radius;
        previous_angle = current_angle;
    };

    for (local_index, global_index) in positive_indices.iter().enumerate() {
        out[*global_index] = boundary_angle - finite_angles[local_index];
    }
    if u_targets[0] == 0.0 {
        out[0] = boundary_angle - PI;
    }

    Ok(out)
}

fn current_radius_and_angle(sim: &Simulation) -> Result<(f64, f64)> {
    let particle = sim
        .get_particle(0)
        .context("missing ray particle in simulation")?;
    let position = particle
        .position()
        .context("ray particle position not available")?;
    let velocity = particle
        .velocity()
        .context("ray particle velocity not available")?;
    Ok((
        (position.0 * position.0 + position.1 * position.1).sqrt(),
        direction_angle(velocity.0, velocity.1),
    ))
}

struct ColumnRangeParams<'a> {
    x_start: usize,
    x_stop: usize,
    width: usize,
    u_values: &'a [f64],
    rs: f64,
    boundary_radius_rs: f64,
    b_crit: f64,
    b_max: f64,
    epsilon: f64,
    cluster_strength: f64,
    ias15_epsilon: f64,
    ias15_initial_dt_scale: f64,
}

fn generate_column_range(params: &ColumnRangeParams<'_>) -> Result<(usize, Vec<f32>)> {
    let local_width = params.x_stop - params.x_start;
    let mut chunk = vec![0.0_f32; params.u_values.len() * local_width];

    for (local_x, x) in (params.x_start..params.x_stop).enumerate() {
        let b = map_pixel_x_to_b(
            x,
            params.width,
            params.b_crit,
            params.b_max,
            params.epsilon,
            params.cluster_strength,
        )?;
        let column = trace_column_deflections(
            params.u_values,
            b,
            params.rs,
            params.boundary_radius_rs,
            params.ias15_epsilon,
            params.ias15_initial_dt_scale,
        )?;
        for (y, value) in column.into_iter().enumerate() {
            chunk[y * local_width + local_x] = value as f32;
        }
    }

    Ok((params.x_start, chunk))
}

fn resolve_worker_count(workers: isize) -> usize {
    if workers == 1 {
        return 1;
    }

    if workers <= 0 {
        return std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
            .max(1);
    }

    workers as usize
}

fn build_column_chunks(width: usize, workers: usize) -> Vec<(usize, usize)> {
    let chunk_count = width.min(workers);
    let base = width / chunk_count;
    let remainder = width % chunk_count;
    let mut chunks = Vec::with_capacity(chunk_count);
    let mut start = 0;

    for index in 0..chunk_count {
        let stop = start + base + usize::from(index < remainder);
        chunks.push((start, stop));
        start = stop;
    }

    chunks
}

fn resolve_b_max(config: &BakeConfig) -> f64 {
    config
        .b_max
        .unwrap_or(config.boundary_radius_rs * config.rs)
}

fn chunk_width(total_width: usize, x_start: usize, chunk: &[f32], height: usize) -> Result<usize> {
    ensure!(height > 0, "height must be > 0");
    ensure!(
        chunk.len().is_multiple_of(height),
        "chunk length must be divisible by height"
    );
    let width = chunk.len() / height;
    ensure!(x_start + width <= total_width, "chunk exceeds lut width");
    Ok(width)
}

#[cfg(test)]
mod tests {
    use super::*;
    use exr::prelude::read_first_flat_layer_from_file;
    use tempfile::tempdir;

    fn initial_ray_state(b: f64, u: f64, rs: f64) -> Result<([f64; 2], [f64; 2])> {
        let r = 1.0 / u;
        let vt = b * u;
        let vr_sq = 1.0 - (b * b * u * u) + (rs * b * b * u * u * u);
        ensure!(
            vr_sq >= 0.0,
            "non-escaping configuration in reference integrator"
        );
        let vr = -vr_sq.sqrt();
        Ok(([r, 0.0], [vr, vt]))
    }

    fn schwarzschild_accel(position: [f64; 2], rs: f64, l2: f64) -> [f64; 2] {
        let radius_sq = position[0] * position[0] + position[1] * position[1];
        let radius = radius_sq.sqrt();
        let scale = -1.5 * rs * l2 / (radius_sq * radius_sq * radius);
        [position[0] * scale, position[1] * scale]
    }

    fn rk4_step(
        position: [f64; 2],
        velocity: [f64; 2],
        dt: f64,
        rs: f64,
        l2: f64,
    ) -> ([f64; 2], [f64; 2]) {
        fn deriv(pos: [f64; 2], vel: [f64; 2], rs: f64, l2: f64) -> ([f64; 2], [f64; 2]) {
            (vel, schwarzschild_accel(pos, rs, l2))
        }

        let (k1x, k1v) = deriv(position, velocity, rs, l2);
        let (k2x, k2v) = deriv(
            [
                position[0] + 0.5 * dt * k1x[0],
                position[1] + 0.5 * dt * k1x[1],
            ],
            [
                velocity[0] + 0.5 * dt * k1v[0],
                velocity[1] + 0.5 * dt * k1v[1],
            ],
            rs,
            l2,
        );
        let (k3x, k3v) = deriv(
            [
                position[0] + 0.5 * dt * k2x[0],
                position[1] + 0.5 * dt * k2x[1],
            ],
            [
                velocity[0] + 0.5 * dt * k2v[0],
                velocity[1] + 0.5 * dt * k2v[1],
            ],
            rs,
            l2,
        );
        let (k4x, k4v) = deriv(
            [position[0] + dt * k3x[0], position[1] + dt * k3x[1]],
            [velocity[0] + dt * k3v[0], velocity[1] + dt * k3v[1]],
            rs,
            l2,
        );

        let new_position = [
            position[0] + (dt / 6.0) * (k1x[0] + 2.0 * k2x[0] + 2.0 * k3x[0] + k4x[0]),
            position[1] + (dt / 6.0) * (k1x[1] + 2.0 * k2x[1] + 2.0 * k3x[1] + k4x[1]),
        ];
        let new_velocity = [
            velocity[0] + (dt / 6.0) * (k1v[0] + 2.0 * k2v[0] + 2.0 * k3v[0] + k4v[0]),
            velocity[1] + (dt / 6.0) * (k1v[1] + 2.0 * k2v[1] + 2.0 * k3v[1] + k4v[1]),
        ];

        (new_position, new_velocity)
    }

    fn reference_far_field_deflection(
        b: f64,
        u: f64,
        rs: f64,
        boundary_radius_rs: f64,
        dt: f64,
    ) -> Result<f64> {
        if rs == 0.0 {
            return Ok(0.0);
        }

        let boundary_radius = boundary_radius_rs * rs;
        let radius = 1.0 / u;
        if radius <= boundary_radius {
            return Ok(0.0);
        }

        let (mut position, mut velocity) = initial_ray_state(b, u, rs)?;
        let initial_theta = velocity[1].atan2(velocity[0]);
        let mut previous_position = position;
        let mut previous_velocity = velocity;

        while (position[0] * position[0] + position[1] * position[1]).sqrt() > boundary_radius {
            previous_position = position;
            previous_velocity = velocity;
            (position, velocity) = rk4_step(position, velocity, dt, rs, b * b);
        }

        let prev_radius = (previous_position[0] * previous_position[0]
            + previous_position[1] * previous_position[1])
            .sqrt();
        let curr_radius = (position[0] * position[0] + position[1] * position[1]).sqrt();
        let blend = (prev_radius - boundary_radius) / (prev_radius - curr_radius);
        let boundary_velocity = [
            previous_velocity[0] + blend * (velocity[0] - previous_velocity[0]),
            previous_velocity[1] + blend * (velocity[1] - previous_velocity[1]),
        ];
        let boundary_theta = boundary_velocity[1].atan2(boundary_velocity[0]);
        Ok(boundary_theta - initial_theta)
    }

    #[test]
    fn y_axis_maps_linearly_to_inverse_distance() {
        let rs = 2.0;
        let height = 4096;
        let max_u = 1.0 / (DEFAULT_BOUNDARY_RADIUS_RS * rs);

        assert_eq!(
            map_pixel_y_to_u(0, height, rs, DEFAULT_BOUNDARY_RADIUS_RS).unwrap(),
            0.0
        );
        assert!(
            (map_pixel_y_to_u(height - 1, height, rs, DEFAULT_BOUNDARY_RADIUS_RS).unwrap() - max_u)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn x_axis_maps_endpoints_correctly() {
        let rs = 1.0;
        let width = 4096;
        let epsilon = 1.0e-6;
        let b_crit = critical_impact_parameter(rs);
        let b_max = DEFAULT_BOUNDARY_RADIUS_RS * rs;

        assert!(
            (map_pixel_x_to_b(0, width, b_crit, b_max, epsilon, DEFAULT_CLUSTER_STRENGTH,)
                .unwrap()
                - (b_crit + epsilon))
                .abs()
                < 1.0e-12
        );
        assert!(
            (map_pixel_x_to_b(
                width - 1,
                width,
                b_crit,
                b_max,
                epsilon,
                DEFAULT_CLUSTER_STRENGTH,
            )
            .unwrap()
                - b_max)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn inverse_x_mapping_round_trips() {
        let rs = 1.0;
        let width = 4096;
        let epsilon = 1.0e-6;
        let b_crit = critical_impact_parameter(rs);
        let b_max = DEFAULT_BOUNDARY_RADIUS_RS * rs;

        for x in [
            0,
            17,
            width / 4,
            (0.8 * (width - 1) as f64).round() as usize,
            width - 1,
        ] {
            let b = map_pixel_x_to_b(x, width, b_crit, b_max, epsilon, DEFAULT_CLUSTER_STRENGTH)
                .unwrap();
            let uv_x = map_b_to_uv_x(b, b_crit, b_max, epsilon, DEFAULT_CLUSTER_STRENGTH).unwrap();
            let expected = x as f64 / (width - 1) as f64;
            assert!((uv_x - expected).abs() < 1.0e-9);
        }
    }

    #[test]
    fn calculate_deflection_is_zero_in_flat_space() {
        let config = BakeConfig {
            rs: 0.0,
            ..BakeConfig::default()
        };
        assert_eq!(calculate_deflection(6.0, 0.04, &config).unwrap(), 0.0);
    }

    #[test]
    fn calculate_deflection_matches_reference_integral() {
        let config = BakeConfig::default();
        let b = 6.5;
        let u = 1.0 / 20.0;
        let rs = 1.0;

        let expected =
            reference_far_field_deflection(b, u, rs, DEFAULT_BOUNDARY_RADIUS_RS, 1.0e-3).unwrap();
        let actual = calculate_deflection(b, u, &config).unwrap();

        assert!(
            (actual - expected).abs() <= 2.0e-4,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn calculate_deflection_is_zero_on_boundary_shell() {
        let config = BakeConfig::default();
        let u_boundary = 1.0 / (DEFAULT_BOUNDARY_RADIUS_RS * config.rs);
        assert_eq!(
            calculate_deflection(10.0, u_boundary, &config).unwrap(),
            0.0
        );
    }

    #[test]
    fn generate_lut_allocates_expected_float_buffer() {
        let config = BakeConfig {
            width: 8,
            height: 4,
            rs: 1.0,
            ..BakeConfig::default()
        };
        let lut = generate_lut(&config).unwrap();

        assert_eq!(lut.width, 8);
        assert_eq!(lut.height, 4);
        assert_eq!(lut.data.len(), 32);
        assert!(lut.data.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn generate_lut_parallel_matches_sequential() {
        let sequential = generate_lut(&BakeConfig {
            width: 8,
            height: 4,
            rs: 1.0,
            workers: 1,
            ..BakeConfig::default()
        })
        .unwrap();
        let parallel = generate_lut(&BakeConfig {
            width: 8,
            height: 4,
            rs: 1.0,
            workers: 2,
            ..BakeConfig::default()
        })
        .unwrap();

        for (left, right) in sequential.data.iter().zip(parallel.data.iter()) {
            assert!((left - right).abs() <= 1.0e-6);
        }
    }

    #[test]
    fn write_exr_outputs_single_channel_float_image() {
        let lut = generate_lut(&BakeConfig {
            width: 8,
            height: 4,
            rs: 1.0,
            ..BakeConfig::default()
        })
        .unwrap();

        let temp_dir = tempdir().unwrap();
        let output_path = temp_dir.path().join("lut.exr");
        write_exr(&output_path, &lut, false).unwrap();

        assert!(output_path.exists());
        let image = read_first_flat_layer_from_file(&output_path).unwrap();
        assert_eq!(image.layer_data.channel_data.list.len(), 1);
        assert_eq!(image.layer_data.channel_data.list[0].name.as_slice(), b"Y");
    }
}
