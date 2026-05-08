use std::f64::consts::PI;

use anyhow::{Context, Result, ensure};
use rebound::{
    create_particle,
    simulation::{
        Integrator, Simulation, SimulationIntegratorWrite, SimulationParticlesRead,
        SimulationParticlesWrite, SimulationSettingsWrite, SimulationStateRead,
        SimulationTransferWrite,
    },
    types::Vec3d,
};

const G: f64 = 4.0 * PI * PI;

#[derive(Debug, Clone, Copy)]
struct Sample {
    time: f64,
    earth_moon_distance: f64,
    emb_sun_distance: f64,
}

#[derive(Debug, Clone, Copy)]
struct TidalKick {
    a_earth: Vec3d,
    a_moon: Vec3d,
    a_emb: Vec3d,
}

fn v_add(a: Vec3d, b: Vec3d) -> Vec3d {
    Vec3d(a.0 + b.0, a.1 + b.1, a.2 + b.2)
}

fn v_sub(a: Vec3d, b: Vec3d) -> Vec3d {
    Vec3d(a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

fn v_scale(v: Vec3d, s: f64) -> Vec3d {
    Vec3d(v.0 * s, v.1 * s, v.2 * s)
}

fn v_norm(v: Vec3d) -> f64 {
    (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
}

fn pos(sim: &Simulation, index: usize) -> Result<Vec3d> {
    sim.get_particle(index)
        .context("missing particle")?
        .position()
        .context("particle position unavailable")
}

fn vel(sim: &Simulation, index: usize) -> Result<Vec3d> {
    sim.get_particle(index)
        .context("missing particle")?
        .velocity()
        .context("particle velocity unavailable")
}

fn mass(sim: &Simulation, index: usize) -> Result<f64> {
    sim.get_particle(index)
        .context("missing particle")?
        .mass()
        .context("particle mass unavailable")
}

fn set_vel(sim: &mut Simulation, index: usize, velocity: Vec3d) -> Result<()> {
    let mut particle = sim.get_particle(index).context("missing particle")?;
    particle
        .set_velocity(velocity.0, velocity.1, velocity.2)
        .context("failed to set particle velocity")
}

fn barycenter(sim: &Simulation) -> Result<Vec3d> {
    ensure!(sim.n() > 0, "simulation has no particles");

    let mut total_mass = 0.0;
    let mut weighted = Vec3d(0.0, 0.0, 0.0);

    for particle in sim.particles() {
        let m = particle.mass().context("particle mass unavailable")?;
        let p = particle.position().context("particle position unavailable")?;
        total_mass += m;
        weighted = v_add(weighted, v_scale(p, m));
    }

    ensure!(total_mass > 0.0, "total mass must be positive");
    Ok(v_scale(weighted, 1.0 / total_mass))
}

fn sun_gravity(pos: Vec3d, sun_pos: Vec3d, sun_mass: f64) -> Vec3d {
    let r = v_sub(pos, sun_pos);
    let inv_r3 = 1.0 / v_norm(r).powi(3);
    v_scale(r, -G * sun_mass * inv_r3)
}

fn tidal_force(sub_sim: &Simulation, main_sim: &Simulation) -> Result<TidalKick> {
    let sun_pos = pos(main_sim, 0)?;
    let emb_pos = pos(main_sim, 1)?;
    let earth_pos = pos(sub_sim, 0)?;
    let moon_pos = pos(sub_sim, 1)?;
    let sub_bc = barycenter(sub_sim)?;

    let a_earth = v_sub(
        sun_gravity(v_add(emb_pos, earth_pos), sun_pos, 1.0),
        sun_gravity(v_add(emb_pos, sub_bc), sun_pos, 1.0),
    );
    let a_moon = v_sub(
        sun_gravity(v_add(emb_pos, moon_pos), sun_pos, 1.0),
        sun_gravity(v_add(emb_pos, sub_bc), sun_pos, 1.0),
    );

    let m_earth = mass(sub_sim, 0)?;
    let m_moon = mass(sub_sim, 1)?;
    let total_sub_mass = m_earth + m_moon;
    let backreaction = v_add(v_scale(a_earth, m_earth), v_scale(a_moon, m_moon));
    let a_emb = v_scale(backreaction, -1.0 / total_sub_mass);

    Ok(TidalKick {
        a_earth,
        a_moon,
        a_emb,
    })
}

fn apply_cross_kick(sub_sim: &mut Simulation, main_sim: &mut Simulation, dt_half: f64) -> Result<()> {
    let kick = tidal_force(sub_sim, main_sim)?;

    let earth_vel = vel(sub_sim, 0)?;
    let moon_vel = vel(sub_sim, 1)?;
    let emb_vel = vel(main_sim, 1)?;

    set_vel(sub_sim, 0, v_add(earth_vel, v_scale(kick.a_earth, dt_half)))?;
    set_vel(sub_sim, 1, v_add(moon_vel, v_scale(kick.a_moon, dt_half)))?;
    set_vel(main_sim, 1, v_add(emb_vel, v_scale(kick.a_emb, dt_half)))?;

    // WHFast caches internal coordinates; synchronize after external velocity edits.
    sub_sim.synchronize();
    main_sim.synchronize();
    Ok(())
}

fn make_main_sim() -> Result<Simulation> {
    let mut sim = Simulation::new();
    sim.set_g(G)
        .set_integrator(Integrator::Whfast)
        .ri_whfast()
        .set_safe_mode(1);

    let m_sun = 1.0;
    let m_emb = 3.0e-6;
    let a = 1.0;
    let v = (G * (m_sun + m_emb) / a).sqrt();

    sim.add_particle(create_particle! {
        mass: m_sun,
        x: 0.0,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: 0.0,
        vz: 0.0,
    })?
    .add_particle(create_particle! {
        mass: m_emb,
        x: a,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: v,
        vz: 0.0,
    })?
    .move_to_com();

    Ok(sim)
}

fn make_sub_sim() -> Result<Simulation> {
    let mut sim = Simulation::new();
    sim.set_g(G)
        .set_integrator(Integrator::Whfast)
        .ri_whfast()
        .set_safe_mode(1);

    let m_earth = 3.0e-6 * 0.987;
    let m_moon = 3.0e-6 * 0.013;
    let separation: f64 = 0.00257;
    let omega = (G * (m_earth + m_moon) / separation.powi(3)).sqrt();

    let x_earth = -m_moon / (m_earth + m_moon) * separation;
    let x_moon = m_earth / (m_earth + m_moon) * separation;
    let vy_earth = -omega * x_earth;
    let vy_moon = omega * x_moon;

    sim.add_particle(create_particle! {
        mass: m_earth,
        x: x_earth,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: vy_earth,
        vz: 0.0,
    })?
    .add_particle(create_particle! {
        mass: m_moon,
        x: x_moon,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: vy_moon,
        vz: 0.0,
    })?
    .move_to_com();

    Ok(sim)
}

fn make_sims(sub_ratio: usize) -> Result<(Simulation, Simulation)> {
    ensure!(sub_ratio > 0, "sub_ratio must be positive");
    let mut main_sim = make_main_sim()?;
    let mut sub_sim = make_sub_sim()?;

    main_sim.set_dt(1.0 / 365.0)?;
    sub_sim.set_dt((1.0 / 365.0) / sub_ratio as f64)?;

    Ok((main_sim, sub_sim))
}

fn bridge_step(main_sim: &mut Simulation, sub_sim: &mut Simulation, dt_outer: f64) -> Result<()> {
    apply_cross_kick(sub_sim, main_sim, 0.5 * dt_outer)?;

    let target_time = main_sim.t() + dt_outer;
    sub_sim.integrate(target_time)?;
    main_sim.integrate(target_time)?;

    apply_cross_kick(sub_sim, main_sim, 0.5 * dt_outer)?;
    Ok(())
}

fn integrate_bridge(
    dt_outer: f64,
    t_end: f64,
    n_samples: usize,
    sub_ratio: usize,
) -> Result<Vec<Sample>> {
    ensure!(n_samples > 0, "n_samples must be positive");

    let (mut main_sim, mut sub_sim) = make_sims(sub_ratio)?;
    main_sim.set_dt(dt_outer)?;
    sub_sim.set_dt(dt_outer / sub_ratio as f64)?;

    let mut samples = Vec::with_capacity(n_samples);

    for sample_index in 1..=n_samples {
        let target_t = t_end * sample_index as f64 / n_samples as f64;

        while main_sim.t() < target_t - 1.0e-14 {
            let step_dt = (target_t - main_sim.t()).min(dt_outer);
            bridge_step(&mut main_sim, &mut sub_sim, step_dt)?;
        }

        let earth = pos(&sub_sim, 0)?;
        let moon = pos(&sub_sim, 1)?;
        let emb = pos(&main_sim, 1)?;

        samples.push(Sample {
            time: target_t,
            earth_moon_distance: v_norm(v_sub(earth, moon)),
            emb_sun_distance: v_norm(emb),
        });
    }

    Ok(samples)
}

fn main() -> Result<()> {
    let samples = integrate_bridge(1.0 / 365.0, 1.0, 50, 50)?;

    for sample in samples.iter().take(5) {
        println!(
            "t={:.6}  d_em={:.6e} AU  d_emb={:.6e} AU",
            sample.time, sample.earth_moon_distance, sample.emb_sun_distance
        );
    }

    Ok(())
}
