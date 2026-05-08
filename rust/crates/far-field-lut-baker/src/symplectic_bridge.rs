use std::f64::consts::PI;

use anyhow::{Context, Result, ensure};
use rebound::{
    bind,
    create_particle,
    simulation::{
        Integrator, Simulation, SimulationIntegratorWrite, SimulationParticlesRead,
        SimulationParticlesWrite, SimulationSettingsWrite, SimulationStateRead,
        SimulationTransferWrite, SimulationWrite,
    },
    types::Vec3d,
};

pub const G: f64 = 4.0 * PI * PI;

#[derive(Debug, Clone, Copy)]
pub struct BodyState {
    pub position: Vec3d,
    pub velocity: Vec3d,
    pub mass: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct BridgeSnapshot {
    pub time: f64,
    pub earth_moon_distance: f64,
    pub emb_sun_distance: f64,
    pub earth: BodyState,
    pub moon: BodyState,
    pub emb: BodyState,
}

#[derive(Debug, Clone, Copy)]
pub struct TidalKick {
    pub a_earth: Vec3d,
    pub a_moon: Vec3d,
    pub a_emb: Vec3d,
}

pub struct SymplecticBridge {
    pub main_sim: Simulation,
    pub sub_sim: Simulation,
    pub dt_outer: f64,
    pub dt_inner: f64,
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

fn body_state(sim: &Simulation, index: usize) -> Result<BodyState> {
    Ok(BodyState {
        position: pos(sim, index)?,
        velocity: vel(sim, index)?,
        mass: mass(sim, index)?,
    })
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

unsafe fn raw_particle_mut(sim: &mut Simulation, index: usize) -> Result<&mut bind::reb_particle> {
    if index >= sim.n() {
        anyhow::bail!("particle index {index} out of bounds");
    }
    let raw_sim = sim.raw_mut();
    Ok(unsafe { &mut *(*raw_sim).particles.add(index) })
}

fn add_velocity_kick(sim: &mut Simulation, index: usize, dv: Vec3d) -> Result<()> {
    let particle = unsafe { raw_particle_mut(sim, index)? };
    particle.vx += dv.0;
    particle.vy += dv.1;
    particle.vz += dv.2;
    Ok(())
}

fn synchronize_after_external_edit(sim: &mut Simulation) {
    unsafe {
        bind::reb_simulation_synchronize(sim.raw_mut());
    }
}

impl SymplecticBridge {
    pub fn new(
        main_sim: Simulation,
        sub_sim: Simulation,
        dt_outer: f64,
        dt_inner: f64,
    ) -> Result<Self> {
        ensure!(dt_outer > 0.0, "dt_outer must be positive");
        ensure!(dt_inner > 0.0, "dt_inner must be positive");
        ensure!(dt_inner <= dt_outer, "dt_inner must be <= dt_outer");

        let mut bridge = Self {
            main_sim,
            sub_sim,
            dt_outer,
            dt_inner,
        };

        bridge.main_sim.set_dt(dt_outer)?;
        bridge.sub_sim.set_dt(dt_inner)?;
        Ok(bridge)
    }

    pub fn new_earth_moon(dt_outer: f64, sub_ratio: usize) -> Result<Self> {
        ensure!(sub_ratio > 0, "sub_ratio must be positive");
        let dt_inner = dt_outer / sub_ratio as f64;
        let main_sim = make_main_sim()?;
        let sub_sim = make_sub_sim()?;
        Self::new(main_sim, sub_sim, dt_outer, dt_inner)
    }

    pub fn tidal_force(&self) -> Result<TidalKick> {
        let sun_pos = pos(&self.main_sim, 0)?;
        let emb_pos = pos(&self.main_sim, 1)?;
        let earth_pos = pos(&self.sub_sim, 0)?;
        let moon_pos = pos(&self.sub_sim, 1)?;
        let sub_bc = barycenter(&self.sub_sim)?;

        let a_earth = v_sub(
            sun_gravity(v_add(emb_pos, earth_pos), sun_pos, 1.0),
            sun_gravity(v_add(emb_pos, sub_bc), sun_pos, 1.0),
        );
        let a_moon = v_sub(
            sun_gravity(v_add(emb_pos, moon_pos), sun_pos, 1.0),
            sun_gravity(v_add(emb_pos, sub_bc), sun_pos, 1.0),
        );

        let m_earth = mass(&self.sub_sim, 0)?;
        let m_moon = mass(&self.sub_sim, 1)?;
        let total_sub_mass = m_earth + m_moon;
        let backreaction = v_add(v_scale(a_earth, m_earth), v_scale(a_moon, m_moon));
        let a_emb = v_scale(backreaction, -1.0 / total_sub_mass);

        Ok(TidalKick {
            a_earth,
            a_moon,
            a_emb,
        })
    }

    pub fn apply_cross_kick(&mut self, dt_half: f64) -> Result<()> {
        let kick = self.tidal_force()?;

        add_velocity_kick(&mut self.sub_sim, 0, v_scale(kick.a_earth, dt_half))?;
        add_velocity_kick(&mut self.sub_sim, 1, v_scale(kick.a_moon, dt_half))?;
        add_velocity_kick(&mut self.main_sim, 1, v_scale(kick.a_emb, dt_half))?;

        synchronize_after_external_edit(&mut self.sub_sim);
        synchronize_after_external_edit(&mut self.main_sim);
        Ok(())
    }

    pub fn step(&mut self, dt: f64) -> Result<()> {
        self.apply_cross_kick(0.5 * dt)?;

        let target_time = self.main_sim.t() + dt;
        self.sub_sim.integrate(target_time)?;
        self.main_sim.integrate(target_time)?;

        self.apply_cross_kick(0.5 * dt)?;
        Ok(())
    }

    pub fn advance(&mut self, duration: f64) -> Result<BridgeSnapshot> {
        ensure!(duration >= 0.0, "duration must be non-negative");

        let target_t = self.main_sim.t() + duration;
        while self.main_sim.t() < target_t - 1.0e-14 {
            let step_dt = (target_t - self.main_sim.t()).min(self.dt_outer);
            self.step(step_dt)?;
        }

        self.snapshot()
    }

    pub fn snapshot(&self) -> Result<BridgeSnapshot> {
        let earth = body_state(&self.sub_sim, 0)?;
        let moon = body_state(&self.sub_sim, 1)?;
        let emb = body_state(&self.main_sim, 1)?;

        Ok(BridgeSnapshot {
            time: self.main_sim.t(),
            earth_moon_distance: v_norm(v_sub(earth.position, moon.position)),
            emb_sun_distance: v_norm(emb.position),
            earth,
            moon,
            emb,
        })
    }
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
