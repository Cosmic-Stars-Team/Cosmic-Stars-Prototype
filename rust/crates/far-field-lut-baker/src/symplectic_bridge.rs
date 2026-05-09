use std::f64::consts::PI;

use anyhow::{Context, Result, ensure};
use rebound::{
    bind, create_particle,
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

#[derive(Debug, Clone)]
pub struct SubsystemSnapshot {
    pub host_main_index: usize,
    pub reaction_main_index: usize,
    pub local_bodies: Vec<BodyState>,
    pub world_bodies: Vec<BodyState>,
}

#[derive(Debug, Clone)]
pub struct BridgeSnapshot {
    pub time: f64,
    pub main_bodies: Vec<BodyState>,
    pub subsystems: Vec<SubsystemSnapshot>,
}

#[derive(Debug, Clone)]
pub struct BridgeKick {
    pub body_accelerations: Vec<Vec3d>,
    pub reaction_acceleration: Vec3d,
}

pub struct BridgeSubsystem {
    pub sim: Simulation,
    pub host_main_index: usize,
    pub reaction_main_index: usize,
    pub perturber_main_indices: Vec<usize>,
    pub dt_inner: f64,
}

pub struct SymplecticBridge {
    pub main_sim: Simulation,
    pub subsystems: Vec<BridgeSubsystem>,
    pub dt_outer: f64,
}

impl BridgeSubsystem {
    pub fn new(
        mut sim: Simulation,
        host_main_index: usize,
        reaction_main_index: usize,
        perturber_main_indices: Vec<usize>,
        dt_inner: f64,
    ) -> Result<Self> {
        ensure!(dt_inner > 0.0, "dt_inner must be positive");
        ensure!(
            !perturber_main_indices.is_empty(),
            "perturber_main_indices must not be empty"
        );
        sim.set_dt(dt_inner)?;
        Ok(Self {
            sim,
            host_main_index,
            reaction_main_index,
            perturber_main_indices,
            dt_inner,
        })
    }
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

fn sim_body_states(sim: &Simulation) -> Result<Vec<BodyState>> {
    (0..sim.n()).map(|index| body_state(sim, index)).collect()
}

fn point_mass_gravity(target_pos: Vec3d, source_pos: Vec3d, source_mass: f64) -> Vec3d {
    let r = v_sub(target_pos, source_pos);
    let inv_r3 = 1.0 / v_norm(r).powi(3);
    v_scale(r, -G * source_mass * inv_r3)
}

fn acceleration_from_main(
    main_sim: &Simulation,
    source_indices: &[usize],
    target_pos: Vec3d,
) -> Result<Vec3d> {
    let mut acceleration = Vec3d(0.0, 0.0, 0.0);
    for &source_index in source_indices {
        let source_pos = pos(main_sim, source_index)?;
        let source_mass = mass(main_sim, source_index)?;
        acceleration = v_add(
            acceleration,
            point_mass_gravity(target_pos, source_pos, source_mass),
        );
    }
    Ok(acceleration)
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
        mut main_sim: Simulation,
        mut subsystems: Vec<BridgeSubsystem>,
        dt_outer: f64,
    ) -> Result<Self> {
        ensure!(dt_outer > 0.0, "dt_outer must be positive");

        main_sim.set_dt(dt_outer)?;
        for subsystem in &mut subsystems {
            ensure!(
                subsystem.dt_inner <= dt_outer,
                "dt_inner must be <= dt_outer for every subsystem"
            );
            let substep_count = dt_outer / subsystem.dt_inner;
            let nearest_substep_count = substep_count.round();
            ensure!(
                (substep_count - nearest_substep_count).abs()
                    <= 1.0e-12 * substep_count.abs().max(1.0),
                "dt_outer must be an integer multiple of every dt_inner"
            );
            subsystem.sim.set_dt(subsystem.dt_inner)?;
        }

        Ok(Self {
            main_sim,
            subsystems,
            dt_outer,
        })
    }

    pub fn new_earth_moon(dt_outer: f64, sub_ratio: usize) -> Result<Self> {
        ensure!(sub_ratio > 0, "sub_ratio must be positive");
        let dt_inner = dt_outer / sub_ratio as f64;

        let main_sim = make_main_sim()?;
        let sub_sim = make_sub_sim()?;
        let subsystem = BridgeSubsystem::new(sub_sim, 1, 1, vec![0], dt_inner)?;

        Self::new(main_sim, vec![subsystem], dt_outer)
    }

    pub fn subsystem_kick(&self, subsystem_index: usize) -> Result<BridgeKick> {
        let subsystem = self
            .subsystems
            .get(subsystem_index)
            .context("subsystem index out of bounds")?;

        let host_position = pos(&self.main_sim, subsystem.host_main_index)?;
        let point_host_acceleration = acceleration_from_main(
            &self.main_sim,
            &subsystem.perturber_main_indices,
            host_position,
        )?;

        let mut total_mass = 0.0;
        let mut direct_accelerations = Vec::with_capacity(subsystem.sim.n());
        let mut masses = Vec::with_capacity(subsystem.sim.n());

        for body_index in 0..subsystem.sim.n() {
            let local_position = pos(&subsystem.sim, body_index)?;
            let body_mass = mass(&subsystem.sim, body_index)?;
            let world_position = v_add(host_position, local_position);
            let body_acceleration = acceleration_from_main(
                &self.main_sim,
                &subsystem.perturber_main_indices,
                world_position,
            )?;

            total_mass += body_mass;
            masses.push(body_mass);
            direct_accelerations.push(body_acceleration);
        }

        ensure!(total_mass > 0.0, "subsystem mass must be positive");
        let weighted_acceleration = direct_accelerations
            .iter()
            .zip(masses.iter())
            .fold(Vec3d(0.0, 0.0, 0.0), |sum, (acceleration, body_mass)| {
                v_add(sum, v_scale(*acceleration, *body_mass))
            });
        let average_acceleration = v_scale(weighted_acceleration, 1.0 / total_mass);
        let body_accelerations = direct_accelerations
            .into_iter()
            .map(|acceleration| v_sub(acceleration, average_acceleration))
            .collect();
        let reaction_acceleration = v_sub(average_acceleration, point_host_acceleration);

        Ok(BridgeKick {
            body_accelerations,
            reaction_acceleration,
        })
    }

    pub fn apply_cross_kick(&mut self, dt_half: f64) -> Result<()> {
        let kicks = (0..self.subsystems.len())
            .map(|index| self.subsystem_kick(index))
            .collect::<Result<Vec<_>>>()?;

        for (subsystem, kick) in self.subsystems.iter_mut().zip(kicks.iter()) {
            for (body_index, acceleration) in kick.body_accelerations.iter().enumerate() {
                add_velocity_kick(
                    &mut subsystem.sim,
                    body_index,
                    v_scale(*acceleration, dt_half),
                )?;
            }
            add_velocity_kick(
                &mut self.main_sim,
                subsystem.reaction_main_index,
                v_scale(kick.reaction_acceleration, dt_half),
            )?;
        }

        synchronize_after_external_edit(&mut self.main_sim);
        for subsystem in &mut self.subsystems {
            synchronize_after_external_edit(&mut subsystem.sim);
        }
        Ok(())
    }

    pub fn step(&mut self, dt: f64) -> Result<()> {
        for (index, subsystem) in self.subsystems.iter().enumerate() {
            ensure!(
                (subsystem.sim.t() - self.main_sim.t()).abs() < 1.0e-12,
                "main_sim and subsystem {index} times are out of sync"
            );
        }

        self.apply_cross_kick(0.5 * dt)?;

        let target_time = self.main_sim.t() + dt;
        self.main_sim.integrate(target_time)?;
        for subsystem in &mut self.subsystems {
            subsystem.sim.integrate(target_time)?;
        }

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
        let main_bodies = sim_body_states(&self.main_sim)?;
        let mut subsystems = Vec::with_capacity(self.subsystems.len());

        for subsystem in &self.subsystems {
            let host = body_state(&self.main_sim, subsystem.host_main_index)?;
            let local_bodies = sim_body_states(&subsystem.sim)?;
            let world_bodies = local_bodies
                .iter()
                .map(|body| BodyState {
                    position: v_add(host.position, body.position),
                    velocity: v_add(host.velocity, body.velocity),
                    mass: body.mass,
                })
                .collect();

            subsystems.push(SubsystemSnapshot {
                host_main_index: subsystem.host_main_index,
                reaction_main_index: subsystem.reaction_main_index,
                local_bodies,
                world_bodies,
            });
        }

        Ok(BridgeSnapshot {
            time: self.main_sim.t(),
            main_bodies,
            subsystems,
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
    let vy_earth = omega * x_earth;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn weighted_average(accelerations: &[Vec3d], masses: &[f64]) -> Vec3d {
        let total_mass = masses.iter().sum::<f64>();
        let weighted = accelerations
            .iter()
            .zip(masses.iter())
            .fold(Vec3d(0.0, 0.0, 0.0), |sum, (acceleration, mass)| {
                v_add(sum, v_scale(*acceleration, *mass))
            });
        v_scale(weighted, 1.0 / total_mass)
    }

    #[test]
    fn earth_moon_demo_keeps_relative_sun_emb_distance_near_one_au() {
        let bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        let snapshot = bridge.snapshot().unwrap();

        let sun = snapshot.main_bodies[0];
        let emb = snapshot.main_bodies[1];
        let relative_distance = v_norm(v_sub(emb.position, sun.position));
        let origin_distance = v_norm(emb.position);

        assert!((relative_distance - 1.0).abs() < 1.0e-9);
        assert!(
            (origin_distance - relative_distance).abs() > 1.0e-12,
            "move_to_com should shift the origin away from the sun"
        );
    }

    #[test]
    fn subsystem_kick_preserves_local_barycenter_and_applies_finite_size_reaction() {
        let bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        let subsystem = &bridge.subsystems[0];
        let kick = bridge.subsystem_kick(0).unwrap();

        let host_position = pos(&bridge.main_sim, subsystem.host_main_index).unwrap();
        let point_acceleration = acceleration_from_main(
            &bridge.main_sim,
            &subsystem.perturber_main_indices,
            host_position,
        )
        .unwrap();

        let masses = (0..subsystem.sim.n())
            .map(|index| mass(&subsystem.sim, index).unwrap())
            .collect::<Vec<_>>();
        let direct_accelerations = (0..subsystem.sim.n())
            .map(|index| {
                acceleration_from_main(
                    &bridge.main_sim,
                    &subsystem.perturber_main_indices,
                    v_add(host_position, pos(&subsystem.sim, index).unwrap()),
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let average_acceleration = weighted_average(&direct_accelerations, &masses);
        let local_average = weighted_average(&kick.body_accelerations, &masses);
        let expected_reaction = v_sub(average_acceleration, point_acceleration);

        let local_average_norm = v_norm(local_average);
        let reaction_error_norm = v_norm(v_sub(kick.reaction_acceleration, expected_reaction));

        assert!(
            local_average_norm < 1.0e-14,
            "local average acceleration norm={local_average_norm:e}"
        );
        assert!(
            reaction_error_norm < 1.0e-14,
            "reaction error norm={reaction_error_norm:e}"
        );
    }

    #[test]
    fn earth_moon_demo_initializes_counter_orbiting_velocities() {
        let bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        let earth = body_state(&bridge.subsystems[0].sim, 0).unwrap();
        let moon = body_state(&bridge.subsystems[0].sim, 1).unwrap();

        assert!(earth.position.0 < 0.0);
        assert!(moon.position.0 > 0.0);
        assert!(earth.velocity.1 < 0.0);
        assert!(moon.velocity.1 > 0.0);
    }

    #[test]
    fn bridge_requires_substeps_to_evenly_divide_outer_step() {
        let main_sim = make_main_sim().unwrap();
        let sub_sim = make_sub_sim().unwrap();
        let subsystem = BridgeSubsystem::new(sub_sim, 1, 1, vec![0], 0.004).unwrap();
        let result = SymplecticBridge::new(main_sim, vec![subsystem], 0.01);

        assert!(result.is_err());
    }

    #[test]
    fn step_rejects_out_of_sync_subsystem_time() {
        let mut bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        bridge.subsystems[0].sim.integrate(1.0 / 1000.0).unwrap();

        let result = bridge.step(1.0 / 365.0);

        assert!(result.is_err());
    }
}
