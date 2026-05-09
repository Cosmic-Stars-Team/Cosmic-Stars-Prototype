use std::f64::consts::PI;

use anyhow::{Context, Result, ensure};
use rebound::{
    create_particle,
    simulation::{
        self, Integrator, Simulation, SimulationIntegratorWrite, SimulationParticlesRead,
        SimulationParticlesWrite, SimulationSettingsWrite, SimulationStateRead,
        SimulationTransferWrite,
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

fn point_mass_gravity(target_pos: Vec3d, source_pos: Vec3d, source_mass: f64) -> Vec3d {
    let r = target_pos - source_pos;
    let inv_r3 = 1.0 / r.length_squared().sqrt().powi(3);
    r * (-G * source_mass * inv_r3)
}

impl SymplecticBridge {
    fn body_state(sim: &Simulation, index: usize) -> Result<BodyState> {
        let particle = sim.get_particle(index).context("missing particle")?;
        Ok(BodyState {
            position: particle
                .position()
                .context("particle position unavailable")?,
            velocity: particle
                .velocity()
                .context("particle velocity unavailable")?,
            mass: particle.mass().context("particle mass unavailable")?,
        })
    }

    fn add_velocity_kick(sim: &mut Simulation, index: usize, dv: Vec3d) -> Result<()> {
        let mut particle = sim.get_particle(index).context("missing particle")?;
        let velocity = particle
            .velocity()
            .context("particle velocity unavailable")?;
        particle
            .set_velocity_vec3d(velocity + dv)
            .context("failed to set particle velocity")?;
        Ok(())
    }

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

        let mut main_sim = Simulation::new();
        main_sim.set_g(G).set_integrator(Integrator::Whfast);

        simulation::set_integrator_config!(main_sim, {
            safe_mode: 1,
        })?;

        let m_sun = 1.0;
        let m_emb = 3.0e-6;
        let a = 1.0;
        let v = (G * (m_sun + m_emb) / a).sqrt();

        main_sim
            .add_particle(create_particle! {
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

        let mut sub_sim = Simulation::new();
        sub_sim.set_g(G).set_integrator(Integrator::Whfast);

        simulation::set_integrator_config!(sub_sim, {
            safe_mode: 1,
        })?;

        let m_earth = 3.0e-6 * 0.987;
        let m_moon = 3.0e-6 * 0.013;
        let separation: f64 = 0.00257;
        let omega = (G * (m_earth + m_moon) / separation.powi(3)).sqrt();

        let x_earth = -m_moon / (m_earth + m_moon) * separation;
        let x_moon = m_earth / (m_earth + m_moon) * separation;
        let vy_earth = omega * x_earth;
        let vy_moon = omega * x_moon;

        sub_sim
            .add_particle(create_particle! {
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

        let subsystem = BridgeSubsystem::new(sub_sim, 1, 1, vec![0], dt_inner)?;

        Self::new(main_sim, vec![subsystem], dt_outer)
    }

    pub fn subsystem_kick(&self, subsystem_index: usize) -> Result<BridgeKick> {
        let subsystem = self
            .subsystems
            .get(subsystem_index)
            .context("subsystem index out of bounds")?;

        let host = Self::body_state(&self.main_sim, subsystem.host_main_index)?;
        let local_bodies = (0..subsystem.sim.n())
            .map(|index| Self::body_state(&subsystem.sim, index))
            .collect::<Result<Vec<_>>>()?;
        let perturbers = subsystem
            .perturber_main_indices
            .iter()
            .map(|&index| Self::body_state(&self.main_sim, index))
            .collect::<Result<Vec<_>>>()?;

        ensure!(!local_bodies.is_empty(), "simulation has no particles");
        let total_mass = local_bodies.iter().map(|body| body.mass).sum::<f64>();
        ensure!(total_mass > 0.0, "total mass must be positive");

        let point_host_acceleration =
            perturbers
                .iter()
                .fold(Vec3d(0.0, 0.0, 0.0), |acceleration, source| {
                    acceleration + point_mass_gravity(host.position, source.position, source.mass)
                });

        let mut direct_accelerations = Vec::with_capacity(local_bodies.len());

        for local_body in &local_bodies {
            let world_position = host.position + local_body.position;
            let body_acceleration =
                perturbers
                    .iter()
                    .fold(Vec3d(0.0, 0.0, 0.0), |acceleration, source| {
                        acceleration
                            + point_mass_gravity(world_position, source.position, source.mass)
                    });

            direct_accelerations.push(body_acceleration);
        }

        let weighted_acceleration = direct_accelerations
            .iter()
            .zip(local_bodies.iter())
            .fold(Vec3d(0.0, 0.0, 0.0), |weighted, (acceleration, body)| {
                weighted + *acceleration * body.mass
            });
        let average_acceleration = weighted_acceleration / total_mass;
        let body_accelerations = direct_accelerations
            .into_iter()
            .map(|acceleration| acceleration - average_acceleration)
            .collect();
        let reaction_acceleration = average_acceleration - point_host_acceleration;

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
                Self::add_velocity_kick(&mut subsystem.sim, body_index, *acceleration * dt_half)?;
            }
            Self::add_velocity_kick(
                &mut self.main_sim,
                subsystem.reaction_main_index,
                kick.reaction_acceleration * dt_half,
            )?;
        }

        self.main_sim.synchronize();
        for subsystem in &mut self.subsystems {
            subsystem.sim.synchronize();
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
        let main_bodies = (0..self.main_sim.n())
            .map(|index| Self::body_state(&self.main_sim, index))
            .collect::<Result<Vec<_>>>()?;
        let mut subsystems = Vec::with_capacity(self.subsystems.len());

        for subsystem in &self.subsystems {
            let host = main_bodies
                .get(subsystem.host_main_index)
                .copied()
                .context("missing particle")?;
            let local_bodies = (0..subsystem.sim.n())
                .map(|index| Self::body_state(&subsystem.sim, index))
                .collect::<Result<Vec<_>>>()?;
            let world_bodies = local_bodies
                .iter()
                .map(|body| BodyState {
                    position: host.position + body.position,
                    velocity: host.velocity + body.velocity,
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

#[cfg(test)]
mod tests {
    use super::*;

    fn weighted_average(accelerations: &[Vec3d], masses: &[f64]) -> Vec3d {
        let total_mass = masses.iter().sum::<f64>();
        let weighted = accelerations
            .iter()
            .zip(masses.iter())
            .fold(Vec3d(0.0, 0.0, 0.0), |sum, (acceleration, mass)| {
                sum + *acceleration * *mass
            });
        weighted / total_mass
    }

    fn acceleration_from_bodies(target_pos: Vec3d, source_bodies: &[BodyState]) -> Vec3d {
        source_bodies
            .iter()
            .fold(Vec3d(0.0, 0.0, 0.0), |acceleration, source| {
                acceleration + point_mass_gravity(target_pos, source.position, source.mass)
            })
    }

    #[test]
    fn earth_moon_demo_keeps_relative_sun_emb_distance_near_one_au() {
        let bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        let snapshot = bridge.snapshot().unwrap();

        let sun = snapshot.main_bodies[0];
        let emb = snapshot.main_bodies[1];
        let relative_distance = (emb.position - sun.position).length_squared().sqrt();
        let origin_distance = emb.position.length_squared().sqrt();

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

        let host =
            SymplecticBridge::body_state(&bridge.main_sim, subsystem.host_main_index).unwrap();
        let perturbers = subsystem
            .perturber_main_indices
            .iter()
            .map(|&index| SymplecticBridge::body_state(&bridge.main_sim, index).unwrap())
            .collect::<Vec<_>>();
        let local_bodies = (0..subsystem.sim.n())
            .map(|index| SymplecticBridge::body_state(&subsystem.sim, index).unwrap())
            .collect::<Vec<_>>();
        let masses = local_bodies
            .iter()
            .map(|body| body.mass)
            .collect::<Vec<_>>();
        let direct_accelerations = local_bodies
            .iter()
            .map(|body| acceleration_from_bodies(host.position + body.position, &perturbers))
            .collect::<Vec<_>>();

        let point_acceleration = acceleration_from_bodies(host.position, &perturbers);
        let average_acceleration = weighted_average(&direct_accelerations, &masses);
        let local_average = weighted_average(&kick.body_accelerations, &masses);
        let expected_reaction = average_acceleration - point_acceleration;

        let local_average_norm = local_average.length_squared().sqrt();
        let reaction_error_norm = (kick.reaction_acceleration - expected_reaction)
            .length_squared()
            .sqrt();

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
        let earth = SymplecticBridge::body_state(&bridge.subsystems[0].sim, 0).unwrap();
        let moon = SymplecticBridge::body_state(&bridge.subsystems[0].sim, 1).unwrap();

        assert!(earth.position.0 < 0.0);
        assert!(moon.position.0 > 0.0);
        assert!(earth.velocity.1 < 0.0);
        assert!(moon.velocity.1 > 0.0);
    }

    #[test]
    fn bridge_requires_substeps_to_evenly_divide_outer_step() {
        let bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50).unwrap();
        let main_sim = bridge.main_sim;
        let sub_sim = bridge.subsystems.into_iter().next().unwrap().sim;
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
