use anyhow::Result;
use symplectic_bridge::SymplecticBridge;

fn main() -> Result<()> {
    let mut bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50)?;

    for _ in 0..5 {
        let snapshot = bridge.advance(1.0 / 50.0)?;
        let sun = snapshot.main_bodies[0];
        let earth = snapshot.subsystems[0].world_bodies[0];
        let moon = snapshot.subsystems[0].world_bodies[1];
        let emb = snapshot.main_bodies[1];
        let d_em = distance(earth.position, moon.position);
        let d_emb = distance(emb.position, sun.position);
        println!(
            "t={:.6}  d_em={:.6e} AU  d_emb={:.6e} AU",
            snapshot.time, d_em, d_emb
        );
    }

    Ok(())
}

fn distance(a: rebound::types::Vec3d, b: rebound::types::Vec3d) -> f64 {
    (a - b).length_squared().sqrt()
}
