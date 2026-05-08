use anyhow::Result;
use far_field_lut_baker::symplectic_bridge::SymplecticBridge;

fn main() -> Result<()> {
    let mut bridge = SymplecticBridge::new_earth_moon(1.0 / 365.0, 50)?;

    for _ in 0..5 {
        let snapshot = bridge.advance(1.0 / 50.0)?;
        println!(
            "t={:.6}  d_em={:.6e} AU  d_emb={:.6e} AU",
            snapshot.time, snapshot.earth_moon_distance, snapshot.emb_sun_distance
        );
    }

    Ok(())
}
