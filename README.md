# Cosmic Stars

Personal project: Solar system dynamics simulation implemented using the Rebound library, along with its general relativistic corrections.

## Features

### N-body Simulation
- High-precision N-body integration using REBOUND
- Support for IAS15 (adaptive) and WHFast (symplectic) integrators
- General relativistic corrections via REBOUNDx
- Flexible time-step control and user-adjustable speed factors (0.1x - 100x)

### Solar System Database
- **Complete Solar System Bodies**:
  - Sun and 8 major planets with accurate orbital elements
  - Major moons (including Galilean moons)
  - Dwarf planets (Pluto, Ceres, Eris, etc.)
- **Accurate Parameters**:
  - Semi-major axis, eccentricity, inclination
  - Mass, orbital periods
  - Verified against JPL horizons data

### Asteroid Belt System
- **Main Belt Asteroids**: 2.0-4.0 AU with configurable number (up to 20,000+)
- **Kirkwood Gaps**: Orbital resonance gaps with Jupiter
- **Hilda Group**: 3:2 resonance asteroids at ~3.97 AU
- **Trojans**: Jupiter's L4/L5 Trojan swarms
- **Realistic Distribution**:
  - Power-law size distribution
  - Proper orbital elements
  - Randomized anomalies

### Simulation Engines
- **simulation_engine.py**: Advanced engine with Unity time synchronization
- **sandbox_engine.py**: Optimized for smooth visual experience and real-time interaction
- **unity_interface.py**: HTTP API server for Unity frontend integration

## Installation

```bash
pip install rebound numpy astropy flask flask-cors
```

### Dependencies
- `rebound`: N-body simulation library
- `reboundx`: Additional physics effects (GR, post-Newtonian)
- `numpy`: Numerical computations
- `astropy`: Astronomical calculations
- `flask`: Web server for Unity integration
- `flask-cors`: Cross-origin resource sharing

## Project Structure

```
ASC/
├── main.py                  # Basic simulation examples and verification
├── celestial_bodies.py      # Solar System celestial body database
├── asteroid_belt.py         # Asteroid belt generation algorithms
├── simulation_engine.py     # Advanced simulation engine with Unity sync
├── sandbox_engine.py        # Sandbox engine for smooth visualization
└── unity_interface.py       # Unity HTTP API interface server
```

## Usage

### Basic Solar System Simulation

```python
from main import create_solar_system_simulation

# Create basic solar system
sim = create_solar_system_simulation()
print(f"Number of bodies: {len(sim.particles)}")
```

### Solar System with Moons

```python
from main import create_solar_system_with_moons

# Include major moons
sim = create_solar_system_with_moons()
print(f"Bodies with moons: {len(sim.particles)}")
```

### Realistic Asteroid Belt

```python
from main import create_realistic_asteroid_system

# Create asteroid system with 20,000 main belt asteroids
# plus 3,000 Hilda group and 5,000 Trojans
sim = create_realistic_asteroid_system(N=20000, seed=42)
```

### Unity Integration (HTTP API)

```python
# Start the Unity interface server
python unity_interface.py

# Server will be available at http://localhost:5000
# API endpoints:
# - POST /api/simulation/init  - Initialize simulation
# - POST /api/simulation/step  - Advance simulation
# - GET  /api/simulation/state - Get current state
# - POST /api/simulation/speed - Adjust time speed
```

### Run Verification Tests

```bash
python main.py
```

This will:
1. List all available celestial bodies in the database
2. Verify orbital elements against expected values
3. Create and demonstrate different simulation configurations

## Technical Details

### Unit System
- **Distance**: Astronomical Units (AU)
- **Time**: Years (yr)
- **Mass**: Solar masses (Msun)

### Integrators
- **IAS15**: Adaptive step-size, 15th order (default for high accuracy)
  - Best for close encounters and highly eccentric orbits
- **WHFast**: Symplectic, 2nd order (faster for large systems)
  - Best for large number of particles with stable orbits

### Time Control
- Adjustable time speed factor: 0.1x to 100x
- Synchronized with Unity at 60 FPS
- Configurable steps per frame for performance tuning
- Real-time orbital mechanics validation

### Performance
- Supports 20,000+ asteroids with WHFast integrator
- Optimized for smooth visualization (sandbox mode)
- Memory-efficient particle management

## Examples

### Custom Planetary System

```python
from main import create_custom_system
from celestial_bodies import SolarSystemBodies

# Create inner planets + Jupiter + Galilean moons
sim = create_custom_system()
```

### Verify Orbital Elements

```python
from main import verify_orbital_elements

# Check orbital element accuracy
verify_orbital_elements()
```

## Recent Updates

- Optimized asteroid belt generation algorithms
- Added Kirkwood gap, Hilda group, and Trojans simulation
- Implemented comprehensive Solar System planet database
- Added dwarf planets and major moons
- Created Unity integration interface with HTTP API
- Implemented time synchronization for real-time visualization

## License

Personal project for educational and research purposes.

## References

- REBOUND: https://github.com/hannorein/rebound
- REBOUNDx: https://github.com/dtamayo/reboundx
- JPL Horizons: https://ssd.jpl.nasa.gov/horizons.cgi
