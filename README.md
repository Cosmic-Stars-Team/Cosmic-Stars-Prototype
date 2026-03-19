# Cosmic Stars Prototype

Cosmic Stars Prototype is the **Python backend prototype** of the Cosmic Stars project.

This stage is focused on **numerical validation**, not final production performance.
The goals are:

- Validate N-body numerical accuracy
- Validate key physical effects (for example Mercury perihelion precession)
- Stabilize the backend output contract (position, velocity, simulation time) before Rust migration

## Scope

- Language: Python
- Core libraries: REBOUND + REBOUNDx
- Integrator policy: **Mercurius only**
- Input: `data/solar_system.json`
- Generated outputs: `data/gen/*`

## Project Layout

```text
.
├── src/
│   ├── main_simulation.py
│   ├── integrator_check.py
│   ├── database.py
│   └── extra_function.py
├── data/
│   ├── solar_system.json
│   └── gen/
├── Justfile
├── pyproject.toml
└── LICENSE
```

## Integrator Choice (Mercurius Only)

This prototype standardizes on **Mercurius** as the only integrator in normal workflow.

Why Mercurius:

- Hybrid strategy for better practical performance/accuracy balance
- Better handling for close-encounter style scenarios than pure symplectic-only setup
- Good fit for prototype-stage validation workloads

Notes:

- `main_simulation.py` defaults to Mercurius
- `integrator_check.py` defaults to Mercurius
- Team policy for this prototype: do not switch to other integrators in routine runs

## Technologies

- **REBOUND**: N-body integration engine
- **REBOUNDx**: extra physical effects (GR term enabled in this project)
- **Python stdlib** (`argparse`, `csv`, `json`): CLI, validation output, data streaming
- **Just**: reproducible local workflow (`venv`, install, run)

## Runtime Environment

Recommended:

- Windows 10/11 (PowerShell)
- Python 3.10+
- Optional: `just` command runner

Install `just` (optional):

```powershell
cargo install just
```

## Quick Start

### Option A: Just (recommended)

Run numerical check with default fallback arguments:

```powershell
just integrator-check
```

Run with custom arguments:

```powershell
just integrator-check "--years 10000 --steps 50000 --sample-bodies Mercury,Earth --output-csv data/gen/custom_check.csv"
```

Run the main simulation stream:

```powershell
just simulate
```

### Option B: Python directly

```powershell
$env:PYTHONPATH="src"
$env:PYTHONPYCACHEPREFIX="__pycache__"
.\.venv\Scripts\python.exe -m integrator_check
.\.venv\Scripts\python.exe -m main_simulation
```

## Units and Data Conventions

- Distance: AU
- Time: yr
- Velocity: AU/yr
- Mass: Msun

In `data/solar_system.json`, angular fields (`inc`, `Omega`, `omega`, `M`) are stored in **degrees** and converted to radians in code before passing into REBOUND.

## Windows Notes: Possible REBOUNDx Adjustments

On Windows, `reboundx` may require environment-specific adjustments depending on Python version and toolchain.
Possible actions:

1. Prefer binary wheels first (`pip install reboundx`)
2. If wheel is unavailable, build from source with a working C/C++ toolchain
3. If import fails, check `.pyd/.dll` dependency resolution (PATH/runtime libraries)
4. Use `.tmp_reboundx/` as local temporary build/debug workspace if needed

These are platform/toolchain concerns, not changes to project physics logic.

## Generated Files

- Simulation stream: `data/gen/simulation_stream.jsonl`
- Validation CSV: `data/gen/integrator_samples.csv` (overridable by CLI args)

## License

This project is released under **GNU GPL v3.0**. See `LICENSE`.
