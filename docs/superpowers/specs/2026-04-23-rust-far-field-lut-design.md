# Rust Far-Field LUT Baker Design

## Goal

Add a Rust implementation of the far-field relativistic ray-deflection LUT baker under `prototype/rust`, using `rebound-rs` from `https://github.com/Cosmic-Stars-Team/rebound-rs` on the `bind/v4.6.0` branch.

## Existing Context

- The current implementation lives in `src/far_field_lut_baker.py`.
- The Python baker computes a 2D LUT by tracing one 2D photon trajectory per X column with REBOUND IAS15 and sampling angular deflection at multiple inverse-radius targets along that trajectory.
- The Python tests in `tests/test_far_field_lut_baker.py` already define the core behavior worth preserving: axis mappings, flat-space zero deflection, agreement with a reference integrator, deterministic LUT shape/dtype, sequential/parallel agreement, and EXR output.

## Approaches Considered

### 1. Direct port into one Rust binary crate

Pros:
- Fastest to build.
- Minimal file count.

Cons:
- Core math, REBOUND integration, CLI, and EXR writing all end up entangled.
- Harder to test without driving the CLI.

### 2. Workspace with one crate exposing both `lib` and `bin`

Pros:
- Keeps the workspace requirement simple while still separating reusable logic from the command-line surface.
- Easy to write unit tests for the baking logic and thin CLI wiring.
- Lets the binary stay close to the Python script behavior.

Cons:
- Slightly more setup than a single `main.rs`.

### 3. Workspace with multiple crates (`core`, `cli`, maybe `ffi`)

Pros:
- Maximum separation.

Cons:
- Unnecessary complexity for the current migration.
- More maintenance surface with no immediate payoff.

## Selected Design

Use approach 2.

Create `prototype/rust` as a Cargo workspace with one member crate, `crates/far-field-lut-baker`. That crate will expose:

- A library with the ported math helpers, REBOUND-backed ray tracing, LUT generation, and EXR writing.
- A binary that mirrors the Python CLI arguments and default output path.

## Architecture

### Workspace Layout

- `rust/Cargo.toml`: workspace definition and shared dependency versions.
- `rust/crates/far-field-lut-baker/Cargo.toml`: crate manifest with `rebound-rs` git dependency and runtime crates.
- `rust/crates/far-field-lut-baker/src/lib.rs`: public API, ported from the Python baker.
- `rust/crates/far-field-lut-baker/src/main.rs`: CLI wrapper.

### REBOUND Integration

The Rust port will use `rebound-rs` traits for:

- Creating a simulation and a test particle.
- Selecting the IAS15 integrator and configuring `epsilon` and `min_dt`.
- Stepping the simulation forward one integrator step at a time.
- Registering an `additional_forces` callback.

The remote `bind/v4.6.0` branch does not expose safe acceleration setters on `ParticleRef`, so the callback will mutate the raw `reb_particle` through `rebound::bind` inside a small, localized `unsafe` block. This keeps the dependency pinned to the requested branch while still matching the Python baker’s force model.

### Parallelism

The library will preserve the Python behavior:

- `workers = 1` runs sequentially.
- `workers = 0` auto-detects CPU count.
- `workers > 1` splits columns into contiguous chunks and uses a dedicated Rayon thread pool to compute them in parallel.

### EXR Output

The library will write scanline EXR files in Rust. For compatibility with the Python script:

- Default output writes a single luminance channel.
- `--rgb` writes the LUT into `R` and zeros into `G/B`.

## Error Handling

- Invalid user inputs become `bail!`/`ensure!` validation errors.
- Impossible photon initial states remain explicit errors.
- EXR I/O errors propagate with context.

## Testing Strategy

Rust tests will cover:

- Axis mappings and inverse X mapping.
- Flat-space and boundary-shell deflection behavior.
- Agreement with a simple reference RK4 integrator for one representative sample.
- LUT allocation shape and sequential/parallel agreement.
- EXR file emission for a tiny image.

## Scope Boundaries

This migration only adds the Rust workspace and baker implementation. It does not remove the Python baker or rewire the existing Python packaging in this pass.
