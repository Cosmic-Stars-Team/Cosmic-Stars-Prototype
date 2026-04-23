# Rust Far-Field LUT Baker Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Rust workspace under `prototype/rust` containing a REBOUND-backed far-field LUT baker and CLI equivalent to the existing Python script.

**Architecture:** Create a single workspace member crate that exposes both library and binary targets. Port the Python baker’s mapping, ray tracing, parallel chunking, and EXR output logic into the library, then keep the CLI as a thin wrapper over that API.

**Tech Stack:** Rust, Cargo workspace, `rebound-rs` git dependency, `clap`, `rayon`, `exr`, `anyhow`

---

### Task 1: Scaffold the Rust workspace

**Files:**
- Create: `rust/Cargo.toml`
- Create: `rust/crates/far-field-lut-baker/Cargo.toml`
- Create: `rust/crates/far-field-lut-baker/src/lib.rs`
- Create: `rust/crates/far-field-lut-baker/src/main.rs`

- [ ] Define the workspace root and shared dependency versions.
- [ ] Add the crate manifest with a git dependency on `https://github.com/Cosmic-Stars-Team/rebound-rs` using branch `bind/v4.6.0`.
- [ ] Add a minimal library and binary stub so `cargo test` can compile the workspace.

### Task 2: Port the baker core

**Files:**
- Modify: `rust/crates/far-field-lut-baker/src/lib.rs`

- [ ] Port the scalar helpers: `critical_impact_parameter`, Y mapping, X mapping, inverse X mapping, angle unwrap, and initial ray state construction.
- [ ] Port the REBOUND-backed ray simulation setup and additional-forces callback.
- [ ] Port the deflection tracer, chunk builder, worker resolution, and LUT generation.
- [ ] Port EXR writing for grayscale and RGB output.

### Task 3: Add CLI parity

**Files:**
- Modify: `rust/crates/far-field-lut-baker/src/main.rs`

- [ ] Add CLI parsing matching the Python baker’s arguments and defaults.
- [ ] Call the library API and print the output path after writing the EXR.

### Task 4: Add regression coverage

**Files:**
- Modify: `rust/crates/far-field-lut-baker/src/lib.rs`

- [ ] Add unit tests for the axis mappings and inverse mapping.
- [ ] Add unit tests for flat-space deflection and boundary behavior.
- [ ] Add a small regression test comparing one sample against an RK4 reference integrator.
- [ ] Add tests for LUT allocation, parallel equivalence, and EXR file output.

### Task 5: Verify and ship

**Files:**
- Modify: repository git state

- [ ] Run `cargo test --manifest-path rust/Cargo.toml`.
- [ ] Run a tiny bake command through the CLI.
- [ ] Review the diff, commit the changes, and push to `origin/main`.
