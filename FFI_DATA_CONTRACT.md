# Cosmic Stars Prototype FFI Data Contract

## 1. Scope

This document defines the **FFI boundary data contract** for Cosmic Stars Prototype.
The current Python prototype may output NDJSON for validation, but the schema below is the contract to preserve across the future Rust FFI boundary.

## 2. Units and Reference Frame

- Distance: `AU`
- Time: `yr`
- Velocity: `AU/yr`
- Mass (input model): `Msun`
- Output reference frame: `barycentric`

## 3. Canonical Frame Order

Canonical sequence:

1. First frame: `meta`
2. Following frames: `snapshot`
3. First snapshot: `tick = 0`
4. `tick` is monotonically increasing

When transport is not NDJSON (for example raw structs across FFI), preserve the same logical ordering.

## 4. Frame Definitions

## 4.1 `meta` frame

Required fields:

- `frame_type` (string): `"meta"`
- `units` (object): `{"distance":"AU","velocity":"AU/yr","time":"yr"}`
- `reference_frame` (string): `"barycentric"`
- `body_names` (string array): body names by index

Example:

```json
{
  "frame_type": "meta",
  "units": {"distance": "AU", "velocity": "AU/yr", "time": "yr"},
  "reference_frame": "barycentric",
  "body_names": ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
}
```

## 4.2 `snapshot` frame

Required fields:

- `frame_type` (string): `"snapshot"`
- `tick` (int): logical simulation step id
- `sim_time_yr` (float): authoritative simulation time (years)
- `time_scale_yr_per_real_sec` (float): current time-scale hint
- `reference_frame` (string): `"barycentric"`
- `bodies` (array): body state records

Optional field:

- `mercury_perihelion_longitude_deg` (float)

Example:

```json
{
  "frame_type": "snapshot",
  "tick": 120,
  "sim_time_yr": 1.2,
  "time_scale_yr_per_real_sec": 1.0,
  "reference_frame": "barycentric",
  "bodies": [
    {
      "id": 0,
      "name": "Sun",
      "position_au": [0.0, 0.0, 0.0],
      "velocity_au_per_yr": [0.0, 0.0, 0.0],
      "distance_from_barycenter_au": 0.0,
      "speed_au_per_yr": 0.0
    }
  ]
}
```

## 4.3 `body` object

Required fields:

- `id` (int): stable index in the frame
- `name` (string): body name
- `position_au` (float[3]): `[x, y, z]`
- `velocity_au_per_yr` (float[3]): `[vx, vy, vz]`
- `distance_from_barycenter_au` (float): `|r|`
- `speed_au_per_yr` (float): `|v|`

## 5. Time Semantics

- `sim_time_yr` is backend-authoritative simulation time.
- `tick` is the deterministic ordering key.
- `time_scale_yr_per_real_sec` is a frontend sync hint (pacing/interpolation).

Frontend/render side should use backend `tick + sim_time_yr` as authority.

## 6. Compatibility Rules (Rust FFI Migration)

- Do not rename/remove existing fields without explicit contract versioning.
- New fields must be additive and optional.
- Keep the same unit system unless introducing a major schema version.
- Preserve logical frame ordering (`meta` then `snapshot`).

## 7. Validation Checklist

Before accepting a new backend implementation:

1. `meta` frame exists and is first
2. `snapshot.tick` is monotonic
3. `body_names` aligns with body indexing
4. Numeric values are finite
5. Unit labels remain `AU`, `AU/yr`, `yr`

