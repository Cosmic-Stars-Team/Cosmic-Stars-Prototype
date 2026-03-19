# Cosmic Stars Prototype Protocol

## 1. Scope

This document defines the backend output contract for **Cosmic Stars Prototype**.
Current transport is file-based NDJSON (`data/gen/simulation_stream.jsonl`), but the same frame schema is intended for future IPC/pipe/shared-memory transport.

## 2. Units and Frame of Reference

- Distance: `AU`
- Time: `yr`
- Velocity: `AU/yr`
- Mass (input model): `Msun`
- Reference frame for state output: `barycentric`

## 3. Stream Format

The stream is NDJSON: **one JSON object per line**.

Required order:

1. First line: `meta` frame
2. Following lines: `snapshot` frames
3. First snapshot must be `tick = 0`
4. `tick` must be monotonically increasing

## 4. Frame Types

## 4.1 `meta` frame

Required fields:

- `frame_type` (string): must be `"meta"`
- `units` (object): `{"distance":"AU","velocity":"AU/yr","time":"yr"}`
- `reference_frame` (string): currently `"barycentric"`
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

- `frame_type` (string): must be `"snapshot"`
- `tick` (int): simulation logical step id
- `sim_time_yr` (float): authoritative simulation time in years
- `time_scale_yr_per_real_sec` (float): current time-scale hint
- `reference_frame` (string): currently `"barycentric"`
- `bodies` (array): body state records

Optional fields:

- `mercury_perihelion_longitude_deg` (float): present when Mercury/Sun orbit can be resolved

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

## 4.3 `body` object in `bodies`

Required fields:

- `id` (int): stable index in this frame
- `name` (string): body name
- `position_au` (float[3]): `[x, y, z]`
- `velocity_au_per_yr` (float[3]): `[vx, vy, vz]`
- `distance_from_barycenter_au` (float): `|r|`
- `speed_au_per_yr` (float): `|v|`

## 5. Time Semantics

- `sim_time_yr` is backend-authoritative simulation time.
- `tick` is an integer logical step counter for deterministic ordering.
- `time_scale_yr_per_real_sec` is a synchronization hint for frontend pacing/interpolation.

Frontend should not infer simulation authority from wall-clock; use `tick + sim_time_yr` from backend.

## 6. Compatibility Rules (for Rust Migration)

- Existing fields must not be renamed or removed without protocol versioning.
- New fields must be additive and optional.
- Unit system must remain unchanged unless a major protocol version is introduced.
- Frame ordering (`meta` first, then snapshots) must be preserved.

## 7. Validation Checklist

Before releasing a new backend implementation:

1. `meta` frame present and first
2. `snapshot.tick` monotonic
3. `body_names` length consistent with body indexing
4. Numeric fields parse to finite values
5. Unit labels unchanged (`AU`, `AU/yr`, `yr`)

