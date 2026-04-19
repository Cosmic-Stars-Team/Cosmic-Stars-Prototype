import argparse
import math
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import rebound

try:
    import OpenEXR
except ImportError:  # pragma: no cover - exercised only when dependency is missing.
    OpenEXR = None


DEFAULT_WIDTH = 4096
DEFAULT_HEIGHT = 4096
DEFAULT_RS = 1.0
DEFAULT_BOUNDARY_RADIUS_RS = 15.0
DEFAULT_CLUSTER_STRENGTH = 8.1
DEFAULT_EPSILON_SCALE = 1.0e-6
DEFAULT_IAS15_EPSILON = 1.0e-10
DEFAULT_IAS15_INITIAL_DT_SCALE = 0.5


def critical_impact_parameter(rs: float) -> float:
    """Photon-sphere critical impact parameter for Schwarzschild in units of Rs."""
    return 1.5 * math.sqrt(3.0) * rs


def map_pixel_y_to_u(
    y: int,
    *,
    height: int,
    rs: float,
    boundary_radius_rs: float = DEFAULT_BOUNDARY_RADIUS_RS,
) -> float:
    """
    Maps texture Y to u = 1 / r.

    v = y / (height - 1)
    u(v) = v * (1 / (boundary_radius_rs * rs))
    """
    if height < 2:
        raise ValueError("height must be at least 2")
    if rs <= 0.0:
        raise ValueError("rs must be > 0")
    if not (0 <= y < height):
        raise ValueError(f"y={y} is out of range for height={height}")

    v = y / (height - 1)
    return v / (boundary_radius_rs * rs)


def map_pixel_x_to_b(
    x: int,
    *,
    width: int,
    b_crit: float,
    b_max: float,
    epsilon: float,
    cluster_strength: float = DEFAULT_CLUSTER_STRENGTH,
) -> float:
    """
    Exponential mapping that clusters samples near b_crit.

    t = x / (width - 1)
    s(t) = (exp(k * t) - 1) / (exp(k) - 1)
    b(t) = (b_crit + epsilon) + (b_max - (b_crit + epsilon)) * s(t)

    With k = 8.1, roughly the first 80% of texels cover only ~19.8% of the
    physical [b_min, b_max] interval.
    """
    if width < 2:
        raise ValueError("width must be at least 2")
    if cluster_strength <= 0.0:
        raise ValueError("cluster_strength must be > 0")
    if epsilon <= 0.0:
        raise ValueError("epsilon must be > 0")
    if b_max <= b_crit + epsilon:
        raise ValueError("b_max must be larger than b_crit + epsilon")
    if not (0 <= x < width):
        raise ValueError(f"x={x} is out of range for width={width}")

    t = x / (width - 1)
    b_min = b_crit + epsilon
    span = b_max - b_min
    scaled = math.expm1(cluster_strength * t) / math.expm1(cluster_strength)
    return b_min + span * scaled


def map_b_to_uv_x(
    b: float,
    *,
    b_crit: float,
    b_max: float,
    epsilon: float,
    cluster_strength: float = DEFAULT_CLUSTER_STRENGTH,
) -> float:
    """
    Inverse of map_pixel_x_to_b in normalized UV space.

    n = (b - (b_crit + epsilon)) / (b_max - (b_crit + epsilon))
    uv.x = ln(1 + n * (exp(k) - 1)) / k
    """
    if cluster_strength <= 0.0:
        raise ValueError("cluster_strength must be > 0")
    if epsilon <= 0.0:
        raise ValueError("epsilon must be > 0")
    b_min = b_crit + epsilon
    if b_max <= b_min:
        raise ValueError("b_max must be larger than b_crit + epsilon")
    if not (b_min <= b <= b_max):
        raise ValueError(f"b={b} is outside [{b_min}, {b_max}]")

    normalized = (b - b_min) / (b_max - b_min)
    return math.log1p(normalized * math.expm1(cluster_strength)) / cluster_strength


def _unwrap_angle(previous_angle: float, current_angle: float) -> float:
    delta = current_angle - previous_angle
    while delta <= -math.pi:
        current_angle += 2.0 * math.pi
        delta = current_angle - previous_angle
    while delta > math.pi:
        current_angle -= 2.0 * math.pi
        delta = current_angle - previous_angle
    return current_angle


def _ray_state_from_b_u(b: float, u: float, rs: float) -> tuple[tuple[float, float], tuple[float, float]]:
    if b <= 0.0:
        raise ValueError("b must be > 0")
    if u <= 0.0:
        raise ValueError("u must be > 0")
    if rs < 0.0:
        raise ValueError("rs must be >= 0")

    radius = 1.0 / u
    v_tangent = b * u
    v_radial_sq = 1.0 - (b * b * u * u) + (rs * b * b * u * u * u)
    if v_radial_sq < 0.0:
        raise ValueError("non-escaping configuration for the supplied (b, u, rs)")

    v_radial = -math.sqrt(v_radial_sq)
    return (radius, 0.0), (v_radial, v_tangent)


def _direction_angle(vx: float, vy: float) -> float:
    return math.atan2(vy, vx)


def _make_ray_simulation(
    *,
    b: float,
    start_u: float,
    rs: float,
    boundary_radius: float,
    ias15_epsilon: float,
    ias15_initial_dt_scale: float,
) -> rebound.Simulation:
    position, velocity = _ray_state_from_b_u(b, start_u, rs)
    l2 = b * b

    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.ri_ias15.epsilon = ias15_epsilon
    sim.ri_ias15.min_dt = 0.0
    sim.dt = max(
        1.0e-3 * max(rs, 1.0),
        ias15_initial_dt_scale * (position[0] - boundary_radius),
    )
    sim.add(m=0.0, x=position[0], y=position[1], vx=velocity[0], vy=velocity[1])

    def schwarzschild_additional_force(sim_pointer):
        local_sim = sim_pointer.contents
        particle = local_sim.particles[0]
        radius_sq = particle.x * particle.x + particle.y * particle.y + particle.z * particle.z
        radius = math.sqrt(radius_sq)
        scale = -1.5 * rs * l2 / (radius_sq * radius_sq * radius)
        particle.ax += scale * particle.x
        particle.ay += scale * particle.y
        particle.az += scale * particle.z

    sim.additional_forces = schwarzschild_additional_force
    return sim


def _trace_column_deflections(
    u_targets: np.ndarray,
    *,
    b: float,
    rs: float,
    boundary_radius_rs: float,
    ias15_epsilon: float = DEFAULT_IAS15_EPSILON,
    ias15_initial_dt_scale: float = DEFAULT_IAS15_INITIAL_DT_SCALE,
) -> np.ndarray:
    if u_targets.ndim != 1:
        raise ValueError("u_targets must be a 1D array")
    if len(u_targets) == 0:
        return np.zeros(0, dtype=np.float64)
    if np.any(u_targets < 0.0):
        raise ValueError("u_targets must be >= 0")
    if np.any(np.diff(u_targets) < 0.0):
        raise ValueError("u_targets must be sorted in ascending order")
    if rs == 0.0:
        return np.zeros_like(u_targets, dtype=np.float64)
    if ias15_epsilon <= 0.0:
        raise ValueError("ias15_epsilon must be > 0")
    if ias15_initial_dt_scale <= 0.0:
        raise ValueError("ias15_initial_dt_scale must be > 0")

    boundary_radius = boundary_radius_rs * rs
    boundary_u = 1.0 / boundary_radius
    if np.any(u_targets > boundary_u):
        raise ValueError("u_targets must stay in the far-field domain outside the boundary sphere")

    out = np.zeros_like(u_targets, dtype=np.float64)
    positive_indices = np.flatnonzero(u_targets > 0.0)
    if len(positive_indices) == 0:
        return out

    finite_u_targets = u_targets[positive_indices]
    start_u = float(finite_u_targets[0])
    start_radius = 1.0 / start_u
    if start_radius <= boundary_radius:
        return out

    sim = _make_ray_simulation(
        b=b,
        start_u=start_u,
        rs=rs,
        boundary_radius=boundary_radius,
        ias15_epsilon=ias15_epsilon,
        ias15_initial_dt_scale=ias15_initial_dt_scale,
    )

    particle = sim.particles[0]
    previous_x = particle.x
    previous_y = particle.y
    previous_vx = particle.vx
    previous_vy = particle.vy
    previous_radius = math.hypot(previous_x, previous_y)
    previous_angle = _direction_angle(previous_vx, previous_vy)

    finite_angles = np.zeros(len(finite_u_targets), dtype=np.float64)
    finite_angles[0] = previous_angle
    next_target = 1
    boundary_angle = previous_angle

    while True:
        sim.step()
        particle = sim.particles[0]
        current_x = particle.x
        current_y = particle.y
        current_vx = particle.vx
        current_vy = particle.vy
        current_radius = math.hypot(current_x, current_y)
        current_angle = _unwrap_angle(previous_angle, _direction_angle(current_vx, current_vy))

        while next_target < len(finite_u_targets):
            target_radius = 1.0 / float(finite_u_targets[next_target])
            if not (current_radius <= target_radius <= previous_radius):
                break
            blend = (previous_radius - target_radius) / (previous_radius - current_radius)
            finite_angles[next_target] = previous_angle + blend * (current_angle - previous_angle)
            next_target += 1

        if current_radius <= boundary_radius:
            if previous_radius == current_radius:
                boundary_angle = current_angle
            else:
                blend = (previous_radius - boundary_radius) / (previous_radius - current_radius)
                boundary_angle = previous_angle + blend * (current_angle - previous_angle)
            break

        previous_x = current_x
        previous_y = current_y
        previous_vx = current_vx
        previous_vy = current_vy
        previous_radius = current_radius
        previous_angle = current_angle

    for local_idx, global_idx in enumerate(positive_indices):
        out[global_idx] = boundary_angle - finite_angles[local_idx]
    if u_targets[0] == 0.0:
        out[0] = boundary_angle - math.pi
    return out


def _generate_column_range(
    *,
    x_start: int,
    x_stop: int,
    width: int,
    u_values: np.ndarray,
    rs: float,
    boundary_radius_rs: float,
    b_crit: float,
    b_max: float,
    epsilon: float,
    cluster_strength: float,
    ias15_epsilon: float,
    ias15_initial_dt_scale: float,
) -> tuple[int, np.ndarray]:
    chunk = np.zeros((len(u_values), x_stop - x_start), dtype=np.float32)
    for local_x, x in enumerate(range(x_start, x_stop)):
        b = map_pixel_x_to_b(
            x,
            width=width,
            b_crit=b_crit,
            b_max=b_max,
            epsilon=epsilon,
            cluster_strength=cluster_strength,
        )
        column = _trace_column_deflections(
            u_values,
            b=b,
            rs=rs,
            boundary_radius_rs=boundary_radius_rs,
            ias15_epsilon=ias15_epsilon,
            ias15_initial_dt_scale=ias15_initial_dt_scale,
        )
        chunk[:, local_x] = column.astype(np.float32, copy=False)
    return x_start, chunk


def _resolve_worker_count(workers: int | None) -> int:
    if workers is None:
        return 1
    if workers <= 0:
        return max(1, os.cpu_count() or 1)
    return workers


def _build_column_chunks(width: int, workers: int) -> list[tuple[int, int]]:
    chunk_count = min(width, workers)
    base = width // chunk_count
    remainder = width % chunk_count
    chunks: list[tuple[int, int]] = []
    start = 0
    for i in range(chunk_count):
        stop = start + base + (1 if i < remainder else 0)
        chunks.append((start, stop))
        start = stop
    return chunks


def calculate_deflection(
    b: float,
    u: float,
    *,
    rs: float,
    ias15_epsilon: float = DEFAULT_IAS15_EPSILON,
    boundary_radius_rs: float = DEFAULT_BOUNDARY_RADIUS_RS,
    ias15_initial_dt_scale: float = DEFAULT_IAS15_INITIAL_DT_SCALE,
) -> float:
    if u < 0.0:
        raise ValueError("u must be >= 0")
    if rs < 0.0:
        raise ValueError("rs must be >= 0")
    if rs == 0.0:
        return 0.0
    boundary_u = 1.0 / (boundary_radius_rs * rs)
    if u >= boundary_u:
        return 0.0

    values = _trace_column_deflections(
        np.array([u, boundary_u], dtype=np.float64),
        b=b,
        rs=rs,
        boundary_radius_rs=boundary_radius_rs,
        ias15_epsilon=ias15_epsilon,
        ias15_initial_dt_scale=ias15_initial_dt_scale,
    )
    return float(values[0])


def generate_lut(
    *,
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
    rs: float = DEFAULT_RS,
    boundary_radius_rs: float = DEFAULT_BOUNDARY_RADIUS_RS,
    b_max: float | None = None,
    cluster_strength: float = DEFAULT_CLUSTER_STRENGTH,
    epsilon_scale: float = DEFAULT_EPSILON_SCALE,
    ias15_epsilon: float = DEFAULT_IAS15_EPSILON,
    ias15_initial_dt_scale: float = DEFAULT_IAS15_INITIAL_DT_SCALE,
    workers: int | None = 1,
) -> np.ndarray:
    if width < 2 or height < 2:
        raise ValueError("width and height must be at least 2")
    if rs <= 0.0:
        raise ValueError("rs must be > 0")
    if epsilon_scale <= 0.0:
        raise ValueError("epsilon_scale must be > 0")
    if ias15_epsilon <= 0.0:
        raise ValueError("ias15_epsilon must be > 0")
    if ias15_initial_dt_scale <= 0.0:
        raise ValueError("ias15_initial_dt_scale must be > 0")
    resolved_workers = _resolve_worker_count(workers)
    if resolved_workers <= 0:
        raise ValueError("workers must resolve to a positive integer")

    b_crit = critical_impact_parameter(rs)
    epsilon = epsilon_scale * rs
    resolved_b_max = b_max if b_max is not None else boundary_radius_rs * rs
    if resolved_b_max <= b_crit + epsilon:
        raise ValueError("resolved b_max must exceed b_crit + epsilon")

    lut = np.zeros((height, width), dtype=np.float32)
    u_values = np.array(
        [
            map_pixel_y_to_u(
                y,
                height=height,
                rs=rs,
                boundary_radius_rs=boundary_radius_rs,
            )
            for y in range(height)
        ],
        dtype=np.float64,
    )

    if resolved_workers == 1:
        _, chunk = _generate_column_range(
            x_start=0,
            x_stop=width,
            width=width,
            u_values=u_values,
            rs=rs,
            boundary_radius_rs=boundary_radius_rs,
            b_crit=b_crit,
            b_max=resolved_b_max,
            epsilon=epsilon,
            cluster_strength=cluster_strength,
            ias15_epsilon=ias15_epsilon,
            ias15_initial_dt_scale=ias15_initial_dt_scale,
        )
        lut[:, :] = chunk
        return lut

    chunk_ranges = _build_column_chunks(width, resolved_workers)
    with ProcessPoolExecutor(
        max_workers=resolved_workers,
        mp_context=get_context("spawn"),
    ) as executor:
        futures = [
            executor.submit(
                _generate_column_range,
                x_start=x_start,
                x_stop=x_stop,
                width=width,
                u_values=u_values,
                rs=rs,
                boundary_radius_rs=boundary_radius_rs,
                b_crit=b_crit,
                b_max=resolved_b_max,
                epsilon=epsilon,
                cluster_strength=cluster_strength,
                ias15_epsilon=ias15_epsilon,
                ias15_initial_dt_scale=ias15_initial_dt_scale,
            )
            for x_start, x_stop in chunk_ranges
        ]
        for future in futures:
            x_start, chunk = future.result()
            x_stop = x_start + chunk.shape[1]
            lut[:, x_start:x_stop] = chunk

    return lut


def write_exr(output_path: str | Path, lut: np.ndarray, *, rgb: bool = False) -> Path:
    if OpenEXR is None:
        raise RuntimeError(
            "OpenEXR is required for EXR export. Install it with `python -m pip install OpenEXR`."
        )

    if lut.ndim != 2:
        raise ValueError("lut must be a 2D array")

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    plane = np.ascontiguousarray(lut.astype(np.float32, copy=False))

    header = {
        "compression": OpenEXR.ZIP_COMPRESSION,
        "type": OpenEXR.scanlineimage,
    }

    if rgb:
        zero = np.zeros_like(plane)
        channels = {"R": plane, "G": zero, "B": zero}
    else:
        channels = {"Y": plane}

    exr_file = OpenEXR.File(header, channels)
    exr_file.write(str(output))
    return output


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Bake a 2D EXR LUT for far-field relativistic ray deflection."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/gen/far_field_ray_deflection_lut_4096.exr",
        help="Output EXR path.",
    )
    parser.add_argument("--width", type=int, default=DEFAULT_WIDTH)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--rs", type=float, default=DEFAULT_RS)
    parser.add_argument(
        "--boundary-radius-rs",
        type=float,
        default=DEFAULT_BOUNDARY_RADIUS_RS,
        help="Bounding sphere radius expressed in Rs.",
    )
    parser.add_argument(
        "--b-max-rs",
        type=float,
        default=DEFAULT_BOUNDARY_RADIUS_RS,
        help="Maximum impact parameter expressed in Rs.",
    )
    parser.add_argument(
        "--cluster-strength",
        type=float,
        default=DEFAULT_CLUSTER_STRENGTH,
        help="Exponential clustering strength for the X axis mapping.",
    )
    parser.add_argument(
        "--epsilon-scale",
        type=float,
        default=DEFAULT_EPSILON_SCALE,
        help="Epsilon multiplier applied to Rs so b starts at b_crit + epsilon.",
    )
    parser.add_argument(
        "--rgb",
        action="store_true",
        help="Write an RGB EXR and place the LUT only in the R channel.",
    )
    parser.add_argument(
        "--ias15-epsilon",
        type=float,
        default=DEFAULT_IAS15_EPSILON,
        help="IAS15 precision control parameter for the geodesic quadrature.",
    )
    parser.add_argument(
        "--ias15-initial-dt-scale",
        type=float,
        default=DEFAULT_IAS15_INITIAL_DT_SCALE,
        help="Scale factor used to seed IAS15's initial timestep from the shell spacing.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Number of worker processes. Use 0 to auto-detect from CPU count.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    lut = generate_lut(
        width=args.width,
        height=args.height,
        rs=args.rs,
        boundary_radius_rs=args.boundary_radius_rs,
        b_max=args.b_max_rs * args.rs,
        cluster_strength=args.cluster_strength,
        epsilon_scale=args.epsilon_scale,
        ias15_epsilon=args.ias15_epsilon,
        ias15_initial_dt_scale=args.ias15_initial_dt_scale,
        workers=args.workers,
    )
    output = write_exr(args.output, lut, rgb=args.rgb)
    print(f"Wrote LUT EXR to {output}")


if __name__ == "__main__":
    main()
