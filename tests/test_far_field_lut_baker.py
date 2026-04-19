import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import OpenEXR

import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from far_field_lut_baker import (
    DEFAULT_BOUNDARY_RADIUS_RS,
    DEFAULT_CLUSTER_STRENGTH,
    calculate_deflection,
    critical_impact_parameter,
    generate_lut,
    map_b_to_uv_x,
    map_pixel_x_to_b,
    map_pixel_y_to_u,
    write_exr,
)


def _initial_ray_state(b: float, u: float, rs: float) -> tuple[np.ndarray, np.ndarray]:
    r = 1.0 / u
    vt = b * u
    vr_sq = 1.0 - (b * b * u * u) + (rs * b * b * u * u * u)
    if vr_sq < 0.0:
        raise ValueError("non-escaping configuration in reference integrator")
    vr = -math.sqrt(vr_sq)
    position = np.array([r, 0.0], dtype=np.float64)
    velocity = np.array([vr, vt], dtype=np.float64)
    return position, velocity


def _schwarzschild_accel(position: np.ndarray, rs: float, l2: float) -> np.ndarray:
    radius_sq = float(np.dot(position, position))
    radius = math.sqrt(radius_sq)
    scale = -1.5 * rs * l2 / (radius_sq * radius_sq * radius)
    return position * scale


def _rk4_step(position: np.ndarray, velocity: np.ndarray, dt: float, rs: float, l2: float) -> tuple[np.ndarray, np.ndarray]:
    def deriv(pos: np.ndarray, vel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return vel, _schwarzschild_accel(pos, rs, l2)

    k1x, k1v = deriv(position, velocity)
    k2x, k2v = deriv(position + 0.5 * dt * k1x, velocity + 0.5 * dt * k1v)
    k3x, k3v = deriv(position + 0.5 * dt * k2x, velocity + 0.5 * dt * k2v)
    k4x, k4v = deriv(position + dt * k3x, velocity + dt * k3v)

    new_position = position + (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
    new_velocity = velocity + (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
    return new_position, new_velocity


def reference_far_field_deflection(
    b: float,
    u: float,
    rs: float,
    *,
    boundary_radius_rs: float = DEFAULT_BOUNDARY_RADIUS_RS,
    dt: float = 1.0e-3,
) -> float:
    if rs == 0.0:
        return 0.0

    boundary_radius = boundary_radius_rs * rs
    radius = 1.0 / u
    if radius <= boundary_radius:
        return 0.0

    position, velocity = _initial_ray_state(b, u, rs)
    initial_theta = math.atan2(velocity[1], velocity[0])
    previous_position = position.copy()
    previous_velocity = velocity.copy()

    while math.hypot(position[0], position[1]) > boundary_radius:
        previous_position = position.copy()
        previous_velocity = velocity.copy()
        position, velocity = _rk4_step(position, velocity, dt, rs, b * b)

    prev_radius = math.hypot(previous_position[0], previous_position[1])
    curr_radius = math.hypot(position[0], position[1])
    blend = (prev_radius - boundary_radius) / (prev_radius - curr_radius)
    boundary_velocity = previous_velocity + blend * (velocity - previous_velocity)
    boundary_theta = math.atan2(boundary_velocity[1], boundary_velocity[0])
    return boundary_theta - initial_theta


class FarFieldLutBakerTests(unittest.TestCase):
    def test_y_axis_maps_linearly_to_inverse_distance(self) -> None:
        rs = 2.0
        height = 4096
        max_u = 1.0 / (DEFAULT_BOUNDARY_RADIUS_RS * rs)

        self.assertEqual(map_pixel_y_to_u(0, height=height, rs=rs), 0.0)
        self.assertAlmostEqual(
            map_pixel_y_to_u(height - 1, height=height, rs=rs),
            max_u,
        )

    def test_x_axis_maps_endpoints_correctly(self) -> None:
        rs = 1.0
        width = 4096
        epsilon = 1.0e-6
        b_crit = critical_impact_parameter(rs)
        b_max = DEFAULT_BOUNDARY_RADIUS_RS * rs

        self.assertAlmostEqual(
            map_pixel_x_to_b(
                0,
                width=width,
                b_crit=b_crit,
                b_max=b_max,
                epsilon=epsilon,
            ),
            b_crit + epsilon,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_b(
                width - 1,
                width=width,
                b_crit=b_crit,
                b_max=b_max,
                epsilon=epsilon,
            ),
            b_max,
        )

    def test_x_axis_maps_linearly_across_the_physical_interval(self) -> None:
        rs = 1.0
        width = 4096
        epsilon = 1.0e-6
        b_crit = critical_impact_parameter(rs)
        b_max = DEFAULT_BOUNDARY_RADIUS_RS * rs
        b_min = b_crit + epsilon
        span = b_max - b_min
        x_80 = int(round(0.8 * (width - 1)))

        mapped_b = map_pixel_x_to_b(
            x_80,
            width=width,
            b_crit=b_crit,
            b_max=b_max,
            epsilon=epsilon,
            cluster_strength=DEFAULT_CLUSTER_STRENGTH,
        )
        normalized_span = (mapped_b - b_min) / span

        self.assertAlmostEqual(normalized_span, 0.8, places=9)

    def test_inverse_x_mapping_round_trips(self) -> None:
        rs = 1.0
        width = 4096
        epsilon = 1.0e-6
        b_crit = critical_impact_parameter(rs)
        b_max = DEFAULT_BOUNDARY_RADIUS_RS * rs

        for x in (0, 17, width // 4, int(0.8 * (width - 1)), width - 1):
            b = map_pixel_x_to_b(
                x,
                width=width,
                b_crit=b_crit,
                b_max=b_max,
                epsilon=epsilon,
            )
            uv_x = map_b_to_uv_x(
                b,
                b_crit=b_crit,
                b_max=b_max,
                epsilon=epsilon,
            )
            expected_uv_x = x / (width - 1)
            self.assertAlmostEqual(uv_x, expected_uv_x, places=9)

    def test_calculate_deflection_is_zero_in_flat_space(self) -> None:
        self.assertAlmostEqual(calculate_deflection(b=6.0, u=0.04, rs=0.0), 0.0, places=8)

    def test_calculate_deflection_matches_reference_integral(self) -> None:
        b = 6.5
        u = 1.0 / 20.0
        rs = 1.0

        expected = reference_far_field_deflection(b=b, u=u, rs=rs)
        actual = calculate_deflection(b=b, u=u, rs=rs)

        self.assertAlmostEqual(actual, expected, delta=2.0e-4)

    def test_calculate_deflection_is_zero_on_boundary_shell(self) -> None:
        rs = 1.0
        u_boundary = 1.0 / (DEFAULT_BOUNDARY_RADIUS_RS * rs)

        self.assertAlmostEqual(
            calculate_deflection(b=10.0, u=u_boundary, rs=rs),
            0.0,
            places=8,
        )

    def test_generate_lut_allocates_expected_float_buffer(self) -> None:
        lut = generate_lut(width=8, height=4, rs=1.0)

        self.assertEqual(lut.shape, (4, 8))
        self.assertEqual(lut.dtype, np.float32)
        self.assertTrue(np.isfinite(lut).all())

    def test_generate_lut_parallel_matches_sequential(self) -> None:
        sequential = generate_lut(width=8, height=4, rs=1.0, workers=1)
        parallel = generate_lut(width=8, height=4, rs=1.0, workers=2)

        self.assertTrue(np.allclose(sequential, parallel, atol=1.0e-6))

    def test_write_exr_outputs_single_channel_float_image(self) -> None:
        lut = generate_lut(width=8, height=4, rs=1.0)

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "lut.exr"
            write_exr(output_path, lut, rgb=False)

            self.assertTrue(output_path.exists())
            exr_file = OpenEXR.File(str(output_path))
            self.assertIn("Y", exr_file.channels())


if __name__ == "__main__":
    unittest.main()
