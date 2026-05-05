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

from blackbody_1d_lut_baker import (
    DEFAULT_MAX_TEMPERATURE_K,
    DEFAULT_PRIMARY_SEGMENT_END_K,
    DEFAULT_PRIMARY_SEGMENT_FRACTION,
    DEFAULT_SECONDARY_SEGMENT_END_K,
    DEFAULT_SECONDARY_SEGMENT_FRACTION,
    DEFAULT_WIDTH,
    generate_lut,
    map_pixel_x_to_temperature,
    map_temperature_to_uv_x,
    temperature_to_linear_srgb,
    write_exr,
)


class Blackbody1DLutBakerTests(unittest.TestCase):
    def test_x_axis_maps_piecewise_to_temperature_range(self) -> None:
        width = 41
        primary_segment_x = int((width - 1) * DEFAULT_PRIMARY_SEGMENT_FRACTION)
        secondary_segment_x = int((width - 1) * DEFAULT_SECONDARY_SEGMENT_FRACTION)
        primary_segment_mid_x = primary_segment_x // 2
        secondary_segment_mid_x = (primary_segment_x + secondary_segment_x) // 2
        tertiary_segment_mid_x = (secondary_segment_x + (width - 1)) // 2

        self.assertEqual(
            map_pixel_x_to_temperature(0, width=width),
            0.0,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(width - 1, width=width),
            DEFAULT_MAX_TEMPERATURE_K,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(primary_segment_x, width=width),
            DEFAULT_PRIMARY_SEGMENT_END_K,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(secondary_segment_x, width=width),
            DEFAULT_SECONDARY_SEGMENT_END_K,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(primary_segment_mid_x, width=width),
            DEFAULT_PRIMARY_SEGMENT_END_K * 0.5,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(secondary_segment_mid_x, width=width),
            DEFAULT_PRIMARY_SEGMENT_END_K + (
                DEFAULT_SECONDARY_SEGMENT_END_K - DEFAULT_PRIMARY_SEGMENT_END_K
            ) * 0.5,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(tertiary_segment_mid_x, width=width),
            DEFAULT_SECONDARY_SEGMENT_END_K + (
                DEFAULT_MAX_TEMPERATURE_K - DEFAULT_SECONDARY_SEGMENT_END_K
            ) * 0.5,
        )

    def test_temperature_to_uv_uses_matching_piecewise_segments(self) -> None:
        self.assertAlmostEqual(map_temperature_to_uv_x(0.0), 0.0)
        self.assertAlmostEqual(
            map_temperature_to_uv_x(DEFAULT_PRIMARY_SEGMENT_END_K),
            DEFAULT_PRIMARY_SEGMENT_FRACTION,
        )
        self.assertAlmostEqual(
            map_temperature_to_uv_x(DEFAULT_SECONDARY_SEGMENT_END_K),
            DEFAULT_SECONDARY_SEGMENT_FRACTION,
        )
        self.assertAlmostEqual(map_temperature_to_uv_x(DEFAULT_MAX_TEMPERATURE_K), 1.0)
        self.assertAlmostEqual(
            map_temperature_to_uv_x(DEFAULT_PRIMARY_SEGMENT_END_K * 0.5),
            DEFAULT_PRIMARY_SEGMENT_FRACTION * 0.5,
        )
        self.assertAlmostEqual(
            map_temperature_to_uv_x(16000.0),
            DEFAULT_PRIMARY_SEGMENT_FRACTION
            + (DEFAULT_SECONDARY_SEGMENT_FRACTION - DEFAULT_PRIMARY_SEGMENT_FRACTION)
            * 0.5,
        )

    def test_zero_kelvin_is_black(self) -> None:
        rgb = temperature_to_linear_srgb(0.0)
        self.assertTrue(np.allclose(rgb, np.zeros(3), atol=1.0e-8))

    def test_sub_visible_temperatures_fade_toward_black(self) -> None:
        faded = temperature_to_linear_srgb(500.0)
        visible = temperature_to_linear_srgb(1000.0)

        self.assertGreater(float(np.max(faded)), 0.0)
        self.assertLess(float(np.max(faded)), 1.0)
        self.assertAlmostEqual(float(np.max(visible)), 1.0, places=6)

    def test_visible_temperatures_are_normalized_to_unit_peak(self) -> None:
        for temperature in (1000.0, 3000.0, 6500.0, 10000.0, 40000.0):
            rgb = temperature_to_linear_srgb(temperature)
            self.assertTrue(np.all(rgb >= 0.0))
            self.assertAlmostEqual(float(np.max(rgb)), 1.0, places=6)

    def test_generate_lut_returns_single_row_rgb_image(self) -> None:
        lut = generate_lut(width=16)

        self.assertEqual(lut.shape, (1, 16, 3))
        self.assertEqual(lut.dtype, np.float32)
        self.assertTrue(np.isfinite(lut).all())

    def test_write_exr_outputs_rgb_channels(self) -> None:
        lut = generate_lut(width=16)

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "blackbody.exr"
            write_exr(output_path, lut)

            self.assertTrue(output_path.exists())
            exr_file = OpenEXR.File(str(output_path))
            channels = exr_file.channels()
            self.assertTrue(
                "RGB" in channels or {"R", "G", "B"}.issubset(channels.keys())
            )


if __name__ == "__main__":
    unittest.main()
