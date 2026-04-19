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
    DEFAULT_WIDTH,
    generate_lut,
    map_pixel_x_to_temperature,
    temperature_to_linear_srgb,
    write_exr,
)


class Blackbody1DLutBakerTests(unittest.TestCase):
    def test_x_axis_maps_quadratically_to_temperature_range(self) -> None:
        width = 5

        self.assertEqual(
            map_pixel_x_to_temperature(0, width=width),
            0.0,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(width - 1, width=width),
            DEFAULT_MAX_TEMPERATURE_K,
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature(1, width=width),
            DEFAULT_MAX_TEMPERATURE_K * (0.25 ** 2),
        )
        self.assertAlmostEqual(
            map_pixel_x_to_temperature((width - 1) // 2, width=width),
            10000.0,
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
