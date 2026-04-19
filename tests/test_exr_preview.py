import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image

import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from exr_preview import load_exr_image, load_single_channel_exr, normalize_for_preview, write_preview_png
from blackbody_1d_lut_baker import write_exr as write_rgb_exr
from far_field_lut_baker import write_exr


class ExrPreviewTests(unittest.TestCase):
    def test_load_exr_image_reads_rgb_lut(self) -> None:
        source = np.array(
            [[[1.0, 0.0, 0.0], [0.0, 0.5, 1.0]]],
            dtype=np.float32,
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            exr_path = Path(tmp_dir) / "source_rgb.exr"
            write_rgb_exr(exr_path, source)

            loaded = load_exr_image(exr_path)

            self.assertEqual(loaded.shape, source.shape)
            self.assertTrue(np.allclose(loaded, source))

    def test_load_single_channel_exr_reads_written_lut(self) -> None:
        source = np.array([[0.0, 1.0], [2.0, 3.0]], dtype=np.float32)

        with tempfile.TemporaryDirectory() as tmp_dir:
            exr_path = Path(tmp_dir) / "source.exr"
            write_exr(exr_path, source, rgb=False)

            loaded = load_single_channel_exr(exr_path)

            self.assertEqual(loaded.shape, source.shape)
            self.assertTrue(np.allclose(loaded, source))

    def test_normalize_for_preview_spreads_values_to_full_range(self) -> None:
        image = np.array([[2.0, 4.0], [6.0, 10.0]], dtype=np.float32)

        normalized = normalize_for_preview(image)

        self.assertEqual(normalized.dtype, np.uint8)
        self.assertEqual(int(normalized.min()), 0)
        self.assertEqual(int(normalized.max()), 255)

    def test_write_preview_png_outputs_png_file(self) -> None:
        image = np.array([[0.0, 1.0], [2.0, 3.0]], dtype=np.float32)

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "preview.png"
            write_preview_png(image, output_path)

            self.assertTrue(output_path.exists())
            with Image.open(output_path) as loaded:
                self.assertEqual(loaded.size, (2, 2))

    def test_write_preview_png_expands_single_row_rgb_lut(self) -> None:
        image = np.array(
            [[[1.0, 0.0, 0.0], [0.0, 0.5, 1.0], [1.0, 1.0, 1.0]]],
            dtype=np.float32,
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "preview_rgb.png"
            write_preview_png(image, output_path, repeat_height=4)

            self.assertTrue(output_path.exists())
            with Image.open(output_path) as loaded:
                self.assertEqual(loaded.size, (3, 4))
                pixels = np.asarray(loaded)
                self.assertTrue(np.all(pixels[0] == pixels[-1]))
                self.assertTrue(np.array_equal(pixels[0, 0], np.array([255, 0, 0], dtype=np.uint8)))


if __name__ == "__main__":
    unittest.main()
