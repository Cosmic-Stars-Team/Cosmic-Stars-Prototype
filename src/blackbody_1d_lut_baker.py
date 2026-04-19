import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np

try:
    import OpenEXR
except ImportError:  # pragma: no cover - exercised only when dependency is missing.
    OpenEXR = None

try:
    import colour
except ImportError:  # pragma: no cover - exercised only when dependency is missing.
    colour = None


DEFAULT_WIDTH = 4096
DEFAULT_HEIGHT = 1
DEFAULT_MIN_TEMPERATURE_K = 0.0
DEFAULT_MAX_TEMPERATURE_K = 40000.0
DEFAULT_FADE_TO_BLACK_END_K = 1000.0
DEFAULT_OUTPUT = "data/gen/blackbody_1d_lut_4k.exr"


def _require_colour() -> None:
    if colour is None:
        raise RuntimeError(
            "colour-science is required for blackbody baking. "
            "Install it with `python -m pip install colour-science`."
        )


@lru_cache(maxsize=1)
def _spectral_shape():
    _require_colour()
    return colour.SpectralShape(360, 780, 1)


@lru_cache(maxsize=1)
def _cie_1931_cmfs():
    _require_colour()
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy()
    return cmfs.align(_spectral_shape())


def smoothstep(edge0: float, edge1: float, x: float) -> float:
    if edge1 <= edge0:
        raise ValueError("edge1 must be greater than edge0")

    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return float(t * t * (3.0 - 2.0 * t))


def map_pixel_x_to_temperature(
    x: int,
    *,
    width: int,
    min_temperature_k: float = DEFAULT_MIN_TEMPERATURE_K,
    max_temperature_k: float = DEFAULT_MAX_TEMPERATURE_K,
) -> float:
    """
    Map texture X to temperature using a quadratic curve so lower temperatures
    receive more texel density than the visually flatter high-temperature tail.
    """
    if width < 2:
        raise ValueError("width must be at least 2")
    if max_temperature_k < min_temperature_k:
        raise ValueError("max_temperature_k must be >= min_temperature_k")
    if not (0 <= x < width):
        raise ValueError(f"x={x} is out of range for width={width}")

    u = x / (width - 1)
    return min_temperature_k + (max_temperature_k - min_temperature_k) * (u * u)


def temperature_to_cie_xy(temperature_k: float) -> np.ndarray:
    """
    Convert a blackbody temperature to exact CIE 1931 xy chromaticity coordinates
    by integrating the Planck spectral distribution against the CIE 1931 2-degree
    standard observer.
    """
    _require_colour()

    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be > 0")

    spectral_distribution = colour.sd_blackbody(temperature_k, _spectral_shape())
    xyz = colour.sd_to_XYZ(
        spectral_distribution,
        cmfs=_cie_1931_cmfs(),
        method="Integration",
    )
    xyz = np.asarray(xyz, dtype=np.float64) / 100.0
    return np.asarray(colour.XYZ_to_xy(xyz), dtype=np.float64)


def temperature_to_linear_srgb(
    temperature_k: float,
    *,
    fade_to_black_end_k: float = DEFAULT_FADE_TO_BLACK_END_K,
) -> np.ndarray:
    """
    Convert temperature to linear sRGB chromaticity, normalized so the peak
    channel equals 1.0. Temperatures below fade_to_black_end_k are smoothly
    attenuated toward black to suppress negligible visible emission.
    """
    if fade_to_black_end_k <= 0.0:
        raise ValueError("fade_to_black_end_k must be > 0")
    if temperature_k <= 0.0:
        return np.zeros(3, dtype=np.float32)

    # Below the visible threshold we only need a smooth fade toward black, not a
    # numerically fragile evaluation of an almost-dark spectrum.
    effective_temperature_k = max(temperature_k, fade_to_black_end_k)
    xy = temperature_to_cie_xy(effective_temperature_k)
    xyz_chromaticity = np.asarray(colour.xy_to_XYZ(xy), dtype=np.float64)
    linear_rgb = np.asarray(
        colour.XYZ_to_RGB(
            xyz_chromaticity,
            "sRGB",
            chromatic_adaptation_transform=None,
            apply_cctf_encoding=False,
        ),
        dtype=np.float64,
    )

    # Negative channels are outside the display gamut; clamp before normalizing
    # so the LUT remains texture-friendly for the Godot shader.
    linear_rgb = np.clip(linear_rgb, 0.0, None)
    peak = float(np.max(linear_rgb))
    if peak > 0.0:
        linear_rgb /= peak
    else:
        linear_rgb[:] = 0.0

    fade = smoothstep(0.0, fade_to_black_end_k, temperature_k)
    return (linear_rgb * fade).astype(np.float32)


def generate_lut(
    *,
    width: int = DEFAULT_WIDTH,
    min_temperature_k: float = DEFAULT_MIN_TEMPERATURE_K,
    max_temperature_k: float = DEFAULT_MAX_TEMPERATURE_K,
    fade_to_black_end_k: float = DEFAULT_FADE_TO_BLACK_END_K,
) -> np.ndarray:
    if width < 2:
        raise ValueError("width must be at least 2")
    if max_temperature_k < min_temperature_k:
        raise ValueError("max_temperature_k must be >= min_temperature_k")

    lut = np.zeros((DEFAULT_HEIGHT, width, 3), dtype=np.float32)
    for x in range(width):
        temperature_k = map_pixel_x_to_temperature(
            x,
            width=width,
            min_temperature_k=min_temperature_k,
            max_temperature_k=max_temperature_k,
        )
        lut[0, x, :] = temperature_to_linear_srgb(
            temperature_k,
            fade_to_black_end_k=fade_to_black_end_k,
        )

    return lut


def write_exr(output_path: str | Path, lut: np.ndarray) -> Path:
    if OpenEXR is None:
        raise RuntimeError(
            "OpenEXR is required for EXR export. Install it with `python -m pip install OpenEXR`."
        )
    if lut.ndim != 3 or lut.shape[2] != 3:
        raise ValueError("lut must be an RGB image with shape (height, width, 3)")

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    image = np.ascontiguousarray(lut.astype(np.float32, copy=False))
    channels = {
        "R": np.ascontiguousarray(image[:, :, 0]),
        "G": np.ascontiguousarray(image[:, :, 1]),
        "B": np.ascontiguousarray(image[:, :, 2]),
    }
    header = {
        "compression": OpenEXR.ZIP_COMPRESSION,
        "type": OpenEXR.scanlineimage,
    }

    exr_file = OpenEXR.File(header, channels)
    exr_file.write(str(output))
    return output


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Bake a 1D EXR LUT that maps blackbody temperature to normalized linear sRGB."
    )
    parser.add_argument("--output", type=str, default=DEFAULT_OUTPUT, help="Output EXR path.")
    parser.add_argument("--width", type=int, default=DEFAULT_WIDTH, help="LUT width in pixels.")
    parser.add_argument(
        "--min-temperature-k",
        type=float,
        default=DEFAULT_MIN_TEMPERATURE_K,
        help="Minimum temperature mapped to the left edge of the LUT.",
    )
    parser.add_argument(
        "--max-temperature-k",
        type=float,
        default=DEFAULT_MAX_TEMPERATURE_K,
        help="Maximum temperature mapped to the right edge of the LUT.",
    )
    parser.add_argument(
        "--fade-to-black-end-k",
        type=float,
        default=DEFAULT_FADE_TO_BLACK_END_K,
        help="Temperature where the black fade reaches full chromaticity.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    lut = generate_lut(
        width=args.width,
        min_temperature_k=args.min_temperature_k,
        max_temperature_k=args.max_temperature_k,
        fade_to_black_end_k=args.fade_to_black_end_k,
    )
    output = write_exr(args.output, lut)
    print(f"Wrote blackbody LUT EXR to {output}")


if __name__ == "__main__":
    main()
