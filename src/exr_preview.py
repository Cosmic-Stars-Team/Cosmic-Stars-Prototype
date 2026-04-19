import argparse
from pathlib import Path

import numpy as np
import OpenEXR
from PIL import Image


def load_exr_image(path: str | Path) -> np.ndarray:
    exr = OpenEXR.File(str(path))
    channels = exr.channels()

    if "Y" in channels:
        data = channels["Y"].pixels
        return np.asarray(data, dtype=np.float32)
    if {"R", "G", "B"}.issubset(channels.keys()):
        return np.stack(
            [
                np.asarray(channels["R"].pixels, dtype=np.float32),
                np.asarray(channels["G"].pixels, dtype=np.float32),
                np.asarray(channels["B"].pixels, dtype=np.float32),
            ],
            axis=-1,
        )
    if "RGB" in channels:
        return np.asarray(channels["RGB"].pixels, dtype=np.float32)
    if "R" in channels:
        data = channels["R"].pixels
        return np.asarray(data, dtype=np.float32)

    raise ValueError(f"Unsupported channel layout: {list(channels.keys())}")


def load_single_channel_exr(path: str | Path) -> np.ndarray:
    image = load_exr_image(path)
    if image.ndim == 2:
        return image
    if image.ndim == 3 and image.shape[2] == 3:
        return np.asarray(image[:, :, 0], dtype=np.float32)
    raise ValueError("single-channel preview expects a 2D image or RGB image")


def normalize_for_preview(image: np.ndarray) -> np.ndarray:
    if image.ndim != 2:
        raise ValueError("preview normalization expects a 2D single-channel image")

    min_value = float(np.min(image))
    max_value = float(np.max(image))
    if max_value <= min_value:
        return np.zeros_like(image, dtype=np.uint8)

    normalized = (image - min_value) / (max_value - min_value)
    return np.clip(normalized * 255.0, 0.0, 255.0).astype(np.uint8)


def normalize_rgb_for_preview(image: np.ndarray) -> np.ndarray:
    if image.ndim != 3 or image.shape[2] != 3:
        raise ValueError("RGB preview normalization expects an image with shape (height, width, 3)")

    clipped = np.clip(image, 0.0, None)
    max_value = float(np.max(clipped))
    if max_value <= 0.0:
        return np.zeros_like(clipped, dtype=np.uint8)

    normalized = clipped / max_value
    return np.clip(normalized * 255.0, 0.0, 255.0).astype(np.uint8)


def write_preview_png(
    image: np.ndarray,
    output_path: str | Path,
    *,
    repeat_height: int = 128,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if repeat_height < 1:
        raise ValueError("repeat_height must be >= 1")

    if image.ndim == 2:
        preview = normalize_for_preview(image)
        Image.fromarray(preview, mode="L").save(output)
        return output

    if image.ndim == 3 and image.shape[2] == 3:
        preview = normalize_rgb_for_preview(image)
        if preview.shape[0] == 1 and repeat_height > 1:
            preview = np.repeat(preview, repeat_height, axis=0)
        Image.fromarray(preview, mode="RGB").save(output)
        return output

    raise ValueError("preview export expects either a 2D single-channel image or a 3D RGB image")

    return output


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert an EXR LUT into a PNG preview. 1D RGB LUTs are expanded into a visible color strip."
    )
    parser.add_argument("input", type=str, help="Input EXR path.")
    parser.add_argument(
        "--output",
        type=str,
        default="",
        help="Optional PNG output path. Defaults to <input>.preview.png",
    )
    parser.add_argument(
        "--repeat-height",
        type=int,
        default=128,
        help="Vertical expansion used when previewing a 1-pixel-tall RGB LUT.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path.with_suffix(".preview.png")

    image = load_exr_image(input_path)
    output = write_preview_png(image, output_path, repeat_height=args.repeat_height)
    print(f"Preview written to {output}")
    print(f"min={float(image.min())} max={float(image.max())}")


if __name__ == "__main__":
    main()
