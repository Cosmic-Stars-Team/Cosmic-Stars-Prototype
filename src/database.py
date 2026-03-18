import json
from pathlib import Path
from typing import Dict, Optional


def _default_solar_system_path() -> Path:
    """默认使用仓库根目录 data/solar_system.json。"""
    return Path(__file__).resolve().parents[1] / "data" / "solar_system.json"


def load_solar_system(filepath: Optional[str] = None) -> Dict:
    path = Path(filepath) if filepath else _default_solar_system_path()
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)
