import json
import math
from typing import Dict, Iterator, List, Optional

import rebound
import reboundx

from database import load_solar_system

# 光速（单位：AU/yr），用于REBOUNDx中的GR项
C_AU_PER_YR = 63239.7263


def enable_mercury_perihelion_precession(sim: rebound.Simulation) -> reboundx.Extras:
    """通过REBOUNDx的广义相对论修正启用水星近日点进动。"""
    rebx = reboundx.Extras(sim)
    gr = rebx.load_force("gr")
    rebx.add_force(gr)
    gr.params["c"] = C_AU_PER_YR
    # 防止rebx对象被回收
    sim.rebx = rebx
    return rebx


def create_simulation(data=None, use_reboundx: bool = True):
    """
    从JSON数据创建rebound模拟

    参数:
        data: 包含太阳系数据的字典，如果为None则从文件加载

    返回:
        (rebound.Simulation对象, 粒子名称列表)
    """
    # 如果没有提供数据，从JSON文件加载
    if data is None:
        data = load_solar_system()

    # 创建模拟对象
    sim = rebound.Simulation()

    # 设置单位系统：年、天文单位、太阳质量
    sim.units = ("yr", "AU", "Msun")
    sim.integrator = "ias15"

    # 添加太阳
    star = data["star"]
    sim.add(
        m=star["mass"],
        hash=star["name"],
        x=star["x"],
        y=star["y"],
        z=star["z"],
        vx=star["vx"],
        vy=star["vy"],
        vz=star["vz"],
    )

    # 添加所有行星
    for planet in data["planets"]:
        sim.add(
            m=planet["mass"],
            a=planet["a"],  # 半长轴
            e=planet["e"],  # 离心率
            inc=planet["inc"],  # 轨道倾角
            Omega=planet["Omega"],  # 升交点经度
            omega=planet["omega"],  # 近心点幅角
            M=planet["M"],  # 平近点角
            hash=planet["name"],  # 行星名称（作为hash值）
        )

    # 移动到质心系（确保系统动量为零）
    sim.move_to_com()

    if use_reboundx:
        enable_mercury_perihelion_precession(sim)

    # 创建名称列表
    names = [star["name"]] + [p["name"] for p in data["planets"]]

    return sim, names


def _mercury_perihelion_longitude_deg(sim: rebound.Simulation) -> Optional[float]:
    try:
        sun = sim.particles["Sun"]
        mercury = sim.particles["Mercury"]
    except Exception:
        return None

    orbit = mercury.orbit(primary=sun)
    return (orbit.Omega + orbit.omega) * 180.0 / math.pi


def _build_body_state(particle, body_id: int, name: str) -> Dict:
    """构建单个天体的状态字典。"""
    r = math.sqrt(particle.x**2 + particle.y**2 + particle.z**2)
    v = math.sqrt(particle.vx**2 + particle.vy**2 + particle.vz**2)
    return {
        "id": body_id,
        "name": name,
        "position_au": [particle.x, particle.y, particle.z],
        "velocity_au_per_yr": [particle.vx, particle.vy, particle.vz],
        "distance_from_barycenter_au": r,
        "speed_au_per_yr": v,
    }


def _build_meta_frame(names: List[str]) -> Dict:
    """构建输出流首帧元信息。"""
    return {
        "frame_type": "meta",
        "units": {"distance": "AU", "velocity": "AU/yr", "time": "yr"},
        "reference_frame": "barycentric",
        "body_names": names,
    }


def _build_snapshot(
    sim: rebound.Simulation,
    names: List[str],
    tick: int,
    time_scale_yr_per_real_sec: float,
) -> Dict:
    """构建前端可消费的一帧快照。"""
    bodies = [
        _build_body_state(
            particle=particle,
            body_id=i,
            name=names[i] if i < len(names) else f"body_{i}",
        )
        for i, particle in enumerate(sim.particles)
    ]

    snapshot = {
        "frame_type": "snapshot",
        "tick": tick,
        "sim_time_yr": sim.t,
        "time_scale_yr_per_real_sec": time_scale_yr_per_real_sec,
        "reference_frame": "barycentric",
        "bodies": bodies,
    }

    varpi = _mercury_perihelion_longitude_deg(sim)
    if varpi is not None:
        snapshot["mercury_perihelion_longitude_deg"] = varpi

    return snapshot


def _iter_snapshot_frames(
    sim: rebound.Simulation,
    names: List[str],
    years: float,
    steps: int,
    snapshot_stride: int,
    time_scale_yr_per_real_sec: float,
) -> Iterator[Dict]:
    """
    只负责积分和采样，不负责IO。
    产出顺序：初始帧 -> 若干采样帧。
    """
    dt = years / steps
    yield _build_snapshot(
        sim,
        names,
        tick=0,
        time_scale_yr_per_real_sec=time_scale_yr_per_real_sec,
    )

    for i in range(steps):
        sim.integrate(sim.t + dt)
        tick = i + 1
        if tick % snapshot_stride == 0 or tick == steps:
            yield _build_snapshot(
                sim,
                names,
                tick=tick,
                time_scale_yr_per_real_sec=time_scale_yr_per_real_sec,
            )


def run_simulation(
    sim,
    names,
    years=100,
    steps=1000,
    snapshot_stride=1,
    output_path="simulation_stream.jsonl",
    time_scale_yr_per_real_sec=1.0,
):
    """
    运行模拟并输出前端可消费的NDJSON快照流

    参数:
        sim: rebound.Simulation对象
        names: 粒子名称列表
        years: 模拟年数
        steps: 积分步数
        snapshot_stride: 每隔多少积分步输出一帧快照
        output_path: 输出jsonl文件路径
        time_scale_yr_per_real_sec: 时间倍率（每1秒真实时间推进多少模拟年）
    """
    if steps <= 0:
        raise ValueError("steps 必须大于 0")
    if snapshot_stride <= 0:
        raise ValueError("snapshot_stride 必须大于 0")

    print(f"模拟开始：{sim.N} 个天体，模拟时长 {years} 年")
    print(f"输出格式：NDJSON -> {output_path}")

    meta = _build_meta_frame(names)
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(json.dumps(meta, ensure_ascii=False) + "\n")
        for frame in _iter_snapshot_frames(
            sim=sim,
            names=names,
            years=years,
            steps=steps,
            snapshot_stride=snapshot_stride,
            time_scale_yr_per_real_sec=time_scale_yr_per_real_sec,
        ):
            f.write(json.dumps(frame, ensure_ascii=False) + "\n")
            written += 1

    print(f"模拟完成：共写入 {written} 帧快照")
    return output_path


if __name__ == "__main__":
    # 从JSON文件加载数据并创建模拟
    data = load_solar_system()
    sim, names = create_simulation(data)

    # 运行10年的模拟，输出供前端消费的数据流
    run_simulation(
        sim,
        names,
        years=10,
        steps=1000,
        snapshot_stride=10,
        output_path="simulation_stream.jsonl",
        time_scale_yr_per_real_sec=1.0,
    )
