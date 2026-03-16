import json
import math

import rebound
from database import load_solar_system

import reboundx

# 光速（单位：AU/yr），用于REBOUNDx中的GR项
C_AU_PER_YR = 63239.7263


def enable_mercury_perihelion_precession(sim):
    """通过REBOUNDx的广义相对论修正启用水星近日点进动。"""
    rebx = reboundx.Extras(sim)
    gr = rebx.load_force("gr")
    rebx.add_force(gr)
    gr.params["c"] = C_AU_PER_YR
    # 防止rebx对象被回收
    sim.rebx = rebx
    return rebx


def create_simulation(data=None, use_reboundx=True):
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


def _print_particle_state(name, particle):
    r = math.sqrt(particle.x**2 + particle.y**2 + particle.z**2)
    v = math.sqrt(particle.vx**2 + particle.vy**2 + particle.vz**2)
    print(
        f"  {name}: "
        f"坐标(x,y,z)=({particle.x:.6f}, {particle.y:.6f}, {particle.z:.6f}) AU, "
        f"速度(vx,vy,vz)=({particle.vx:.6f}, {particle.vy:.6f}, {particle.vz:.6f}) AU/yr, "
        f"|r|={r:.6f} AU, |v|={v:.6f} AU/yr"
    )


def _print_mercury_perihelion(sim):
    sun = sim.particles["Sun"]
    mercury = sim.particles["Mercury"]
    orbit = mercury.orbit(primary=sun)
    varpi_deg = (orbit.Omega + orbit.omega) * 180.0 / math.pi
    print(f"  Mercury近日点经度(varpi): {varpi_deg:.6f} deg")


def _mercury_perihelion_longitude_deg(sim):
    try:
        sun = sim.particles["Sun"]
        mercury = sim.particles["Mercury"]
    except Exception:
        return None

    orbit = mercury.orbit(primary=sun)
    return (orbit.Omega + orbit.omega) * 180.0 / math.pi


def _build_snapshot(sim, names, tick, time_scale_yr_per_real_sec):
    bodies = []
    for i, particle in enumerate(sim.particles):
        name = names[i] if i < len(names) else f"body_{i}"
        r = math.sqrt(particle.x**2 + particle.y**2 + particle.z**2)
        v = math.sqrt(particle.vx**2 + particle.vy**2 + particle.vz**2)
        bodies.append(
            {
                "id": i,
                "name": name,
                "position_au": [particle.x, particle.y, particle.z],
                "velocity_au_per_yr": [particle.vx, particle.vy, particle.vz],
                "distance_from_barycenter_au": r,
                "speed_au_per_yr": v,
            }
        )

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

    # 进行积分
    dt = years / steps
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        meta = {
            "frame_type": "meta",
            "units": {"distance": "AU", "velocity": "AU/yr", "time": "yr"},
            "reference_frame": "barycentric",
            "body_names": names,
        }
        f.write(json.dumps(meta, ensure_ascii=False) + "\n")

        # 初始帧
        f.write(
            json.dumps(
                _build_snapshot(
                    sim,
                    names,
                    tick=0,
                    time_scale_yr_per_real_sec=time_scale_yr_per_real_sec,
                ),
                ensure_ascii=False,
            )
            + "\n"
        )
        written += 1

        for i in range(steps):
            sim.integrate(sim.t + dt)
            tick = i + 1
            if tick % snapshot_stride == 0 or tick == steps:
                snapshot = _build_snapshot(
                    sim,
                    names,
                    tick=tick,
                    time_scale_yr_per_real_sec=time_scale_yr_per_real_sec,
                )
                f.write(json.dumps(snapshot, ensure_ascii=False) + "\n")
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
