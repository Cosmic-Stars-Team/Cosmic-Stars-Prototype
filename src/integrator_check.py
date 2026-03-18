import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List

from main_simulation import create_simulation


def angular_momentum_norm(sim) -> float:
    """计算系统总角动量模长（质心系下应近似守恒）。"""
    lx = 0.0
    ly = 0.0
    lz = 0.0
    for p in sim.particles:
        lx += p.m * (p.y * p.vz - p.z * p.vy)
        ly += p.m * (p.z * p.vx - p.x * p.vz)
        lz += p.m * (p.x * p.vy - p.y * p.vx)
    return math.sqrt(lx * lx + ly * ly + lz * lz)


def body_state(sim, name: str) -> Dict[str, float]:
    """
    读取单个天体状态:
    - 坐标/速度: 质心系
    - 轨道根数: 相对太阳（若该天体不是太阳）
    """
    p = sim.particles[name]
    r = math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z)
    v = math.sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz)
    row = {
        "name": name,
        "x": p.x,
        "y": p.y,
        "z": p.z,
        "vx": p.vx,
        "vy": p.vy,
        "vz": p.vz,
        "r": r,
        "v": v,
        "a": float("nan"),
        "e": float("nan"),
        "inc_deg": float("nan"),
    }

    if name != "Sun":
        sun = sim.particles["Sun"]
        orb = p.orbit(primary=sun)
        row["a"] = orb.a
        row["e"] = orb.e
        row["inc_deg"] = orb.inc * 180.0 / math.pi
    return row


def resolve_sample_bodies(names: List[str], raw: str) -> List[str]:
    """解析命令行传入的天体名称列表。"""
    if not raw.strip():
        return [n for n in names if n != "Sun"]

    selected = [x.strip() for x in raw.split(",") if x.strip()]
    valid = set(names)
    bad = [n for n in selected if n not in valid]
    if bad:
        raise ValueError(f"无效天体名称: {bad}，可选: {names}")
    return selected


def print_brief(title: str, rows: Iterable[Dict[str, float]]) -> None:
    """打印一组天体的关键参数，便于快速目视检查。"""
    print(f"\n{title}")
    print("name        a(AU)         e           x(AU)         y(AU)         z(AU)")
    for row in rows:
        a_str = "nan" if math.isnan(row["a"]) else f"{row['a']:.9f}"
        e_str = "nan" if math.isnan(row["e"]) else f"{row['e']:.9f}"
        print(
            f"{row['name']:<10}  {a_str:>12}  {e_str:>11}  "
            f"{row['x']:>12.6f}  {row['y']:>12.6f}  {row['z']:>12.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="快速检查积分器精度：输出半长轴/坐标/速度，并统计守恒量漂移。"
    )
    parser.add_argument("--years", type=float, default=10.0, help="模拟总时长（年）")
    parser.add_argument("--steps", type=int, default=20000, help="积分步数")
    parser.add_argument(
        "--integrator",
        type=str,
        default="mercurius",
        help="积分器名称，例如 mercurius / ias15 / whfast",
    )
    parser.add_argument(
        "--use-reboundx",
        action="store_true",
        help="启用REBOUNDx（例如水星GR近日点进动）",
    )
    parser.add_argument(
        "--sample-bodies",
        type=str,
        default="Mercury,Earth",
        help="要采样输出的天体名，逗号分隔；留空表示全部行星",
    )
    parser.add_argument(
        "--sample-stride",
        type=int,
        default=200,
        help="每多少积分步记录一次采样",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default="data/gen/integrator_samples.csv",
        help="采样输出CSV路径",
    )
    args = parser.parse_args()

    if args.steps <= 0:
        raise ValueError("steps 必须 > 0")
    if args.sample_stride <= 0:
        raise ValueError("sample_stride 必须 > 0")

    sim, names = create_simulation(use_reboundx=args.use_reboundx)
    sim.integrator = args.integrator

    # 固定步长积分器（如whfast）需要显式设置dt。
    dt = args.years / args.steps
    if args.integrator.lower() in {"mercurius", "whfast", "leapfrog", "sei"}:
        sim.dt = dt

    sample_bodies = resolve_sample_bodies(names, args.sample_bodies)
    print(
        f"Integrator={sim.integrator}, years={args.years}, steps={args.steps}, dt={dt:.6e} yr"
    )
    print(f"Sample bodies={sample_bodies}")
    print(f"Output CSV={args.output_csv}")

    init_rows = [body_state(sim, n) for n in sample_bodies]
    print_brief("Initial snapshot", init_rows)

    e0 = sim.energy()
    l0 = angular_momentum_norm(sim)

    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "tick",
                "sim_time_yr",
                "name",
                "x",
                "y",
                "z",
                "vx",
                "vy",
                "vz",
                "r",
                "v",
                "a",
                "e",
                "inc_deg",
            ],
        )
        writer.writeheader()

        # 写入初始帧
        for n in sample_bodies:
            row = body_state(sim, n)
            row["tick"] = 0
            row["sim_time_yr"] = sim.t
            writer.writerow(row)

        # 主积分循环
        for i in range(args.steps):
            target_t = sim.t + dt
            sim.integrate(target_t)
            tick = i + 1

            if tick % args.sample_stride == 0 or tick == args.steps:
                for n in sample_bodies:
                    row = body_state(sim, n)
                    row["tick"] = tick
                    row["sim_time_yr"] = sim.t
                    writer.writerow(row)

    e1 = sim.energy()
    l1 = angular_momentum_norm(sim)
    rel_e = abs((e1 - e0) / e0) if e0 != 0 else float("nan")
    rel_l = abs((l1 - l0) / l0) if l0 != 0 else float("nan")

    final_rows = [body_state(sim, n) for n in sample_bodies]
    print_brief("Final snapshot", final_rows)

    print("\nConservation drift (越小通常越好):")
    print(f"relative_energy_drift     = {rel_e:.3e}")
    print(f"relative_ang_mom_drift    = {rel_l:.3e}")


if __name__ == "__main__":
    main()
