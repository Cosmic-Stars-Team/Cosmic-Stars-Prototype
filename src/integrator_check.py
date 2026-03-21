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
        "varpi_deg": float("nan"),
    }

    if name != "Sun":
        sun = sim.particles["Sun"]
        orb = p.orbit(primary=sun)
        row["a"] = orb.a
        row["e"] = orb.e
        row["inc_deg"] = orb.inc * 180.0 / math.pi
        row["varpi_deg"] = (orb.Omega + orb.omega) * 180.0 / math.pi
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


def perihelion_rate_arcsec_per_century(
    times_yr: List[float], angles_rad: List[float]
) -> float:
    """
    根据时间序列估计近日点进动速率（arcsec/century）。
    angles_rad 允许存在 2π 跳变，本函数会先做展开再拟合斜率。
    """
    n = len(times_yr)
    if n < 2 or n != len(angles_rad):
        return float("nan")

    # unwrap angles in-place into a new list
    unwrapped = [angles_rad[0]]
    for a in angles_rad[1:]:
        current = a
        prev = unwrapped[-1]
        delta = current - prev
        while delta > math.pi:
            current -= 2.0 * math.pi
            delta = current - prev
        while delta < -math.pi:
            current += 2.0 * math.pi
            delta = current - prev
        unwrapped.append(current)

    # least-squares slope
    x_mean = sum(times_yr) / n
    y_mean = sum(unwrapped) / n
    num = 0.0
    den = 0.0
    for x, y in zip(times_yr, unwrapped):
        dx = x - x_mean
        num += dx * (y - y_mean)
        den += dx * dx
    if den == 0.0:
        return float("nan")

    slope_rad_per_yr = num / den
    return slope_rad_per_yr * 180.0 / math.pi * 3600.0 * 100.0


def build_compare_paths(output_csv: Path) -> tuple[Path, Path]:
    """根据输出路径生成开/关 REBOUNDx 的对照文件名。"""
    parent = output_csv.parent
    stem = output_csv.stem
    suffix = output_csv.suffix or ".csv"
    no_rebx = parent / f"{stem}_no_reboundx{suffix}"
    with_rebx = parent / f"{stem}_with_reboundx{suffix}"
    return no_rebx, with_rebx


def run_case(
    args: argparse.Namespace, use_reboundx: bool, output_csv: Path, case_name: str
) -> float:
    sim, names = create_simulation(use_reboundx=use_reboundx)
    sim.integrator = args.integrator

    dt = args.years / args.steps
    if args.integrator.lower() in {"mercurius", "whfast", "leapfrog", "sei"}:
        sim.dt = dt

    sample_bodies = resolve_sample_bodies(names, args.sample_bodies)
    print(
        f"\n=== {case_name} ===\n"
        f"Integrator={sim.integrator}, years={args.years}, steps={args.steps}, dt={dt:.6e} yr"
    )
    print(f"Sample bodies={sample_bodies}")
    print(f"Output CSV={output_csv}")

    init_rows = [body_state(sim, n) for n in sample_bodies]
    print_brief("Initial snapshot", init_rows)

    e0 = sim.energy()
    l0 = angular_momentum_norm(sim)

    times_yr: List[float] = []
    varpi_rad: List[float] = []
    get_varpi = None
    if args.report_mercury_perihelion:
        sun = sim.particles["Sun"]
        mercury = sim.particles["Mercury"]

        def _get_varpi() -> float:
            orb = mercury.orbit(primary=sun)
            return orb.Omega + orb.omega

        get_varpi = _get_varpi
        times_yr.append(sim.t)
        varpi_rad.append(get_varpi())

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
                "varpi_deg",
            ],
        )
        writer.writeheader()

        for n in sample_bodies:
            row = body_state(sim, n)
            row["tick"] = 0
            row["sim_time_yr"] = sim.t
            writer.writerow(row)

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

                if get_varpi is not None:
                    times_yr.append(sim.t)
                    varpi_rad.append(get_varpi())

    e1 = sim.energy()
    l1 = angular_momentum_norm(sim)
    rel_e = abs((e1 - e0) / e0) if e0 != 0 else float("nan")
    rel_l = abs((l1 - l0) / l0) if l0 != 0 else float("nan")

    final_rows = [body_state(sim, n) for n in sample_bodies]
    print_brief("Final snapshot", final_rows)

    print("\nConservation drift (越小通常越好):")
    print(f"relative_energy_drift     = {rel_e:.3e}")
    print(f"relative_ang_mom_drift    = {rel_l:.3e}")

    rate = float("nan")
    if get_varpi is not None:
        rate = perihelion_rate_arcsec_per_century(times_yr, varpi_rad)
        print(f"mercury_perihelion_rate_arcsec_per_century = {rate:.6f}")
    return rate


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
    parser.add_argument(
        "--report-mercury-perihelion",
        action="store_true",
        help="输出水星近日点长期漂移率（arcsec/century）",
    )
    parser.add_argument(
        "--compare-reboundx",
        action="store_true",
        help="在相同参数下分别运行开/关 REBOUNDx，并打印差值",
    )
    args = parser.parse_args()

    if args.steps <= 0:
        raise ValueError("steps 必须 > 0")
    if args.sample_stride <= 0:
        raise ValueError("sample_stride 必须 > 0")

    output_csv = Path(args.output_csv)
    if args.compare_reboundx:
        no_path, with_path = build_compare_paths(output_csv)
        no_rate = run_case(
            args=args,
            use_reboundx=False,
            output_csv=no_path,
            case_name="WITHOUT REBOUNDx",
        )
        with_rate = run_case(
            args=args,
            use_reboundx=True,
            output_csv=with_path,
            case_name="WITH REBOUNDx",
        )
        if args.report_mercury_perihelion:
            delta = with_rate - no_rate
            print("\nPerihelion compare:")
            print(f"delta_mercury_perihelion_arcsec_per_century = {delta:.6f}")
    else:
        run_case(
            args=args,
            use_reboundx=args.use_reboundx,
            output_csv=output_csv,
            case_name="SINGLE RUN",
        )


if __name__ == "__main__":
    main()
