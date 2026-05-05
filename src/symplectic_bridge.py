"""
Symplectic Bridge Integrator — 路线 3

H = H_main(Q,P) + H_sub(q,p) + H_cross(q, Q)

二阶 Strang splitting:
  1. exp(dt/2 · H_cross)  — kick: modify p and P from full state
  2. exp(dt · (H_main + H_sub))  — drift: integrate both REBOUND sims
  3. exp(dt/2 · H_cross)  — kick again at new state

关键性质:
  - H_cross 取决于完整相空间 (q, Q)，系统自治
  - kick 写回 main_sim.EMB 动量，闭环
  - 输出采样对齐到全局参考时间，保证收敛测试不污染
"""
import numpy as np
import rebound

G = 4.0 * np.pi ** 2

# ──────────────────────────────────────────────
# 工具
# ──────────────────────────────────────────────


def _pos(p):
    return np.array([p.x, p.y, p.z])


def _vel(p):
    return np.array([p.vx, p.vy, p.vz])


def _barycenter(sim):
    mt = 0.0
    r = np.zeros(3)
    for p in sim.particles:
        mt += p.m
        r += p.m * _pos(p)
    return r / mt


# ──────────────────────────────────────────────
# 潮汐力（H_cross 的梯度）
# ──────────────────────────────────────────────


def _sun_gravity(pos, sun_pos, M_sun=1.0):
    rv = pos - sun_pos
    return -G * M_sun * rv / np.linalg.norm(rv) ** 3


def tidal_force(sub_sim, main_sim):
    """
    当前状态下的潮汐力。

    Returns
    -------
    a_earth, a_moon : (3,) ndarray
        地球 / 月球受的潮汐加速度 (AU/yr²)
    a_emb : (3,) ndarray
        EMB 的反作用加速度 (动量守恒)
    """
    sp = _pos(main_sim.particles[0])
    ep = _pos(main_sim.particles[1])
    e = _pos(sub_sim.particles[0])
    m = _pos(sub_sim.particles[1])
    bc = _barycenter(sub_sim)

    a_e = _sun_gravity(ep + e, sp) - _sun_gravity(ep + bc, sp)
    a_m = _sun_gravity(ep + m, sp) - _sun_gravity(ep + bc, sp)

    mt = sub_sim.particles[0].m + sub_sim.particles[1].m
    a_emb = -(sub_sim.particles[0].m * a_e + sub_sim.particles[1].m * a_m) / mt

    return a_e, a_m, a_emb


# ──────────────────────────────────────────────
# Kick 算子: exp(Δt · H_cross)
# ──────────────────────────────────────────────


def _apply_cross_kick(sub_sim, main_sim, dt_half):
    a_e, a_m, a_emb = tidal_force(sub_sim, main_sim)

    sub_sim.particles[0].vx += a_e[0] * dt_half
    sub_sim.particles[0].vy += a_e[1] * dt_half
    sub_sim.particles[0].vz += a_e[2] * dt_half

    sub_sim.particles[1].vx += a_m[0] * dt_half
    sub_sim.particles[1].vy += a_m[1] * dt_half
    sub_sim.particles[1].vz += a_m[2] * dt_half

    main_sim.particles[1].vx += a_emb[0] * dt_half
    main_sim.particles[1].vy += a_emb[1] * dt_half
    main_sim.particles[1].vz += a_emb[2] * dt_half


# ──────────────────────────────────────────────
# 初始化
# ──────────────────────────────────────────────


def make_sims(sub_ratio=50):
    """
    创建 main_sim (Sun + EMB) 和 sub_sim (Earth + Moon).

    Parameters
    ----------
    sub_ratio : int
        子系步长 = 外步长 / sub_ratio
    """
    main_sim = rebound.Simulation()
    main_sim.units = ("AU", "yr", "Msun")
    main_sim.add(m=1.0)
    main_sim.add(m=3e-6, a=1.0)
    main_sim.move_to_com()
    main_sim.integrator = "wh"

    sub_sim = rebound.Simulation()
    sub_sim.units = ("AU", "yr", "Msun")
    sub_sim.add(m=3.0e-6 * 0.987)
    sub_sim.add(m=3.0e-6 * 0.013, a=0.00257)
    sub_sim.move_to_com()
    sub_sim.integrator = "wh"

    return main_sim, sub_sim


# ──────────────────────────────────────────────
# 单步推进
# ──────────────────────────────────────────────


def bridge_step(main_sim, sub_sim, dt_outer):
    """一次完整的二阶 KDK 步进 (dt_outer)."""
    _apply_cross_kick(sub_sim, main_sim, dt_outer * 0.5)

    target = main_sim.t + dt_outer
    sub_sim.integrate(target)
    main_sim.integrate(target)

    _apply_cross_kick(sub_sim, main_sim, dt_outer * 0.5)


# ──────────────────────────────────────────────
# 全文积分（对齐采样）
# ──────────────────────────────────────────────


def integrate(
    dt_outer,
    t_end=1.0,
    n_samples=50,
    sub_ratio=50,
    return_diagnostics=False,
):
    """
    Symplectic bridge 全文积分。

    Parameters
    ----------
    dt_outer : float
        外步长 (yr)
    t_end : float
        积分总时长 (yr)
    n_samples : int
        输出采样点数
    sub_ratio : int
        子系内步数 = dt_outer / sub_ratio
    return_diagnostics : bool
        返回时间 / 距离 / 能量等诊断数据

    Returns
    -------
    distances : ndarray (n_samples,)
        地月距离 (AU)
    若 return_diagnostics=True:
        (times, distances, emb_distances, energies)
    """
    main_sim, sub_sim = make_sims(sub_ratio)
    main_sim.dt = dt_outer
    sub_sim.dt = dt_outer / sub_ratio

    # 精确对齐的采样时间
    sample_times = np.linspace(0, t_end, n_samples + 1)[1:]

    distances = np.empty(n_samples)
    times = sample_times.copy()
    emb_dists = np.empty(n_samples) if return_diagnostics else None
    energies = np.empty(n_samples) if return_diagnostics else None

    for i, target_t in enumerate(sample_times):
        # 推进到精确目标时间
        while main_sim.t < target_t - 1e-14:
            step_dt = min(dt_outer, target_t - main_sim.t)
            bridge_step(main_sim, sub_sim, step_dt)

        e = _pos(sub_sim.particles[0])
        m = _pos(sub_sim.particles[1])
        distances[i] = np.linalg.norm(e - m)

        if return_diagnostics:
            emb = _pos(main_sim.particles[1])
            emb_dists[i] = np.linalg.norm(emb)
            # 简易能量诊断 (仅地月部分)
            v_e = _vel(sub_sim.particles[0])
            v_m = _vel(sub_sim.particles[1])
            me, mm = sub_sim.particles[0].m, sub_sim.particles[1].m
            K = 0.5 * me * np.dot(v_e, v_e) + 0.5 * mm * np.dot(v_m, v_m)
            U = -G * me * mm / distances[i]
            energies[i] = K + U

    if return_diagnostics:
        return times, distances, emb_dists, energies
    return distances


# ──────────────────────────────────────────────
# 原始 KDK（对照组，对齐采样）
# ──────────────────────────────────────────────


def integrate_original(dt_outer, t_end=1.0, n_samples=50):
    """
    原始单向 KDK 方案：潮汐力只写 sub_sim，读 main_sim 为外部驱动。
    使用对齐采样以保证比较公平。
    """
    main_sim, sub_sim = make_sims()
    main_sim.dt = dt_outer
    sub_sim.dt = dt_outer / 50

    sample_times = np.linspace(0, t_end, n_samples + 1)[1:]
    distances = np.empty(n_samples)

    for i, target_t in enumerate(sample_times):
        while main_sim.t < target_t - 1e-14:
            step_dt = min(dt_outer, target_t - main_sim.t)

            # --- 原始 KDK (单向写) ---
            sp = _pos(main_sim.particles[0])
            ep = _pos(main_sim.particles[1])
            bc = _barycenter(sub_sim)
            a_e = _sun_gravity(ep + _pos(sub_sim.particles[0]), sp) - \
                  _sun_gravity(ep + bc, sp)
            a_m = _sun_gravity(ep + _pos(sub_sim.particles[1]), sp) - \
                  _sun_gravity(ep + bc, sp)
            sub_sim.particles[0].vx += a_e[0] * step_dt * 0.5
            sub_sim.particles[0].vy += a_e[1] * step_dt * 0.5
            sub_sim.particles[0].vz += a_e[2] * step_dt * 0.5
            sub_sim.particles[1].vx += a_m[0] * step_dt * 0.5
            sub_sim.particles[1].vy += a_m[1] * step_dt * 0.5
            sub_sim.particles[1].vz += a_m[2] * step_dt * 0.5

            tgt = main_sim.t + step_dt
            sub_sim.integrate(tgt)
            main_sim.integrate(tgt)

            sp = _pos(main_sim.particles[0])
            ep = _pos(main_sim.particles[1])
            bc = _barycenter(sub_sim)
            a_e = _sun_gravity(ep + _pos(sub_sim.particles[0]), sp) - \
                  _sun_gravity(ep + bc, sp)
            a_m = _sun_gravity(ep + _pos(sub_sim.particles[1]), sp) - \
                  _sun_gravity(ep + bc, sp)
            sub_sim.particles[0].vx += a_e[0] * step_dt * 0.5
            sub_sim.particles[0].vy += a_e[1] * step_dt * 0.5
            sub_sim.particles[0].vz += a_e[2] * step_dt * 0.5
            sub_sim.particles[1].vx += a_m[0] * step_dt * 0.5
            sub_sim.particles[1].vy += a_m[1] * step_dt * 0.5
            sub_sim.particles[1].vz += a_m[2] * step_dt * 0.5

        e = _pos(sub_sim.particles[0])
        m = _pos(sub_sim.particles[1])
        distances[i] = np.linalg.norm(e - m)

    return distances


# ──────────────────────────────────────────────
# 收敛测试 + 对比
# ──────────────────────────────────────────────

if __name__ == "__main__":
    dts = [1 / 365, 1 / 730, 1 / 1460, 1 / 2920]

    # ── Bridge ──
    print("=" * 60)
    print("Symplectic Bridge — 收敛测试 (对齐采样)")
    print("=" * 60)

    results = {}
    for dt in dts:
        results[dt] = integrate(dt)
        print(f"  dt={int(1/dt):>4d}d  done")

    print()
    print(f'  {"dt":>8s}  {"RMS vs finest":>14s}  {"ratio":>8s}  {"order":>6s}')
    print(f'  {"-" * 40}')
    for i in range(len(dts) - 1):
        e1 = np.linalg.norm(results[dts[i]] - results[dts[-1]])
        e1 /= np.sqrt(len(results[dts[-1]]))
        e2 = np.linalg.norm(results[dts[i + 1]] - results[dts[-1]])
        e2 /= np.sqrt(len(results[dts[-1]]))
        if e2 > 0:
            r = e1 / e2
            o = np.log2(r)
            print(f"  {1/dts[i]:>4.0f}d -> {1/dts[i+1]:>4.0f}d"
                  f"  {e1:>14.3e}  {r:>8.3f}  {o:>6.3f}")

    print()
    print(f"  理论预期二阶: ratio ≈ 4.0")

    # ── 与原始 KDK 对比 ──
    print()
    print("=" * 60)
    print("对比: Bridge vs 原始 KDK (对齐采样, dt=1d)")
    print("=" * 60)

    b = integrate(1 / 365)
    o = integrate_original(1 / 365)
    diff = np.linalg.norm(b - o) / np.sqrt(len(b))
    print(f"  Bridge L2 norm           = {np.linalg.norm(b)/np.sqrt(len(b)):.6e}")
    print(f"  Original L2 norm         = {np.linalg.norm(o)/np.sqrt(len(o)):.6e}")
    print(f"  L2 difference            = {diff:.3e}")
    print(f"  (两者应非常接近，因 EMB 反作用极小)")

    # ── 能量诊断 ──
    print()
    print("=" * 60)
    print("Energy diagnostics (Bridge, dt=1d)")
    print("=" * 60)

    t, d, emb, E = integrate(1 / 365, return_diagnostics=True)
    dE = E - E[0]
    print(f"  dE_max  = {np.max(np.abs(dE)):.3e}")
    print(f"  dE_rms  = {np.std(dE):.3e}")
    print(f"  d_mean  = {np.mean(d):.6f} AU")
    print(f"  d_std   = {np.std(d):.3e} AU")
