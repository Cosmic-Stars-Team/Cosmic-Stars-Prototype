import numpy as np
import rebound

G = 4.0 * np.pi ** 2  # (AU, yr, Msun) 单位制下，G = 4π²

# -----------------------------
# 工具函数
# -----------------------------


def compute_com(sim):
    m = 0.0
    r = np.zeros(3)
    v = np.zeros(3)
    for p in sim.particles:
        m += p.m
        r += p.m * np.array([p.x, p.y, p.z])
        v += p.m * np.array([p.vx, p.vy, p.vz])
    return r / m, v / m


def tidal_acc(body_pos, com_pos, sun_pos, M_sun=1.0, G=G):
    r1 = body_pos - sun_pos
    r2 = com_pos - sun_pos

    a1 = -G * M_sun * r1 / np.linalg.norm(r1) ** 3
    a2 = -G * M_sun * r2 / np.linalg.norm(r2) ** 3

    return a1 - a2


def apply_kick(sub_sim, emb_pos, sun_pos, dt_main):
    com_pos, _ = compute_com(sub_sim)

    for p in sub_sim.particles:
        # 子系坐标 = 地月质心系坐标，加上 EMB 日心位置得到日心坐标
        body_rel = np.array([p.x, p.y, p.z])
        body_helio = emb_pos + body_rel
        com_helio = emb_pos + com_pos
        a = tidal_acc(body_helio, com_helio, sun_pos)
        p.vx += a[0] * (dt_main * 0.5)
        p.vy += a[1] * (dt_main * 0.5)
        p.vz += a[2] * (dt_main * 0.5)


# -----------------------------
# 主系统（Sun + EMB）
# -----------------------------

main_sim = rebound.Simulation()
main_sim.units = ("AU", "yr", "Msun")

# Sun
main_sim.add(m=1.0)

# EMB（近似地球轨道）
main_sim.add(m=3e-6, a=1.0)  # 地月总质量

main_sim.move_to_com()
main_sim.dt = 1 / 365  # 1 天

# -----------------------------
# 子系统（Earth + Moon）
# -----------------------------

sub_sim = rebound.Simulation()
sub_sim.units = ("AU", "yr", "Msun")

# Earth
sub_sim.add(m=3.0e-6 * 0.987)

# Moon
sub_sim.add(m=3.0e-6 * 0.013, a=0.00257)  # ~384000 km ≈ 0.00257 AU

sub_sim.move_to_com()
sub_sim.dt = 1 / 365 / 50  # ~200 秒

# -----------------------------
# 主循环
# -----------------------------

t = 0.0
t_end = 1.0  # 1 年
dt_main = main_sim.dt

print("Start simulation...")

while t < t_end:
    # 当前太阳和 EMB 位置（日心系）
    sun = main_sim.particles[0]
    sun_pos = np.array([sun.x, sun.y, sun.z])
    emb = main_sim.particles[1]
    emb_pos = np.array([emb.x, emb.y, emb.z])

    # 1️⃣ Kick（开始）
    apply_kick(sub_sim, emb_pos, sun_pos, dt_main)

    # 2️⃣ 子系统推进
    sub_sim.integrate(t + dt_main)

    # 3️⃣ 主系统推进
    main_sim.integrate(t + dt_main)

    # 更新太阳和 EMB 位置
    sun = main_sim.particles[0]
    sun_pos_new = np.array([sun.x, sun.y, sun.z])
    emb = main_sim.particles[1]
    emb_pos_new = np.array([emb.x, emb.y, emb.z])

    # 4️⃣ Kick（结束）
    apply_kick(sub_sim, emb_pos_new, sun_pos_new, dt_main)

    t += dt_main

    # -----------------------------
    # 输出（简单检查）
    # -----------------------------
    earth = sub_sim.particles[0]
    moon = sub_sim.particles[1]

    r = np.linalg.norm(
        np.array([earth.x, earth.y, earth.z]) - np.array([moon.x, moon.y, moon.z])
    )

    print(f"t={t:.3f} yr  Earth-Moon distance={r:.6f} AU")
