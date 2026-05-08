## 辛桥 (Symplectic Bridge) — 完整技术方案

### 1. 问题：多时间尺度

太阳系中卫星轨道周期比行星短 2–3 个数量级：

| 系统 | 行星周期 | 卫星周期 | 步长需求 |
|------|----------|----------|----------|
| 地月 | 1 年 | 27 天 | 统一 dt ≈ 1 天可接受 |
| 木卫 | 12 年 | 1.8 天 | 若统一 dt，计算量 x200 |

若对全太阳系采用统一小步长来分辨木卫一，八大行星的积分成本将增加两个数量级。必须做**多速率积分**：主系统大步长，卫星系小步长。

### 2. 哈密顿量分裂

考虑日-行星-卫星三层系统，总哈密顿量为：

$$H = H_{\text{main}}(Q,P) + H_{\text{sub}}(q,p) + H_{\text{cross}}(q,Q)$$

其中：

- **$H_{\text{main}}$**：日-行星二体（大步长 WHFast 积分）
  $$H_{\text{main}} = \frac{P^2}{2M_{\text{pl}}} - \frac{GM_{\text{sun}}M_{\text{pl}}}{|\mathbf{Q}|}$$

- **$H_{\text{sub}}$**：行星-卫星多体（小步长 WHFast 积分）
  $$H_{\text{sub}} = \sum_i \frac{p_i^2}{2m_i} - \sum_i \frac{GM_{\text{pl}}m_i}{|\mathbf{q}_i|} + \sum_{i<j} \frac{G m_i m_j}{|\mathbf{q}_i - \mathbf{q}_j|}$$

- **$H_{\text{cross}}$**：太阳潮汐力（耦合项，大步长 Kick）
  $$H_{\text{cross}} = \sum_i \left[ \frac{GM_{\text{sun}} m_i}{|\mathbf{Q} + \mathbf{q}_i|} - \frac{GM_{\text{sun}} m_i}{|\mathbf{Q}|} - \frac{GM_{\text{sun}} m_i \cdot \mathbf{q}_i \cdot \mathbf{Q}}{|\mathbf{Q}|^3} \right]$$

  （括号内正是潮汐势的标准展开：直接引力减质心引力减常数项）

### 3. 算子分裂：二阶 Strang

$$e^{\Delta t(\hat{H}_{\text{main}} + \hat{H}_{\text{sub}} + \hat{H}_{\text{cross}})} \approx e^{\frac{\Delta t}{2}\hat{H}_{\text{cross}}} \cdot e^{\Delta t(\hat{H}_{\text{main}} + \hat{H}_{\text{sub}})} \cdot e^{\frac{\Delta t}{2}\hat{H}_{\text{cross}}}$$

**物理步骤**：

```
Step 1: 半部 Kick (Δt/2)
  计算 H_cross 的梯度（潮汐力），修改子系统和主系统的动量

Step 2: 完整 Drift (Δt)
  两个 REBOUND 模拟独立推进到同一目标时间：
    - main_sim 以大步长 Δt 积分 H_main
    - sub_sim  以小步长 δt = Δt/N 积分 H_sub

Step 3: 半部 Kick (Δt/2)
  在新位置重新计算潮汐力，再次修改动量
```

**BCH 公式误差分析**：

$$e^{\frac{\Delta t}{2}\hat{C}} e^{\Delta t(\hat{A}+\hat{B})} e^{\frac{\Delta t}{2}\hat{C}} = e^{\Delta t(\hat{A}+\hat{B}+\hat{C}) + O(\Delta t^3)}$$

其中 $\hat{A}=\hat{H}_{\text{main}}$, $\hat{B}=\hat{H}_{\text{sub}}$, $\hat{C}=\hat{H}_{\text{cross}}$。由于 Strang 分裂对称，$\Delta t$ 奇次项抵消，整体达到**二阶精度**。

### 4. 潮汐力计算

Kick 阶段需要 $H_{\text{cross}}$ 对全相空间的梯度：

$$\frac{\partial H_{\text{cross}}}{\partial \mathbf{q}_i} = m_i \left[ \mathbf{a}_{\text{sun}}(\mathbf{r}_i) - \mathbf{a}_{\text{sun}}(\mathbf{R}_{\text{bary}}) \right]$$

$$\frac{\partial H_{\text{cross}}}{\partial \mathbf{Q}} = -\sum_i \frac{\partial H_{\text{cross}}}{\partial \mathbf{q}_i} \quad \text{(动量守恒)}$$

其中：

$$\mathbf{a}_{\text{sun}}(\mathbf{r}) = -\frac{GM_{\text{sun}}}{|\mathbf{r}-\mathbf{r}_{\text{sun}}|^3}(\mathbf{r}-\mathbf{r}_{\text{sun}})$$

算法中的具体实现：

```python
def tidal_force(sub_sim, main_sim):
    sun_pos = pos(main_sim[0])          # 太阳位置
    bc_helio = pos(main_sim[1])         # 行星质心位置（日心）
    sub_bc = barycenter(sub_sim)        # 子系统质心（≈0）

    for each particle p in sub_sim:
        r_helio = bc_helio + pos(p)     # 天体日心坐标
        a_tidal = a_sun(r_helio) - a_sun(bc_helio + sub_bc)
        kick(p, a_tidal * dt/2)

    # 反作用力：回写主系统质心
    a_back = -sum(m_i * a_tidal_i) / M_total
    kick(main_sim[1], a_back * dt/2)
```

**关键性质**：

- 系统**自治**：$H_{\text{cross}}$ 仅依赖当前全相空间 $(q, Q)$，不显式依赖时间 $t$
- **动量守恒**：子系统和主系统之间通过反作用力闭环
- **时间可逆**：Kick-Drift-Kick 序列对称，二阶精度

### 5. 输出采样对齐陷阱

**曾踩过的坑**：收敛测试显示 ~1 阶而非预期 2 阶。

**原因**：不同 dt 的积分在不等价的时刻采样。若每固定步数记录一次（如每 10 步），dt 不同则采样时间也不同，输出混入了 O(dt) 的采样误差。

**修复**：对齐到全局参考时间序列：

```python
# 错误做法：每 N 步记录（采样时间随 dt 漂移）
while t < t_end:
    step()                   # 步长 = dt
    if step_count % 10 == 0:
        record()             # 不同 dt 在不同时刻记录

# 正确做法：对齐到固定时间网格
sample_times = linspace(0, t_end, n_samples)[1:]

for target_t in sample_times:
    while t < target_t:
        step(min(dt, target_t - t))   # 可变最后一步
    record()                           # 永远在精确的 target_t 记录
```

修复后收敛测试结果（地月系距离精度）：

| dt 对比 | RMS 误差比 | 收敛阶 |
|---------|-----------|--------|
| 365d → 730d | 4.3× | 2.1 |
| 730d → 1460d | 4.9× | 2.3 |
| 1460d → 2920d | 5.2× | 2.4 |

理论期望：二阶精度 → 步长减半误差降 4×。实测 4.3–5.2×，**确认二阶收敛**。

### 6. 三组验证测试

#### 测试 1：倒置日心坐标（"杀死引力"）

将主系统日心坐标中的太阳移到 $10^6$ AU 外，使潮汐力 ≈ 0。地月系应退化为孤立二体，WHFast 保持轨道不变。

- 输入：潮汐力 ≈ $10^{-6}$ 量级
- 输出：地月距离恒定，能量守恒
- 结果：通过 ✓

#### 测试 2：收敛阶测试（dt 扫描）

在 365d–2920d 范围内扫描外步长，以最细步长为参考计算 RMS 误差。确认二阶收敛。

- 结果：通过 ✓（见上表）

#### 测试 3：能量诊断

1 年积分，地月系能量涨落 $\Delta E_{\text{rms}} \approx 10^{-8}$，长期无漂移。

- 结果：通过 ✓

### 7. 木卫系扩展

架构与地月系完全同构：

| 层次 | 地月系 | 木卫系 |
|------|--------|--------|
| H_main | 太阳 + 地月质心 (a=1 AU) | 太阳 + 木星系质心 (a=5.2 AU) |
| H_sub | 地球 (p[0]) + 月球 (p[1]) | 木星 (p[0]) + 四颗伽利略卫星 |
| H_cross | 太阳潮汐力 | 太阳潮汐力 |

**易错点**：木卫质量通常以 $M_{\text{JUP}}$ 为单位，需换算为 $M_{\text{sun}}$：

```
Io 质量 = 4.70e-5 × M_JUP = 4.70e-5 × (1/1047.35) Msun ≈ 4.49e-8 Msun
```

若直接使用 4.70e-5 Msun，Io 质量偏大 1000×，WHFast 的"摄动 << 中心势"前提崩溃，卫星被弹射。

**验证结果**（50 年积分）：

| 卫星 | 周期 | 误差 |
|------|------|------|
| Io | 1.769 天 | < 0.1% |
| Europa | 3.552 天 | < 0.1% |
| Ganymede | 7.154 天 | < 0.1% |
| Callisto | 16.67 天 | < 0.1% |

谐振角 $\phi = \lambda_{\text{Io}} - 3\lambda_{\text{Europa}} + 2\lambda_{\text{Ganymede}}$ 锚定在 180° 附近做大幅天平动，dφ/dt 仅 0.05 rad/yr，与 µ 值完全自洽。

### 8. 与原始 KDK 对比

原始 KDK（单向 Kick，仅子系受力，无反馈）：

$$\text{Kick}_{\text{sub}} \to \text{Drift}_{\text{sub}} \to \text{Drift}_{\text{main}} \to \text{Kick}_{\text{sub}}$$

问题：动量不守恒，主系不感知子系。

辛桥（双向 Kick，动量闭环）：

$$\text{Kick}_{\text{sub+main}} \to \text{Drift}_{\text{sub}} + \text{Drift}_{\text{main}} \to \text{Kick}_{\text{sub+main}}$$

对比测试（dt = 1 天，1 年积分）：

| 指标 | 原始 KDK | 辛桥 |
|------|---------|------|
| 地月距离 RMS | 0.002571 AU | 0.002571 AU |
| 二者差异 | 2.2e-10 AU | — |
| 动量守恒 | ✗ (漂移) | ✓ (闭环) |

二者在短期模拟中结果几乎相同（地月系 EMB 反作用极小），但辛桥保证了长期动量守恒和可扩展性。
