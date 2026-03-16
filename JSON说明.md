# 太阳系数据文件说明 (solar_system.json)

## 文件结构概述

存储了太阳系八大行星的轨道参数和质量数据，用于rebound天体力学模拟。

## 数据结构

```json
{
  "star": { ... },      // 太阳数据
  "planets": [ ... ]    // 行星数据数组
}
```

---

## 参数说明

### 恒星 (star) 参数

| 参数名 | 类型 | 单位 | 说明 |
|--------|------|------|------|
| name | string | - | 天体名称（"Sun"） |
| mass | float | 太阳质量 | 质量，太阳为1.0 |
| radius | float | AU | 半径（约0.00465 AU） |
| x, y, z | float | AU | 位置坐标，太阳位于原点(0,0,0) |
| vx, vy, vz | float | AU/yr | 速度分量，太阳速度为(0,0,0) |

---

### 行星 (planets) 参数

每个行星对象包含以下轨道要素：

| 参数名 | 类型 | 单位 | 说明 |
|--------|------|------|------|
| name | string | - | 行星名称（英文） |
| mass | float | 太阳质量 | 行星质量（相对于太阳） |
| a | float | AU | 半长轴（semimajor axis） |
| e | float | - | 离心率（eccentricity），0≤e<1 |
| inc | float | 度(°) | 轨道倾角（inclination） |
| Omega | float | 度(°) | 升交点经度（longitude of ascending node） |
| omega | float | 度(°) | 近心点幅角（argument of perihelion） |
| M | float | 度(°) | 平近点角（mean anomaly） |

---

## 八大行星详细数据

### 1. 水星 (Mercury)
- **质量**: 1.660136×10⁻⁷ 太阳质量
- **轨道半径**: 0.387 AU
- **特点**: 离心率最大(0.206)，轨道最扁

### 2. 金星 (Venus)
- **质量**: 2.447838×10⁻⁶ 太阳质量
- **轨道半径**: 0.723 AU
- **特点**: 轨道最接近圆形(离心率0.00677)

### 3. 地球 (Earth)
- **质量**: 3.003489×10⁻⁶ 太阳质量
- **轨道半径**: 1.000 AU（参考标准）
- **特点**: 轨道倾角接近0°

### 4. 火星 (Mars)
- **质量**: 3.227151×10⁻⁷ 太阳质量
- **轨道半径**: 1.524 AU
- **特点**: 离心率较大(0.093)

### 5. 木星 (Jupiter)
- **质量**: 9.547919×10⁻⁴ 太阳质量（最重）
- **轨道半径**: 5.204 AU
- **特点**: 质量最大，约为其他行星总和的2.5倍

### 6. 土星 (Saturn)
- **质量**: 2.858859×10⁻⁴ 太阳质量
- **轨道半径**: 9.582 AU
- **特点**: 著名的环系统

### 7. 天王星 (Uranus)
- **质量**: 4.366244×10⁻⁵ 太阳质量
- **轨道半径**: 19.229 AU
- **特点**: 轨道倾角较小(0.773°)

### 8. 海王星 (Neptune)
- **质量**: 5.151389×10⁻⁵ 太阳质量
- **轨道半径**: 30.104 AU
- **特点**: 离心率很小(0.00946)，轨道接近圆形

---

## 使用示例

### Python读取JSON数据

```python
import json
import rebound

# 读取JSON文件
with open('solar_system.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

# 创建rebound模拟
sim = rebound.Simulation()
sim.units = ('yr', 'AU', 'Msun')

# 添加太阳
sun = data['star']
sim.add(
    name=sun['name'],
    mass=sun['mass'],
    x=sun['x'], y=sun['y'], z=sun['z'],
    vx=sun['vx'], vy=sun['vy'], vz=sun['vz']
)

# 添加行星
for planet in data['planets']:
    sim.add(
        name=planet['name'],
        mass=planet['mass'],
        a=planet['a'],
        e=planet['e'],
        inc=planet['inc'],
        Omega=planet['Omega'],
        omega=planet['omega'],
        M=planet['M']
    )

# 移动到质心系
sim.move_to_com()

# 进行模拟
sim.integrate(100)  # 模拟100年
```

---

## 单位说明

- **距离**: 天文单位 (AU)
  - 1 AU = 地球到太阳的平均距离 ≈ 1.496×10⁸ km

- **时间**: 年 (yr)
  - 1 yr = 地球公转周期 ≈ 365.25 天

- **质量**: 太阳质量 (Msun)
  - 1 Msun = 1.989×10³⁰ kg

- **角度**: 度 (°)

---

## 注意事项

1. **角度单位**: 所有角度参数使用度数，非弧度
2. **质量比例**: 质量以太阳质量为单位，都是小数值
3. **轨道参数**: 使用开普勒轨道要素，rebound会自动转换为笛卡尔坐标
4. **参考历元**: 这些轨道要素基于J2000历元
5. **精度**: 数据精度对小行星模拟非常重要

---

## 扩展说明

如需添加更多天体（如小行星、卫星），请遵循相同的参数格式：

```json
{
  "name": "天体名称",
  "mass": 质量值,
  "a": 半长轴,
  "e": 离心率,
  "inc": 轨道倾角,
  "Omega": 升交点经度,
  "omega": 近心点幅角,
  "M": 平近点角
}
```
