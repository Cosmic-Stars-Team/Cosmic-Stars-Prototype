import math
import time
from dataclasses import dataclass
from typing import Callable, Dict

# 时间换算辅助常量
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def sim_seconds_to_years(sim_seconds: float) -> float:
    """把模拟秒数转换为模拟年数。"""
    return sim_seconds / SECONDS_PER_YEAR


def years_to_sim_seconds(years: float) -> float:
    """把模拟年数转换为模拟秒数。"""
    return years * SECONDS_PER_YEAR


def ratio_to_scale_yr_per_sec(sim_seconds_per_real_second: float) -> float:
    """
    把“1秒真实时间 -> X秒模拟时间”的比例
    转成内部统一单位“模拟年/真实秒”。
    """
    return sim_seconds_to_years(sim_seconds_per_real_second)


# 常用速度档位（内部单位：模拟年/真实秒）
SCALE_PRESETS: Dict[str, float] = {
    "1s/1ms": ratio_to_scale_yr_per_sec(0.001),
    "1s/1sec": ratio_to_scale_yr_per_sec(1.0),
    "1s/1day": ratio_to_scale_yr_per_sec(24.0 * 3600.0),
    "1s/1month": ratio_to_scale_yr_per_sec(30.0 * 24.0 * 3600.0),
    "1s/1year": ratio_to_scale_yr_per_sec(365.25 * 24.0 * 3600.0),
}


@dataclass
class TimeScaleStepResult:
    # 本次更新使用的真实时间间隔
    real_dt_sec: float
    # 本次更新平滑前后的时间倍率
    scale_before_yr_per_sec: float
    scale_after_yr_per_sec: float
    # 本次更新推进前后的“模拟时间预算”
    budget_before_yr: float
    budget_after_yr: float
    # 本次更新执行的固定子步数量
    substeps: int
    # 本次更新实际推进的模拟年数
    simulated_yr: float
    # 更新后的控制器状态
    sim_time_yr: float
    tick: int


class TimeScaleSynchronizer:
    """
    后端时间尺度同步器。

    它维护三条时间线的一致性：
    1) 真实时间（墙钟时间）
    2) 模拟时间（后端权威时间）
    3) 固定物理子步（保证积分稳定）
    """

    def __init__(
        self,
        physics_dt_yr: float = 1.0 / (365.25 * 24.0),  # 约1小时模拟时间
        initial_scale_yr_per_sec: float = SCALE_PRESETS["1s/1sec"],
        max_substeps: int = 256,
    ) -> None:
        if physics_dt_yr <= 0.0:
            raise ValueError("physics_dt_yr 必须大于 0")
        if initial_scale_yr_per_sec < 0.0:
            raise ValueError("initial_scale_yr_per_sec 必须大于等于 0")
        if max_substeps <= 0:
            raise ValueError("max_substeps 必须大于 0")

        # 固定物理步长（单位：模拟年），不要每帧改变
        self.physics_dt_yr = physics_dt_yr
        # 当前与目标播放速度（单位：模拟年/真实秒）
        self.current_scale_yr_per_sec = initial_scale_yr_per_sec
        self.target_scale_yr_per_sec = initial_scale_yr_per_sec
        # 速度变化的平滑时间常数
        self.ramp_tau_sec = 0.0
        # 待积分的模拟时间预算
        self._budget_yr = 0.0
        # 权威模拟时钟与逻辑tick
        self.sim_time_yr = 0.0
        self.tick = 0
        # 安全上限：防止单帧追赶步数过多
        self.max_substeps = max_substeps

    def set_time_scale(self, target_yr_per_sec: float, ramp_sec: float = 0.3) -> None:
        """
        请求新的时间倍率。
        ramp_sec 控制加减速是否平滑。
        """
        if target_yr_per_sec < 0.0:
            raise ValueError("target_yr_per_sec 必须大于等于 0")
        if ramp_sec < 0.0:
            raise ValueError("ramp_sec 必须大于等于 0")

        self.target_scale_yr_per_sec = target_yr_per_sec
        self.ramp_tau_sec = ramp_sec

        # 如果禁用平滑，立即切换到目标倍率
        if ramp_sec == 0.0:
            self.current_scale_yr_per_sec = target_yr_per_sec

    def set_time_scale_by_ratio(
        self, sim_seconds_per_real_second: float, ramp_sec: float = 0.3
    ) -> None:
        """辅助方法：按“1秒真实时间 -> X秒模拟时间”设置倍率。"""
        self.set_time_scale(
            ratio_to_scale_yr_per_sec(sim_seconds_per_real_second), ramp_sec=ramp_sec
        )

    def set_time_scale_by_preset(self, preset: str, ramp_sec: float = 0.3) -> None:
        """辅助方法：使用预设档位（如 `1s/1ms`、`1s/1year`）。"""
        if preset not in SCALE_PRESETS:
            raise KeyError(f"未知预设 '{preset}'，可选：{list(SCALE_PRESETS)}")
        self.set_time_scale(SCALE_PRESETS[preset], ramp_sec=ramp_sec)

    def _smooth_scale(self, real_dt_sec: float) -> None:
        """使用指数平滑把当前倍率推进到目标倍率。"""
        if self.current_scale_yr_per_sec == self.target_scale_yr_per_sec:
            return
        if self.ramp_tau_sec <= 0.0:
            self.current_scale_yr_per_sec = self.target_scale_yr_per_sec
            return

        # alpha 在 (0,1)：dt 越大，本帧收敛越快
        alpha = 1.0 - math.exp(-real_dt_sec / self.ramp_tau_sec)
        self.current_scale_yr_per_sec += (
            self.target_scale_yr_per_sec - self.current_scale_yr_per_sec
        ) * alpha

    def advance(
        self,
        real_dt_sec: float,
        step_simulation: Callable[[float], None],
    ) -> TimeScaleStepResult:
        """
        按真实时间 real_dt_sec 推进一步同步器。

        step_simulation 回调会收到固定模拟步长（单位：年）。
        REBOUND 示例：
            lambda dt_yr: sim.integrate(sim.t + dt_yr)
        """
        if real_dt_sec < 0.0:
            raise ValueError("real_dt_sec 必须大于等于 0")

        scale_before = self.current_scale_yr_per_sec
        budget_before = self._budget_yr

        # 1) 先平滑更新时间倍率
        self._smooth_scale(real_dt_sec)

        # 2) 把真实时间换算成“可推进的模拟时间预算”
        self._budget_yr += real_dt_sec * self.current_scale_yr_per_sec

        # 3) 限制预算上限，避免帧卡顿后发生“死亡螺旋”
        max_budget_yr = self.physics_dt_yr * self.max_substeps
        if self._budget_yr > max_budget_yr:
            self._budget_yr = max_budget_yr

        # 4) 用固定步长循环积分，直到预算耗尽
        substeps = 0
        simulated_yr = 0.0
        while self._budget_yr >= self.physics_dt_yr and substeps < self.max_substeps:
            step_simulation(self.physics_dt_yr)
            self._budget_yr -= self.physics_dt_yr
            self.sim_time_yr += self.physics_dt_yr
            self.tick += 1
            simulated_yr += self.physics_dt_yr
            substeps += 1

        # 5) 返回本次推进统计，便于日志和监控
        return TimeScaleStepResult(
            real_dt_sec=real_dt_sec,
            scale_before_yr_per_sec=scale_before,
            scale_after_yr_per_sec=self.current_scale_yr_per_sec,
            budget_before_yr=budget_before,
            budget_after_yr=self._budget_yr,
            substeps=substeps,
            simulated_yr=simulated_yr,
            sim_time_yr=self.sim_time_yr,
            tick=self.tick,
        )

    def build_time_sync_payload(self) -> Dict[str, float]:
        """
        生成给前端的最小时间同步数据包。
        前端可用这些字段做插值和时钟对齐。
        """
        return {
            "tick": self.tick,
            "sim_time_yr": self.sim_time_yr,
            "time_scale_yr_per_real_sec": self.current_scale_yr_per_sec,
            "time_scale_target_yr_per_real_sec": self.target_scale_yr_per_sec,
            "server_unix_time_sec": time.time(),
        }


def make_rebound_stepper(sim) -> Callable[[float], None]:
    """
    REBOUND 适配器，返回固定步长积分回调。
    用法：
        stepper = make_rebound_stepper(sim)
        controller.advance(real_dt_sec=0.016, step_simulation=stepper)
    """

    def _step(dt_yr: float) -> None:
        sim.integrate(sim.t + dt_yr)

    return _step
