"""确定性统计工具。

只依赖标准库。Student t 分布用正则化不完全贝塔函数计算，避免引入
第三方运行依赖；所有函数对相同输入返回相同浮点结果。
"""

import math


_MAX_ITER = 300
_EPS = 3.0e-14
_FPMIN = 1.0e-300


def _betacf(a, b, x):
    """不完全贝塔函数的连分式（Numerical Recipes 形式）。"""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _FPMIN:
        d = _FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, _MAX_ITER + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = 1.0 + aa / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < _FPMIN:
            d = _FPMIN
        c = 1.0 + aa / c
        if abs(c) < _FPMIN:
            c = _FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return h


def _betai(a, b, x):
    """正则化不完全贝塔函数 I_x(a, b)。"""
    if not (a > 0.0 and b > 0.0):
        raise ValueError('贝塔函数参数必须为正')
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def student_t_cdf(value, degrees_of_freedom):
    """Student t 分布的累积分布函数。"""
    if degrees_of_freedom <= 0:
        raise ValueError('自由度必须为正')
    if not math.isfinite(value):
        if value > 0:
            return 1.0
        if value < 0:
            return 0.0
        raise ValueError('t 值必须有限')
    tail = 0.5 * _betai(degrees_of_freedom / 2.0, 0.5,
                        degrees_of_freedom / (degrees_of_freedom + value * value))
    return 1.0 - tail if value >= 0 else tail


def student_t_ppf(probability, degrees_of_freedom):
    """Student t 分布的分位数，二分法保证确定性。"""
    if not 0.0 < probability < 1.0:
        raise ValueError('概率必须在 (0, 1) 内')
    if degrees_of_freedom <= 0:
        raise ValueError('自由度必须为正')
    if probability == 0.5:
        return 0.0
    if probability < 0.5:
        return -student_t_ppf(1.0 - probability, degrees_of_freedom)

    low = 0.0
    high = 1.0
    while student_t_cdf(high, degrees_of_freedom) < probability:
        high *= 2.0
        if high > 1.0e12:
            raise ValueError('t 分位数超出可计算范围')
    for _ in range(200):
        middle = (low + high) / 2.0
        if student_t_cdf(middle, degrees_of_freedom) < probability:
            low = middle
        else:
            high = middle
        if high - low <= 1.0e-13 * max(1.0, high):
            break
    return (low + high) / 2.0


def normal_cdf(value):
    """标准正态分布累积分布函数。"""
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def tost(estimate, standard_error, degrees_of_freedom, lower, upper, alpha):
    """配对 TOST：两个单侧检验与对应的 1-2α 双侧区间。

    下界检验 H0: μ <= lower，上界检验 H0: μ >= upper。
    两个单侧 p 值取较大者作为 TOST 的 p 值；区间使用
    t_{1-α, df} 的半宽，即置信水平 1-2α。
    """
    if not (math.isfinite(estimate) and math.isfinite(standard_error)):
        raise ValueError('估计值与标准误必须有限')
    if standard_error <= 0.0:
        raise ValueError('标准误必须为正')
    if not lower < upper:
        raise ValueError('等效上下界顺序错误')
    if not 0.0 < alpha < 0.5:
        raise ValueError('alpha 必须在 (0, 0.5) 内')

    lower_statistic = (estimate - lower) / standard_error
    upper_statistic = (estimate - upper) / standard_error
    p_lower = 1.0 - student_t_cdf(lower_statistic, degrees_of_freedom)
    p_upper = student_t_cdf(upper_statistic, degrees_of_freedom)
    p_value = max(p_lower, p_upper)
    critical = student_t_ppf(1.0 - alpha, degrees_of_freedom)
    interval_lower = estimate - critical * standard_error
    interval_upper = estimate + critical * standard_error
    return {
        'p_lower': min(1.0, max(0.0, p_lower)),
        'p_upper': min(1.0, max(0.0, p_upper)),
        'p_value': min(1.0, max(0.0, p_value)),
        'critical_value': critical,
        'interval': [interval_lower, interval_upper],
        'interval_contains_zero': interval_lower <= 0.0 <= interval_upper,
        'interval_within_bounds': (
            interval_lower >= lower and interval_upper <= upper
        ),
    }


def holm_adjust(p_values, alpha):
    """Holm 逐步法；None 表示该假设不可估计，仍占一个家族位置。"""
    if not 0.0 < alpha < 0.5:
        raise ValueError('alpha 必须在 (0, 0.5) 内')
    total = len(p_values)
    if total == 0:
        return {'adjusted': [], 'reject': []}
    numeric = [1.0 if value is None else value for value in p_values]
    for value in numeric:
        if not 0.0 <= value <= 1.0:
            raise ValueError('p 值必须在 [0, 1] 内')
    order = sorted(range(total), key=lambda index: (numeric[index], index))
    adjusted = [None] * total
    running = 0.0
    for position, index in enumerate(order):
        candidate = min(1.0, (total - position) * numeric[index])
        running = max(running, candidate)
        adjusted[index] = running
    reject = [
        p_values[index] is not None and adjusted[index] <= alpha
        for index in range(total)
    ]
    return {'adjusted': adjusted, 'reject': reject}
