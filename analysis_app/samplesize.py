"""样本量探索：基于先导数据方差与用户给定真实偏移的功效曲线。

这里只回答“若真实偏移为 δ，计划规模能多大概率通过 TOST”，
不把观测到的差异或事后功效当作等效证据。
"""

import math
from collections import defaultdict

from .schema import UNIT
from .stats import normal_cdf, student_t_ppf


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} 必须是数值')
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{name} 必须是正的有限数')
    return value


def _unit_variance(records):
    """先导数据的 sample 级差值方差。

    valid=0 不参与；重复送样先按 sample 合并，避免把重复测量当作
    独立单位而低估所需样本量。
    """
    groups = defaultdict(list)
    for row in records:
        if row['valid'] not in (0, 1):
            raise ValueError('valid 只能是 0 或 1')
        if row['valid'] == 0 or row['old_value'] is None:
            continue
        value = row['new_value'] - row['old_value']
        if not math.isfinite(value):
            raise ValueError('配对差值必须有限')
        groups[row['sample']].append(value)
    means = [math.fsum(sorted(values)) / len(values)
             for _, values in sorted(groups.items())]
    if len(means) < 2:
        raise ValueError('先导数据至少需要两个有效 sample')
    mean = math.fsum(sorted(means)) / len(means)
    variance = math.fsum((value - mean) ** 2 for value in means) / (len(means) - 1)
    if variance <= 0.0:
        raise ValueError('先导数据方差为零，无法计算样本量')
    return variance, len(means)


def power_curve(records, lower, upper, effect, alpha, sample_sizes,
                design_effect=1.0):
    """给定真实偏移 effect 的功效曲线。

    variance 是 sample 级差值方差；design_effect 允许把 lot/operator
    聚类造成的方差放大计入；检验自由度为 n - 1。
    """
    variance, pilot_units = _unit_variance(records)
    if isinstance(effect, bool) or not isinstance(effect, (int, float)) \
            or not math.isfinite(effect):
        raise ValueError('真实偏移必须是有限数')
    if isinstance(lower, bool) or isinstance(upper, bool) \
            or not isinstance(lower, (int, float)) \
            or not isinstance(upper, (int, float)) \
            or not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError('等效界值必须是有限数')
    if not lower < upper:
        raise ValueError('等效下界必须小于上界')
    if not 0.0 < alpha < 0.5:
        raise ValueError('alpha 必须在 (0, 0.5) 内')
    design_effect = _positive(design_effect, '设计效应')
    effective_variance = variance * design_effect

    points = []
    for size in sample_sizes:
        if isinstance(size, bool) or not isinstance(size, int) or size < 2:
            raise ValueError('样本量必须是不少于 2 的整数')
        standard_error = math.sqrt(effective_variance / size)
        critical = student_t_ppf(1.0 - alpha, size - 1)
        power = (normal_cdf((upper - effect) / standard_error - critical)
                 - normal_cdf((lower - effect) / standard_error + critical))
        points.append({
            'sample_size': size,
            'standard_error': standard_error,
            'critical_value': critical,
            'power': min(1.0, max(0.0, power)),
        })
    reaching = [point['sample_size'] for point in points
                if point['power'] >= 0.8]
    return {
        'unit': UNIT,
        'pilot_variance': variance,
        'pilot_units': pilot_units,
        'design_effect': design_effect,
        'effective_variance': effective_variance,
        'assumed_true_difference': float(effect),
        'bounds': {'lower': float(lower), 'upper': float(upper)},
        'alpha': alpha,
        'target_power': 0.8,
        'smallest_size_at_target': min(reaching) if reaching else None,
        'curve': points,
        'note': '功效基于假定的真实偏移，不是观测证据；不得作为等效结论。',
    }
