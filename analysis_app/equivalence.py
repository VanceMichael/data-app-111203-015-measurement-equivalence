"""等效性研究：配对 TOST 与聚类稳健推断。

分析单位是 sample 的配对差值均值，重复送样不会重复计数。lot 与
operator 造成的相关性通过多向聚类稳健方差（CR1）同时处理：
先按 sample 合并重复测量，再把每个观测的残差影响按 lot、operator
及其交集聚类求和。

计划内的总体与全部预设产线共享一个 Holm 家族；只有家族内全部
判定通过且研究范围满足，才声明等效。
"""

import math
from collections import defaultdict
from itertools import combinations

from .io import fingerprint
from .plans import verify
from .schema import UNIT
from .stats import holm_adjust, tost


VERDICT_EQUIVALENT = 'equivalent'
VERDICT_NOT_PROVEN = 'not_proven'
VERDICT_UNAVAILABLE = 'unavailable'
VERDICTS = (VERDICT_EQUIVALENT, VERDICT_NOT_PROVEN, VERDICT_UNAVAILABLE)


def _usable_pairs(records):
    """valid=0 与缺失旧值都不进入配对推断。"""
    pairs = []
    for row in records:
        if row['valid'] not in (0, 1):
            raise ValueError('valid 只能是 0 或 1')
        if row['valid'] == 0 or row['old_value'] is None:
            continue
        value = row['new_value'] - row['old_value']
        if not math.isfinite(value):
            raise ValueError('配对差值必须有限')
        pairs.append(row)
    return pairs


def _sample_units(pairs):
    """把重复送样折叠成 sample 级分析单位。

    每个 sample 的差值均值是分析单位；同时记录观测数、单一 lot/line
    标签，以及是否跨 lot/line。operator 可在 sample 内交叉，由聚类
    方差处理；lot/line 必须嵌套在 sample 内，否则无法定义聚类。
    """
    by_sample = defaultdict(list)
    for row in pairs:
        by_sample[row['sample']].append(row)
    units = []
    for sample in sorted(by_sample):
        members = by_sample[sample]
        differences = [row['new_value'] - row['old_value']
                       for row in members]
        lots = sorted({row['lot'] for row in members})
        lines = sorted({row['line'] for row in members})
        units.append({
            'sample': sample,
            'members': members,
            'differences': differences,
            'mean': math.fsum(sorted(differences)) / len(differences),
            'lot': lots[0] if len(lots) == 1 else None,
            'line': lines[0] if len(lines) == 1 else None,
        })
    return units


def _operator_count(units):
    return len({row['operator'] for unit in units for row in unit['members']})


def _unavailable(units):
    """返回不可估计原因；可估计时返回 None。

    单一批次不在此处拒绝：批次方差无法估计时该维度会被省略并
    在结果中标注，研究范围判定负责阻止跨批声明。跨线 sample 仍可
    进入总体（line 不是聚类维度），但不会进入任何产线亚组，并由
    范围核对记录。
    """
    if len(units) < 2:
        return '分析单位少于两个'
    if len({unit['mean'] for unit in units}) == 1:
        return '配对差值方差为零'
    if any(unit['lot'] is None for unit in units):
        return 'sample 跨多个 lot，重复测量未嵌套在批次内'
    return None


def _cluster_counts(units, dimensions):
    counts = {}
    for dimension in dimensions:
        if dimension == 'lot':
            counts[dimension] = len({unit['lot'] for unit in units})
        elif dimension == 'operator':
            counts[dimension] = _operator_count(units)
        else:
            raise ValueError(f'未知聚类维度：{dimension}')
    return counts


def _cluster_variance(units, dimensions):
    """多向聚类稳健方差（CR1）。

    每个观测的残差影响 (x - μ) / (单位数 * 该 sample 的观测数) 使
    每个 sample 的总影响恰为 (y_s - μ) / 单位数。聚类按维度组合求和，
    每个组合用 G/(G-1) 修正，交集项以交替符号相减。
    """
    total = len(units)
    estimate = math.fsum(sorted(unit['mean'] for unit in units)) / total
    sample_lot = {unit['sample']: unit['lot'] for unit in units}
    influence = []
    for unit in units:
        weight = 1.0 / (total * len(unit['differences']))
        for row, difference in zip(unit['members'], unit['differences']):
            influence.append((row['observation_id'], row,
                              sample_lot[row['sample']],
                              (difference - estimate) * weight))
    influence.sort(key=lambda item: item[0])

    for dimension, count in _cluster_counts(units, dimensions).items():
        if count < 2:
            raise ValueError(f'{dimension} 只有一个取值，无法估计方差')

    variance = 0.0
    for size in range(1, len(dimensions) + 1):
        for combo in combinations(dimensions, size):
            sums = defaultdict(float)
            for _, row, lot, value in influence:
                key = tuple(row['operator'] if field == 'operator' else lot
                            for field in combo)
                sums[key] += value
            count = len(sums)
            if count < 2:
                raise ValueError('聚类数量不足，无法估计方差')
            pieces = [sums[key] for key in sorted(sums, key=str)]
            variance += (-1.0) ** (size + 1) * count / (count - 1.0) \
                * math.fsum(value * value for value in pieces)
    if variance < 0.0 and variance > -1.0e-15:
        variance = 0.0
    return estimate, variance


def _estimate(units, plan):
    """执行 TOST；返回 (结果, 不可估计原因)。

    聚类维度中若某维度只有一个取值（例如只有一批材料），该维度
    的方差分量不可估计，会被省略并在结果中列出，而不是当作零。
    """
    reason = _unavailable(units)
    if reason is not None:
        return None, reason
    requested = plan['clustering']['dimensions']
    counts = _cluster_counts(units, requested)
    usable = [dimension for dimension in requested if counts[dimension] >= 2]
    omitted = [dimension for dimension in requested if counts[dimension] < 2]
    if not usable:
        return None, '聚类维度均只有一个取值，无法估计方差'
    try:
        estimate, variance = _cluster_variance(units, usable)
        degrees_of_freedom = min(
            _cluster_counts(units, usable).values()) - 1
    except ValueError as error:
        return None, str(error)
    if not math.isfinite(variance) or variance < 0.0:
        return None, '聚类方差为负或非有限'
    standard_error = math.sqrt(variance)
    if standard_error == 0.0:
        return None, '聚类标准误为零'
    bounds = plan['bounds']
    result = tost(estimate, standard_error, degrees_of_freedom,
                  bounds['lower'], bounds['upper'], plan['alpha'])
    result['estimate'] = estimate
    result['standard_error'] = standard_error
    result['degrees_of_freedom'] = degrees_of_freedom
    result['variance'] = variance
    result['clusters'] = {dimension: counts[dimension] for dimension in usable}
    result['omitted_clusters'] = omitted
    return result, None


def _scope(pairs, plan):
    """研究范围核对：计划产线、跨批适用与 sample 归属。

    范围以进入配对推断的数据为准；valid=0 的行不构成证据，也不会
    被用来宣称或阻止等效。计划产线与数据中出现但未登记的产线都会
    被核对，缺失的计划产线会使对应判定不可估计。
    """
    planned = set(plan['subgroups']['lines'])
    observed = {row['line'] for row in pairs}
    lots = {row['lot'] for row in pairs}
    lots_per_sample = defaultdict(set)
    lines_per_sample = defaultdict(set)
    for row in pairs:
        lots_per_sample[row['sample']].add(row['lot'])
        lines_per_sample[row['sample']].add(row['line'])
    issues = {
        'unknown_lines': sorted(observed - planned),
        'missing_planned_lines': sorted(planned - observed),
        'single_lot': len(lots) < 2,
        'lots': len(lots),
        'samples_spanning_lots': sorted(
            sample for sample, values in lots_per_sample.items() if len(values) > 1),
        'samples_spanning_lines': sorted(
            sample for sample, values in lines_per_sample.items() if len(values) > 1),
    }
    ok = (not issues['unknown_lines']
          and not issues['missing_planned_lines']
          and not issues['samples_spanning_lots']
          and not issues['samples_spanning_lines'])
    if plan['scope']['require_multiple_lots'] and issues['single_lot']:
        ok = False
    return ok, issues


def _comparison(target, result, reason, index, family, plan, scope_issues):
    """组装单个判定的输出；未通过范围核对时不宣称等效。"""
    group_units = target['units']
    comparison = {
        'group': target['group'],
        'scope': target['scope'],
        'value': target['value'],
        'status': 'estimated',
        'verdict': VERDICT_NOT_PROVEN,
        'reason': reason,
        'estimate': None,
        'standard_error': None,
        'degrees_of_freedom': None,
        'p_value': None,
        'p_adjusted': None,
        'p_lower': None,
        'p_upper': None,
        'interval_1_minus_2alpha': None,
        'within_bounds': None,
        'units': len(group_units),
        'observations': sum(len(unit['differences']) for unit in group_units),
        'lots': len({unit['lot'] for unit in group_units if unit['lot']}),
        'operators': _operator_count(group_units),
        'repeated_samples': sum(1 for unit in group_units
                                if len(unit['differences']) > 1),
        'maximum_repeats': max((len(unit['differences']) for unit in group_units),
                               default=0),
        'clusters': {},
        'omitted_clusters': [],
        'multiple_testing': {
            'method': 'holm',
            'family_size': len(family['adjusted']),
            'alpha': plan['alpha'],
            'p_adjusted': None,
            'reject': False,
        },
    }
    if result is None:
        comparison['status'] = 'unavailable'
        comparison['verdict'] = VERDICT_UNAVAILABLE
        reason_text = reason or '不可估计'
        comparison['reason'] = reason_text
        return comparison

    comparison.update({
        'estimate': result['estimate'],
        'standard_error': result['standard_error'],
        'degrees_of_freedom': result['degrees_of_freedom'],
        'p_value': result['p_value'],
        'p_adjusted': family['adjusted'][index],
        'p_lower': result['p_lower'],
        'p_upper': result['p_upper'],
        'interval_1_minus_2alpha': result['interval'],
        'within_bounds': result['interval_within_bounds'],
        'clusters': result['clusters'],
        'omitted_clusters': result['omitted_clusters'],
    })
    comparison['multiple_testing']['p_adjusted'] = family['adjusted'][index]
    comparison['multiple_testing']['reject'] = family['reject'][index]
    if plan['scope']['require_multiple_lots'] and comparison['lots'] < 2:
        comparison['reason'] = '只有一批材料，不能跨批外推'
    elif comparison['scope'] is None and scope_issues['samples_spanning_lines']:
        comparison['reason'] = 'sample 跨多条产线'
    elif family['reject'][index] and result['interval_within_bounds']:
        comparison['verdict'] = VERDICT_EQUIVALENT
    return comparison


def analyze(records, plan):
    """按固化计划分析观测，返回判定、区间与核对信息。"""
    plan = verify(plan)
    pairs = _usable_pairs(records)
    units = _sample_units(pairs)
    scope_ok, scope_issues = _scope(pairs, plan)

    targets = [{'group': 'overall', 'scope': None, 'value': None, 'units': units}]
    for line in plan['subgroups']['lines']:
        members = [row for row in pairs if row['line'] == line]
        targets.append({'group': f'line:{line}', 'scope': 'line',
                        'value': line, 'units': _sample_units(members)})

    estimates = [(target, *_estimate(target['units'], plan))
                 for target in targets]
    p_values = [None if result is None else result['p_value']
                for _, result, _ in estimates]
    family = holm_adjust(p_values, plan['alpha'])
    all_rejected = bool(family['reject']) and all(family['reject'])
    all_unavailable = all(result is None for _, result, _ in estimates)

    comparisons = [_comparison(target, result, reason, index, family, plan,
                               scope_issues)
                   for index, (target, result, reason) in enumerate(estimates)]

    counts = {name: 0 for name in VERDICTS}
    for comparison in comparisons:
        counts[comparison['verdict']] += 1

    if not pairs or all_unavailable:
        study_verdict = VERDICT_UNAVAILABLE
    elif scope_ok and all(comparison['verdict'] == VERDICT_EQUIVALENT
                          for comparison in comparisons):
        study_verdict = VERDICT_EQUIVALENT
    else:
        study_verdict = VERDICT_NOT_PROVEN

    return {
        'product': plan['product'],
        'study_type': plan['study_type'],
        'unit': UNIT,
        'input_fingerprint': fingerprint(records),
        'plan_fingerprint': plan['plan_fingerprint'],
        'bounds': dict(plan['bounds']),
        'alpha': plan['alpha'],
        'clustering': {
            'unit': 'sample',
            'dimensions': list(plan['clustering']['dimensions']),
            'method': plan['clustering']['method'],
            'independent_unit_note': '重复送样按 sample 合并，不作为独立单位',
        },
        'multiple_testing': {
            'method': 'holm',
            'family': [comparison['group'] for comparison in comparisons],
            'alpha': plan['alpha'],
            'all_rejected': all_rejected,
        },
        'valid_pairs': len(pairs),
        'excluded_rows': len(records) - len(pairs),
        'samples': len(units),
        'lots': scope_issues['lots'],
        'operators': len({row['operator'] for row in pairs}),
        'cross_lot': not scope_issues['single_lot'],
        'scope_ok': scope_ok,
        'scope_issues': scope_issues,
        'verdict': study_verdict,
        'verdict_counts': counts,
        'comparisons': comparisons,
        'definitions': {
            'equivalent': '总体与全部预设产线在 Holm 家族内均拒绝两个单侧原假设，'
                          '1-2α 区间落在界值内，且研究范围满足',
            'not_proven': '至少一个计划判定未通过，或只有一批材料等范围前提不满足',
            'unavailable': '有效配对不足、方差为零或聚类不足，不能给出判定',
            'interval': '与配对 TOST 对应的双侧 1-2α 区间',
        },
    }
