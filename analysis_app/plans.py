"""等效性研究计划。

计划在分析前固化：界值以原始微米单位登记，分析只能读取带指纹的计划。
任何字段变化都会改变指纹，从而使旧分析不可被覆盖。
"""

import hashlib
import json
import math

from .schema import PRODUCT, UNIT


SCHEMA_VERSION = 1
STUDY_TYPE = 'paired_tost_equivalence'
DIMENSIONS = ('lot', 'operator')
PLAN_FIELDS = {
    'schema_version', 'product', 'study_type', 'bounds', 'alpha',
    'clustering', 'subgroups', 'scope',
}
CONFIG_FIELDS = {
    'lower', 'upper', 'alpha', 'lines', 'cluster_dimensions',
    'require_multiple_lots',
}


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'计划字段必须是数值：{name}')
    converted = float(value)
    if not math.isfinite(converted):
        raise ValueError(f'计划字段必须是有限数：{name}')
    return converted


def _boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f'计划字段必须是布尔值：{name}')
    return value


def _string_list(value, name):
    if not isinstance(value, list) or not value:
        raise ValueError(f'计划字段必须是非空字符串列表：{name}')
    result = []
    for item in value:
        if not isinstance(item, str) or item == '':
            raise ValueError(f'计划字段必须是非空字符串列表：{name}')
        if item not in result:
            result.append(item)
    return sorted(result)


def canonical(payload):
    """计划指纹使用的规范 JSON；键排序、无多余空白。"""
    return json.dumps(payload, ensure_ascii=False, sort_keys=True,
                      separators=(',', ':'), allow_nan=False)


def fingerprint(plan):
    content = {key: value for key, value in plan.items()
               if key != 'plan_fingerprint'}
    return hashlib.sha256(canonical(content).encode('utf-8')).hexdigest()


def normalize(plan):
    """校验计划内容并返回规范化副本，不处理指纹字段。"""
    if not isinstance(plan, dict):
        raise ValueError('计划必须是 JSON 对象')
    missing = PLAN_FIELDS - set(plan)
    extra = set(plan) - PLAN_FIELDS - {'plan_fingerprint'}
    if missing:
        raise ValueError('计划缺少字段：' + ','.join(sorted(missing)))
    if extra:
        raise ValueError('计划包含未知字段：' + ','.join(sorted(extra)))
    if plan['schema_version'] != SCHEMA_VERSION:
        raise ValueError('计划版本不受支持')
    if plan['product'] != PRODUCT:
        raise ValueError('计划产品标识不符')
    if plan['study_type'] != STUDY_TYPE:
        raise ValueError('计划研究类型不符')

    bounds = plan['bounds']
    if not isinstance(bounds, dict) or set(bounds) != {'lower', 'upper', 'unit'}:
        raise ValueError('计划界值必须包含 lower、upper、unit')
    lower = _number(bounds['lower'], 'bounds.lower')
    upper = _number(bounds['upper'], 'bounds.upper')
    if not lower < upper:
        raise ValueError('等效下界必须小于上界')
    if bounds['unit'] != UNIT:
        raise ValueError(f'界值单位必须是{UNIT}')
    alpha = _number(plan['alpha'], 'alpha')
    if not 0.0 < alpha < 0.5:
        raise ValueError('alpha 必须在 (0, 0.5) 内')

    clustering = plan['clustering']
    if not isinstance(clustering, dict) or set(clustering) != {
            'unit', 'dimensions', 'method'}:
        raise ValueError('计划聚类字段不完整')
    if clustering['unit'] != 'sample':
        raise ValueError('重复测量必须以 sample 为单位')
    if clustering['method'] != 'multiway_cluster_robust_cr1':
        raise ValueError('聚类方法必须是 multiway_cluster_robust_cr1')
    dimensions = _string_list(clustering['dimensions'], 'clustering.dimensions')
    unknown = set(dimensions) - set(DIMENSIONS)
    if unknown:
        raise ValueError('未知聚类维度：' + ','.join(sorted(unknown)))

    subgroups = plan['subgroups']
    if not isinstance(subgroups, dict) or set(subgroups) != {
            'dimension', 'lines', 'include_overall', 'multiplicity'}:
        raise ValueError('计划亚组字段不完整')
    if subgroups['dimension'] != 'line':
        raise ValueError('亚组维度必须是 line')
    lines = _string_list(subgroups['lines'], 'subgroups.lines')
    if subgroups['include_overall'] is not True:
        raise ValueError('计划必须同时包含总体判定')
    if subgroups['multiplicity'] != 'holm':
        raise ValueError('多重判定方法必须是 holm')

    scope = plan['scope']
    if not isinstance(scope, dict) or set(scope) != {'require_multiple_lots'}:
        raise ValueError('计划范围字段不完整')
    require_multiple_lots = _boolean(scope['require_multiple_lots'],
                                     'scope.require_multiple_lots')

    return {
        'schema_version': SCHEMA_VERSION,
        'product': PRODUCT,
        'study_type': STUDY_TYPE,
        'bounds': {'lower': lower, 'upper': upper, 'unit': UNIT},
        'alpha': alpha,
        'clustering': {
            'unit': 'sample',
            'dimensions': dimensions,
            'method': 'multiway_cluster_robust_cr1',
        },
        'subgroups': {
            'dimension': 'line',
            'lines': lines,
            'include_overall': True,
            'multiplicity': 'holm',
        },
        'scope': {'require_multiple_lots': require_multiple_lots},
    }


def freeze(plan):
    """返回带指纹的固化计划。"""
    normalized = normalize(plan)
    normalized['plan_fingerprint'] = fingerprint(normalized)
    return normalized


def verify(plan):
    """校验内存中的计划指纹，返回规范化副本。"""
    if not isinstance(plan, dict) or 'plan_fingerprint' not in plan:
        raise ValueError('计划缺少指纹，不能用于分析')
    claimed = plan['plan_fingerprint']
    normalized = normalize(plan)
    actual = fingerprint(normalized)
    if not isinstance(claimed, str) or claimed != actual:
        raise ValueError('计划指纹不匹配，计划在登记后被修改')
    normalized['plan_fingerprint'] = actual
    return normalized


def build_plan(lower, upper, lines, alpha=0.05,
               cluster_dimensions=None, require_multiple_lots=True):
    """从配置值构造固化计划。界值为原始微米单位。"""
    return freeze({
        'schema_version': SCHEMA_VERSION,
        'product': PRODUCT,
        'study_type': STUDY_TYPE,
        'bounds': {'lower': lower, 'upper': upper, 'unit': UNIT},
        'alpha': alpha,
        'clustering': {
            'unit': 'sample',
            'dimensions': list(cluster_dimensions or DIMENSIONS),
            'method': 'multiway_cluster_robust_cr1',
        },
        'subgroups': {
            'dimension': 'line',
            'lines': list(lines),
            'include_overall': True,
            'multiplicity': 'holm',
        },
        'scope': {'require_multiple_lots': require_multiple_lots},
    })


def from_config(config):
    """从登记配置构造计划；配置值同样使用原始微米单位。"""
    if not isinstance(config, dict):
        raise ValueError('配置必须是 JSON 对象')
    unknown = set(config) - CONFIG_FIELDS
    if unknown:
        raise ValueError('配置包含未知字段：' + ','.join(sorted(unknown)))
    missing = {'lower', 'upper', 'lines'} - set(config)
    if missing:
        raise ValueError('配置缺少字段：' + ','.join(sorted(missing)))
    return build_plan(
        config['lower'],
        config['upper'],
        config['lines'],
        alpha=config.get('alpha', 0.05),
        cluster_dimensions=config.get('cluster_dimensions', DIMENSIONS),
        require_multiple_lots=config.get('require_multiple_lots', True),
    )
