"""可复算的边界例子。

每个例子包含合成观测、内联计划配置和预期判定。运行 examples 命令会
重新固化计划、重新分析，并把结论、估计值、区间、聚类信息与 Holm 校正
值逐项与清单登记的期望值比对；任一项不符即失败。
"""

import json
import math
from pathlib import Path

from .equivalence import analyze
from .io import load_records
from .plans import freeze


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = ROOT / 'examples' / 'cases.json'


def _close(actual, expected, tolerance):
    if expected is None:
        return actual is None
    if actual is None:
        return False
    return math.isclose(actual, expected, rel_tol=tolerance,
                        abs_tol=tolerance)


def _compare(actual, expected, tolerance, prefix=''):
    """递归比较结果子集；清单只声明关心的字段。"""
    mismatches = []
    if isinstance(expected, dict):
        for key, value in expected.items():
            if not isinstance(actual, dict) or key not in actual:
                mismatches.append(f'{prefix}{key} 缺失')
                continue
            mismatches.extend(_compare(actual[key], value, tolerance,
                                       f'{prefix}{key}.'))
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            mismatches.append(f'{prefix} 列表长度不符')
        else:
            for index, value in enumerate(expected):
                mismatches.extend(_compare(actual[index], value, tolerance,
                                           f'{prefix}{index}.'))
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        if not _close(actual, expected, tolerance):
            mismatches.append(f'{prefix} 期望 {expected} 实际 {actual}')
    elif actual != expected:
        mismatches.append(f'{prefix} 期望 {expected!r} 实际 {actual!r}')
    return mismatches


def _resolve(result, path):
    """按清单路径取值，支持 comparisons[group]。"""
    actual = result
    for part in path.split('.'):
        if part.startswith('comparisons['):
            group = part[len('comparisons['):-1]
            found = None
            for comparison in result['comparisons']:
                if comparison['group'] == group:
                    found = comparison
                    break
            actual = found
        elif isinstance(actual, dict):
            actual = actual.get(part)
        else:
            return None
    return actual


def run_examples(manifest_path=None):
    """重新固化并分析清单中的全部例子。"""
    manifest_path = Path(manifest_path or DEFAULT_MANIFEST)
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    base = manifest_path.parent
    results = []
    for case in manifest['cases']:
        source = (base / case['input']).resolve()
        records = load_records(source)
        plan = freeze(case['plan'])
        result = analyze(records, plan)
        observed = {path: _resolve(result, path) for path in case['expect']}
        tolerance = case.get('tolerance', manifest.get('tolerance', 1e-9))
        mismatches = []
        for path, expected in case['expect'].items():
            mismatches.extend(_compare(observed[path], expected, tolerance,
                                       f'{path} = '))
        results.append({
            'case': case['name'],
            'description': case['description'],
            'input': case['input'],
            'plan_fingerprint': plan['plan_fingerprint'],
            'observed': observed,
            'mismatches': mismatches,
            'passed': not mismatches,
        })
    return results
