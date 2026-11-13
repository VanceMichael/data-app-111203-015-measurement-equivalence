import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

from .equivalence import analyze
from .examples import run_examples
from .io import load_records
from .plans import from_config, verify
from .samplesize import power_curve
from .summary import summarize


def save_json(path, result):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.analysis-', dir=destination.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, sort_keys=True, allow_nan=False, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_json(path):
    with Path(path).open(encoding='utf-8') as handle:
        return json.load(handle)


def _print(result):
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False))


def _plan_command(args):
    destination = Path(args.output)
    if destination.exists():
        raise ValueError('计划文件已存在；登记后的计划不可覆盖，'
                         '如需新计划请使用新路径')
    plan = from_config(read_json(args.config))
    save_json(destination, plan)
    return {'plan': str(destination), 'plan_fingerprint': plan['plan_fingerprint']}


def _analyze_command(args):
    plan = verify(read_json(args.plan))
    destination = Path(args.output) if args.output else None
    if destination is not None and destination.exists():
        previous = read_json(destination)
        if previous.get('plan_fingerprint') != plan['plan_fingerprint']:
            raise ValueError('既有分析基于不同计划指纹，拒绝覆盖；请写到新路径')
    return analyze(load_records(args.input), plan)


def _power_command(args):
    plan = verify(read_json(args.plan))
    sizes = [int(item) for item in args.sizes.split(',') if item]
    if args.effect is None:
        raise ValueError('必须给出用户假定的真实偏移 --effect；'
                         '不得用观测到的差值代替')
    result = power_curve(load_records(args.input), plan['bounds']['lower'],
                         plan['bounds']['upper'], args.effect, plan['alpha'],
                         sizes, design_effect=args.design_effect)
    result['plan_fingerprint'] = plan['plan_fingerprint']
    return result


def _check_plan_command(args):
    plan = verify(read_json(args.plan))
    return {'plan_fingerprint': plan['plan_fingerprint'], 'verified': True}


def _examples_command(args):
    results = run_examples(args.manifest)
    failed = [item for item in results if not item['passed']]
    summary = {'cases': results, 'passed': len(results) - len(failed),
               'failed': len(failed)}
    if failed:
        raise ValueError('例子与预期不符：' + ','.join(item['case'] for item in failed))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='本地观测数据分析')
    commands = parser.add_subparsers(dest='command', required=True)

    summary = commands.add_parser('summary', help='按业务维度核对观测口径')
    summary.add_argument('--input', type=Path, default=Path('data/observations.csv'))
    summary.add_argument('--output', type=Path)

    plan = commands.add_parser('plan', help='登记并固化等效性研究计划')
    plan.add_argument('--config', type=Path, required=True,
                      help='原始微米单位的界值与亚组配置')
    plan.add_argument('--output', type=Path, required=True)

    analyze_command = commands.add_parser('analyze', help='按固化计划做配对 TOST')
    analyze_command.add_argument('--plan', type=Path, required=True)
    analyze_command.add_argument('--input', type=Path,
                                 default=Path('data/observations.csv'))
    analyze_command.add_argument('--output', type=Path)

    power = commands.add_parser('power', help='基于先导方差的功效曲线')
    power.add_argument('--plan', type=Path, required=True)
    power.add_argument('--input', type=Path,
                       default=Path('data/observations.csv'))
    power.add_argument('--effect', type=float, required=True,
                       help='用户假定的真实偏移（微米），不是观测结果')
    power.add_argument('--sizes', default='12,18,24,30,36,48,60,90,120')
    power.add_argument('--design-effect', type=float, default=1.0,
                       help='可选的方差放大系数，用于计入聚类设计效应')
    power.add_argument('--output', type=Path)

    check = commands.add_parser('check-plan', help='核对计划指纹')
    check.add_argument('--plan', type=Path, required=True)

    examples = commands.add_parser('examples', help='复算边界例子并核对预期')
    examples.add_argument('--manifest', type=Path,
                         default=Path('examples/cases.json'))

    args = parser.parse_args(argv)
    try:
        if args.command == 'summary':
            result = summarize(load_records(args.input))
        elif args.command == 'plan':
            result = _plan_command(args)
        elif args.command == 'analyze':
            result = _analyze_command(args)
        elif args.command == 'power':
            result = _power_command(args)
        elif args.command == 'check-plan':
            result = _check_plan_command(args)
        else:
            result = _examples_command(args)
        if args.command in ('summary', 'analyze', 'power') and args.output:
            save_json(args.output, result)
        _print(result)
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
