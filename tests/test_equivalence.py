import copy
import csv
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from analysis_app.equivalence import analyze
from analysis_app.examples import run_examples
from analysis_app.io import load_records
from analysis_app.plans import (
    build_plan, fingerprint, from_config, normalize, verify,
)
from analysis_app.samplesize import power_curve
from analysis_app.schema import FIELDS
from analysis_app.stats import (
    holm_adjust, student_t_cdf, student_t_ppf, tost,
)


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data' / 'observations.csv'


def write_csv(directory, rows, name='input.csv'):
    path = Path(directory) / name
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    return path


def make_row(observation_id, sample, line='L0', lot=None, operator='O0',
             difference=0.1, valid=1, old_value=10.0, day=0):
    return {
        'observation_id': observation_id,
        'sample': sample,
        'line': line,
        'day': day,
        'new_value': 10.0 + difference,
        'old_value': old_value,
        'operator': operator,
        'lot': lot or sample,
        'valid': valid,
    }


class StatsTest(unittest.TestCase):
    def test_student_t_quantiles_match_reference(self):
        cases = [(0.975, 10, 2.228138852), (0.975, 1, 12.70620474),
                 (0.95, 2, 2.91998558), (0.995, 5, 4.032142984)]
        for probability, degrees, expected in cases:
            self.assertAlmostEqual(student_t_ppf(probability, degrees),
                                   expected, places=7)

    def test_student_t_cdf_is_symmetric(self):
        self.assertAlmostEqual(student_t_cdf(0.0, 7), 0.5)
        self.assertAlmostEqual(student_t_cdf(-1.5, 7),
                               1.0 - student_t_cdf(1.5, 7))

    def test_tost_interval_matches_decision(self):
        result = tost(0.1, 0.01, 10, -0.5, 0.5, 0.05)
        critical = student_t_ppf(0.95, 10)
        self.assertAlmostEqual(result['interval'][0], 0.1 - critical * 0.01)
        self.assertAlmostEqual(result['interval'][1], 0.1 + critical * 0.01)
        self.assertTrue(result['interval_within_bounds'])
        self.assertEqual(result['p_value'], max(result['p_lower'],
                                                result['p_upper']))

    def test_tost_at_bound_is_not_significant(self):
        result = tost(0.5, 0.01, 10, -0.5, 0.5, 0.05)
        self.assertAlmostEqual(result['p_upper'], 0.5)
        self.assertGreater(result['p_value'], 0.05)
        self.assertFalse(result['interval_within_bounds'])

    def test_tost_rejects_bad_bounds(self):
        with self.assertRaisesRegex(ValueError, '顺序'):
            tost(0.1, 0.01, 10, 0.5, -0.5, 0.05)
        with self.assertRaisesRegex(ValueError, '标准误'):
            tost(0.1, 0.0, 10, -0.5, 0.5, 0.05)

    def test_holm_is_monotone_and_controls_family(self):
        result = holm_adjust([0.01, 0.04, 0.03], 0.05)
        self.assertEqual(result['adjusted'], [0.03, 0.06, 0.06])
        self.assertEqual(result['reject'], [True, False, False])

    def test_holm_keeps_unestimable_hypotheses_in_family(self):
        result = holm_adjust([0.01, None], 0.05)
        self.assertEqual(result['adjusted'][1], 1.0)
        self.assertFalse(result['reject'][1])


class PlanTest(unittest.TestCase):
    def test_fingerprint_is_stable_and_order_independent(self):
        first = build_plan(-0.5, 0.5, ['L0', 'L1'])
        second = build_plan(-0.5, 0.5, ['L1', 'L0'])
        self.assertEqual(first['plan_fingerprint'],
                         second['plan_fingerprint'])
        self.assertEqual(verify(first), first)

    def test_bounds_are_in_raw_micrometres(self):
        plan = build_plan(-0.5, 0.5, ['L0'])
        self.assertEqual(plan['bounds'], {'lower': -0.5, 'upper': 0.5,
                                          'unit': '微米'})

    def test_bounds_order_and_unit_enforced(self):
        with self.assertRaisesRegex(ValueError, '下界必须小于上界'):
            build_plan(0.5, -0.5, ['L0'])
        config = {'lower': -0.5, 'upper': 0.5, 'lines': ['L0']}
        plan = from_config(config)
        changed = copy.deepcopy(plan)
        changed['bounds']['unit'] = '毫米'
        with self.assertRaisesRegex(ValueError, '单位'):
            normalize(changed)

    def test_tampered_plan_is_rejected(self):
        plan = build_plan(-0.5, 0.5, ['L0'])
        changed = copy.deepcopy(plan)
        changed['bounds']['upper'] = 0.9
        with self.assertRaisesRegex(ValueError, '指纹'):
            verify(changed)

    def test_plan_without_fingerprint_is_rejected(self):
        plan = build_plan(-0.5, 0.5, ['L0'])
        del plan['plan_fingerprint']
        with self.assertRaisesRegex(ValueError, '指纹'):
            verify(plan)

    def test_clustering_unit_must_be_sample(self):
        plan = build_plan(-0.5, 0.5, ['L0'])
        changed = copy.deepcopy(plan)
        changed['clustering']['unit'] = 'observation'
        with self.assertRaisesRegex(ValueError, 'sample'):
            normalize(changed)

    def test_unknown_plan_field_rejected(self):
        plan = build_plan(-0.5, 0.5, ['L0'])
        plan['extra'] = 1
        with self.assertRaisesRegex(ValueError, '未知字段'):
            normalize(plan)

    def test_config_requires_lines(self):
        with self.assertRaisesRegex(ValueError, '缺少字段'):
            from_config({'lower': -0.5, 'upper': 0.5})


class EquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = load_records(DATA)

    def plan(self, lower=-0.5, upper=0.5, lines=('L0', 'L1', 'L2')):
        return build_plan(lower, upper, list(lines))

    def test_fixture_declares_equivalence_with_all_lines(self):
        result = analyze(self.records, self.plan())
        self.assertEqual(result['verdict'], 'equivalent')
        self.assertEqual([c['group'] for c in result['comparisons']],
                         ['overall', 'line:L0', 'line:L1', 'line:L2'])
        self.assertEqual(result['verdict_counts']['equivalent'], 4)
        self.assertTrue(result['scope_ok'])
        for comparison in result['comparisons']:
            self.assertTrue(comparison['within_bounds'])
            self.assertTrue(comparison['multiple_testing']['reject'])

    def test_valid_zero_rows_are_excluded(self):
        changed = copy.deepcopy(self.records)
        for row in changed:
            if row['valid'] == 0:
                row['old_value'] = 999.0
        baseline = analyze(self.records, self.plan())
        altered = analyze(changed, self.plan())
        baseline.pop('input_fingerprint')
        altered.pop('input_fingerprint')
        self.assertEqual(altered, baseline)

    def test_valid_values_outside_binary_rejected(self):
        changed = copy.deepcopy(self.records)
        changed[0]['valid'] = 2
        with self.assertRaisesRegex(ValueError, 'valid'):
            analyze(changed, self.plan())

    def test_repeated_samples_collapse_to_units(self):
        rows = []
        for index, (sample, lot, difference) in enumerate([
                ('S0', 'A', 0.09), ('S0', 'A', 0.11),
                ('S1', 'A', 0.08), ('S1', 'A', 0.12),
                ('S2', 'B', 0.10), ('S2', 'B', 0.10),
                ('S3', 'B', 0.09), ('S3', 'B', 0.11)]):
            rows.append(make_row(f'o{index}', sample, lot=lot,
                                 operator=f'O{index % 2}',
                                 difference=difference, day=index))
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0']))
        overall = result['comparisons'][0]
        self.assertEqual(result['valid_pairs'], 8)
        self.assertEqual(result['samples'], 4)
        self.assertEqual(overall['units'], 4)
        self.assertEqual(overall['observations'], 8)
        self.assertEqual(overall['repeated_samples'], 4)
        self.assertEqual(overall['maximum_repeats'], 2)

    def test_no_valid_pairs_is_unavailable(self):
        rows = [make_row('o0', 'S0', valid=0, old_value=None),
                make_row('o1', 'S1', valid=0, old_value=None)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0']))
        self.assertEqual(result['verdict'], 'unavailable')
        self.assertEqual(result['valid_pairs'], 0)
        self.assertEqual(result['comparisons'][0]['status'], 'unavailable')

    def test_zero_variance_is_unavailable_not_equivalent(self):
        rows = [make_row(f'o{index}', f'S{index}', lot='ABCD'[index % 4],
                         operator=f'O{index % 3}') for index in range(8)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0']))
        self.assertEqual(result['verdict'], 'unavailable')
        self.assertEqual(result['comparisons'][0]['reason'],
                         '配对差值方差为零')
        self.assertIsNone(result['comparisons'][0]['p_value'])

    def test_mean_exactly_at_bound_is_not_proven(self):
        differences = [0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.45, 0.55]
        rows = [make_row(f'o{index}', f'S{index}', lot='ABCD'[index % 4],
                         operator=f'O{index % 3}',
                         difference=differences[index])
                for index in range(8)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0']))
        overall = result['comparisons'][0]
        self.assertEqual(result['verdict'], 'not_proven')
        self.assertAlmostEqual(overall['estimate'], 0.5)
        self.assertFalse(overall['within_bounds'])
        self.assertGreater(overall['p_value'], 0.05)

    def test_asymmetric_bounds_are_honoured(self):
        differences = [0.08, 0.12, 0.09, 0.11, 0.10, 0.095, 0.105, 0.1]
        rows = [make_row(f'o{index}', f'S{index}', lot='ABCD'[index % 4],
                         operator=f'O{index % 3}',
                         difference=differences[index])
                for index in range(8)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        wide = analyze(records, self.plan(-0.05, 0.3, ['L0']))
        narrow = analyze(records, self.plan(-0.05, 0.099, ['L0']))
        self.assertEqual(wide['verdict'], 'equivalent')
        self.assertEqual(narrow['verdict'], 'not_proven')
        self.assertNotEqual(wide['plan_fingerprint'],
                            narrow['plan_fingerprint'])

    def test_single_lot_cannot_claim_cross_lot(self):
        differences = [0.08, 0.12, 0.09, 0.11, 0.10, 0.095, 0.105, 0.1]
        rows = [make_row(f'o{index}', f'S{index}', lot='A',
                         operator=f'O{index % 3}',
                         difference=differences[index])
                for index in range(8)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0']))
        overall = result['comparisons'][0]
        self.assertEqual(result['verdict'], 'not_proven')
        self.assertFalse(result['scope_ok'])
        self.assertTrue(result['scope_issues']['single_lot'])
        self.assertFalse(result['cross_lot'])
        self.assertIn('lot', overall['omitted_clusters'])
        self.assertEqual(overall['reason'], '只有一批材料，不能跨批外推')

    def test_one_failing_line_blocks_equivalence(self):
        rows = []
        good = [0.08, 0.12, 0.09, 0.11, 0.10, 0.095, 0.105, 0.1]
        bad = [1.9, 2.1, 2.0, 1.95, 2.05, 2.0, 1.9, 2.1]
        for index in range(8):
            rows.append(make_row(f'g{index}', f'G{index}', line='L0',
                                 lot='ABCD'[index % 4],
                                 operator=f'O{index % 3}',
                                 difference=good[index]))
        for index in range(8):
            rows.append(make_row(f'b{index}', f'B{index}', line='L1',
                                 lot='ABCD'[index % 4],
                                 operator=f'O{index % 3}',
                                 difference=bad[index]))
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0', 'L1']))
        verdicts = {c['group']: c['verdict'] for c in result['comparisons']}
        self.assertEqual(verdicts['line:L0'], 'equivalent')
        self.assertEqual(verdicts['line:L1'], 'not_proven')
        self.assertEqual(verdicts['overall'], 'not_proven')
        self.assertEqual(result['verdict'], 'not_proven')
        self.assertEqual(result['multiple_testing']['family'],
                         ['overall', 'line:L0', 'line:L1'])

    def test_holm_adjusts_line_p_values(self):
        result = analyze(self.records, self.plan())
        for comparison in result['comparisons']:
            self.assertGreaterEqual(comparison['p_adjusted'],
                                    comparison['p_value'])
            self.assertEqual(comparison['multiple_testing']['family_size'], 4)

    def test_unknown_line_is_a_scope_problem(self):
        changed = copy.deepcopy(self.records)
        target = next(row for row in changed if row['valid'] == 1)
        target['line'] = 'L9'
        result = analyze(changed, self.plan())
        self.assertIn('L9', result['scope_issues']['unknown_lines'])
        self.assertFalse(result['scope_ok'])
        self.assertEqual(result['verdict'], 'not_proven')

    def test_invalid_unknown_line_is_not_evidence(self):
        changed = copy.deepcopy(self.records)
        target = next(row for row in changed if row['valid'] == 0)
        target['line'] = 'L9'
        result = analyze(changed, self.plan())
        self.assertEqual(result['scope_issues']['unknown_lines'], [])
        self.assertTrue(result['scope_ok'])
        self.assertEqual(result['verdict'], 'equivalent')

    def test_missing_planned_line_blocks_equivalence(self):
        result = analyze(self.records, self.plan(lines=('L0', 'L1', 'L9')))
        self.assertIn('L9', result['scope_issues']['missing_planned_lines'])
        self.assertFalse(result['scope_ok'])
        self.assertEqual(result['verdict'], 'not_proven')
        verdicts = {c['group']: c['verdict'] for c in result['comparisons']}
        self.assertEqual(verdicts['line:L9'], 'unavailable')

    def test_sample_spanning_lines_is_flagged(self):
        rows = []
        differences = [0.08, 0.12, 0.09, 0.11, 0.10, 0.095, 0.105, 0.1]
        for index, difference in enumerate(differences):
            rows.append(make_row(f'a{index}', f'A{index}', line='L0',
                                 lot='ABCD'[index % 4],
                                 operator=f'O{index % 3}',
                                 difference=difference))
            rows.append(make_row(f'b{index}', f'B{index}', line='L1',
                                 lot='ABCD'[index % 4],
                                 operator=f'O{index % 3}',
                                 difference=difference))
        rows.append(make_row('x0', 'X0', line='L0', lot='A', difference=0.1))
        rows.append(make_row('x1', 'X0', line='L1', lot='A', difference=0.1))
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        result = analyze(records, self.plan(-0.5, 0.5, ['L0', 'L1']))
        self.assertEqual(result['scope_issues']['samples_spanning_lines'], ['X0'])
        self.assertFalse(result['scope_ok'])
        self.assertEqual(result['verdict'], 'not_proven')
        verdicts = {c['group']: c['verdict'] for c in result['comparisons']}
        self.assertEqual(verdicts['line:L0'], 'equivalent')
        self.assertEqual(verdicts['line:L1'], 'equivalent')
        self.assertEqual(verdicts['overall'], 'not_proven')

    def test_plan_must_be_frozen_before_analysis(self):
        plan = self.plan()
        del plan['plan_fingerprint']
        with self.assertRaisesRegex(ValueError, '指纹'):
            analyze(self.records, plan)

    def test_input_order_does_not_change_analysis(self):
        shuffled = copy.deepcopy(self.records)
        random.Random(7).shuffle(shuffled)
        self.assertEqual(analyze(shuffled, self.plan()),
                         analyze(self.records, self.plan()))

    def test_result_carries_both_fingerprints(self):
        result = analyze(self.records, self.plan())
        self.assertEqual(result['plan_fingerprint'],
                         fingerprint(normalize(self.plan())))
        self.assertEqual(len(result['input_fingerprint']), 64)


class SampleSizeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = load_records(DATA)

    def test_power_increases_with_sample_size(self):
        result = power_curve(self.records, -0.2, 0.2, 0.17, 0.05,
                             [12, 24, 36, 60, 90, 120])
        powers = [point['power'] for point in result['curve']]
        self.assertEqual(powers, sorted(powers))
        self.assertEqual(result['smallest_size_at_target'], 36)
        self.assertIn('不是观测证据', result['note'])

    def test_larger_true_effect_reduces_power(self):
        close = power_curve(self.records, -0.2, 0.2, 0.05, 0.05, [24])
        far = power_curve(self.records, -0.2, 0.2, 0.19, 0.05, [24])
        self.assertGreater(close['curve'][0]['power'],
                           far['curve'][0]['power'])

    def test_design_effect_reduces_power(self):
        plain = power_curve(self.records, -0.2, 0.2, 0.1, 0.05, [24])
        clustered = power_curve(self.records, -0.2, 0.2, 0.1, 0.05, [24],
                                design_effect=9.0)
        self.assertGreater(plain['curve'][0]['power'],
                           clustered['curve'][0]['power'])

    def test_invalid_sample_sizes_rejected(self):
        with self.assertRaisesRegex(ValueError, '样本量'):
            power_curve(self.records, -0.2, 0.2, 0.1, 0.05, [1])

    def test_zero_pilot_variance_rejected(self):
        rows = [make_row(f'o{index}', f'S{index}', lot='AB'[index % 2],
                         operator=f'O{index % 2}')
                for index in range(4)]
        with tempfile.TemporaryDirectory() as directory:
            records = load_records(write_csv(directory, rows))
        with self.assertRaisesRegex(ValueError, '方差为零'):
            power_curve(records, -0.5, 0.5, 0.1, 0.05, [12])


class CliTest(unittest.TestCase):
    def run_cli(self, *arguments):
        return subprocess.run([sys.executable, '-m', 'analysis_app',
                               *map(str, arguments)],
                              cwd=ROOT, text=True, capture_output=True)

    def test_plan_analyze_and_check_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'config.json'
            config.write_text(json.dumps({
                'lower': -0.5, 'upper': 0.5, 'lines': ['L0', 'L1', 'L2'],
            }), encoding='utf-8')
            plan_path = Path(directory) / 'plan.json'
            made = self.run_cli('plan', '--config', config, '--output', plan_path)
            self.assertEqual(made.returncode, 0, made.stderr)
            plan = json.loads(plan_path.read_text(encoding='utf-8'))
            checked = self.run_cli('check-plan', '--plan', plan_path)
            self.assertEqual(checked.returncode, 0, checked.stderr)
            self.assertEqual(json.loads(checked.stdout)['plan_fingerprint'],
                             plan['plan_fingerprint'])
            output = Path(directory) / 'result.json'
            analyzed = self.run_cli('analyze', '--plan', plan_path,
                                    '--output', output)
            self.assertEqual(analyzed.returncode, 0, analyzed.stderr)
            result = json.loads(output.read_text(encoding='utf-8'))
            self.assertEqual(result['plan_fingerprint'],
                             plan['plan_fingerprint'])
            self.assertEqual(result['verdict'], 'equivalent')

    def test_existing_plan_is_not_overwritten_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'config.json'
            config.write_text(json.dumps({
                'lower': -0.5, 'upper': 0.5, 'lines': ['L0'],
            }), encoding='utf-8')
            plan_path = Path(directory) / 'plan.json'
            first = self.run_cli('plan', '--config', config, '--output', plan_path)
            self.assertEqual(first.returncode, 0, first.stderr)
            before = plan_path.read_text(encoding='utf-8')
            second = self.run_cli('plan', '--config', config, '--output', plan_path)
            self.assertEqual(second.returncode, 2)
            self.assertEqual(plan_path.read_text(encoding='utf-8'), before)

    def test_tampered_plan_cannot_be_analyzed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'plan.json'
            plan = build_plan(-0.5, 0.5, ['L0'])
            plan['bounds']['upper'] = 0.9
            path.write_text(json.dumps(plan), encoding='utf-8')
            run = self.run_cli('analyze', '--plan', path,
                               '--output', Path(directory) / 'result.json')
            self.assertEqual(run.returncode, 2)
            self.assertIn('指纹', run.stderr)
            self.assertFalse((Path(directory) / 'result.json').exists())

    def test_new_plan_does_not_overwrite_old_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.json'
            old_plan_path = Path(directory) / 'old.json'
            new_plan_path = Path(directory) / 'new.json'
            old_plan_path.write_text(json.dumps(build_plan(-0.5, 0.5, ['L0'])),
                                     encoding='utf-8')
            new_plan_path.write_text(json.dumps(build_plan(-0.9, 0.9, ['L0'])),
                                     encoding='utf-8')
            first = self.run_cli('analyze', '--plan', old_plan_path,
                                 '--output', output)
            self.assertEqual(first.returncode, 0, first.stderr)
            before = output.read_text(encoding='utf-8')
            second = self.run_cli('analyze', '--plan', new_plan_path,
                                  '--output', output)
            self.assertEqual(second.returncode, 2)
            self.assertIn('指纹', second.stderr)
            self.assertEqual(output.read_text(encoding='utf-8'), before)

    def test_same_plan_may_refresh_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'result.json'
            plan_path = Path(directory) / 'plan.json'
            plan_path.write_text(json.dumps(build_plan(-0.5, 0.5, ['L0'])),
                                 encoding='utf-8')
            for _ in range(2):
                run = self.run_cli('analyze', '--plan', plan_path,
                                   '--output', output)
                self.assertEqual(run.returncode, 0, run.stderr)

    def test_power_requires_explicit_effect(self):
        with tempfile.TemporaryDirectory() as directory:
            plan_path = Path(directory) / 'plan.json'
            plan_path.write_text(json.dumps(build_plan(-0.2, 0.2, ['L0'])),
                                 encoding='utf-8')
            missing = self.run_cli('power', '--plan', plan_path)
            self.assertEqual(missing.returncode, 2)
            given = self.run_cli('power', '--plan', plan_path,
                                 '--effect', '0.17', '--sizes', '12,24')
            self.assertEqual(given.returncode, 0, given.stderr)
            result = json.loads(given.stdout)
            self.assertEqual(result['assumed_true_difference'], 0.17)
            self.assertEqual(len(result['curve']), 2)

    def test_examples_reproduce_expected_boundaries(self):
        run = self.run_cli('examples')
        self.assertEqual(run.returncode, 0, run.stderr)
        result = json.loads(run.stdout)
        self.assertEqual(result['failed'], 0)
        self.assertEqual(result['passed'], 7)
        names = {case['case'] for case in result['cases']}
        self.assertTrue({'no-valid-pairs', 'zero-variance', 'mean-at-bound',
                         'asymmetric-bounds', 'repeated-samples',
                         'single-lot', 'line-multiplicity'} <= names)


if __name__ == '__main__':
    unittest.main()
