import math
import statistics
from collections import Counter


def describe(records):
    pairs = [row for row in records if row['valid'] and row['old_value'] is not None]
    return {'valid_pairs': len(pairs), 'difference_mean': statistics.fmean(row['new_value'] - row['old_value'] for row in pairs) if pairs else None}
