import sys
import numpy as np

"""
Calculates the boolean intersection of two or more filters
"""


def err(*args):
    print(*args, file=sys.stderr)
    exit(1)


filter_files = sys.argv[1:-1]
outputfile = sys.argv[-1]

assert len(filter_files) >= 2

filters = []
for filter_file in filter_files:
    filter = np.load(filter_file)
    if not filter.ndim == 1:
        err(f"'{filter_file}' has wrong dimensions")
    if sorted(np.unique(filter).tolist()) not in ([0], [1], [0, 1]):
        err(f"'{filter_file}' must contain ones and zeroes")
    if filters and filter.shape != filters[0].shape:
        err("Filters must have the same shape")
    filters.append(filter)

out_filter = filters[0].astype(bool)
for filter in filters[1:]:
    out_filter &= filter.astype(bool)

np.save(outputfile, out_filter)
