import argparse
import sys
import numpy as np

p = argparse.ArgumentParser("Applies a filter mask to an input file")
p.add_argument("input_file", help="Input file, either a .npy or a text file")
p.add_argument(
    "filter_file",
    help="""Filter mask file.
    In .npy format, as a 1D bool array with the same number of elements as the input file.
    1 indicates to keep, 0 indicates to discard.""",
)
p.add_argument(
    "output_file",
    help="Output file. Same format as the input file, with the mask applued.",
)


def err(*args):
    print(*args, file=sys.stderr)
    exit(1)


args = p.parse_args()
is_numpy = False
if args.input_file.endswith(".npy"):
    is_numpy = True
    inp = np.load(args.input_file)
    kind = "elements"
else:
    inp = []
    with open(args.input_file) as f:
        for l in f.readlines():
            inp.append(l)
    kind = "lines"
filter = np.load(args.filter_file)
if filter.ndim != 1:
    err("Filter must be a one-dimensional array")

if sorted(np.unique(filter).tolist()) not in ([0], [1], [0, 1]):
    err("Filter must contain ones and zeroes")

if len(filter) != len(inp):
    err(
        f"Filter does not match input: there are {len(inp)} {kind}, while the filter has {len(filter)} elements"
    )

if is_numpy:
    outp = inp[filter]
    np.save(args.output_file, outp)
else:
    outp = []
    for keep, line in zip(filter, inp):
        if keep:
            outp.append(line)
    with open(args.output_file, "w") as f:
        for l in outp:
            print(l, file=f, end="")
