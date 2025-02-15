#!/usr/bin/env python3
# Authors: Isaure Chauvot de Beauchene (CNRS), Sjoerd de Vries (CNRS)

import numpy as np
import sys, argparse


def err(*args):
    print(*args, file=sys.stderr)
    exit(1)


p = argparse.ArgumentParser(
    description="Remove di/trinucleotides with inter-nucleotide clashes"
)
p.add_argument(
    "fragments.npy", help="Input fragment coordinates, shape=(X,Y,3) floats."
)
p.add_argument(
    "template.npy",
    help="Di/trinucleotide template structure in parsed PDB (ppdb) numpy format, with Y atoms",
)
p.add_argument(
    "filter.npy",
    help="""Output mask of X bools.
               Indicates if the fragment is kept (1) or rejected because of clashes (0).""",
)
p.add_argument("--threshold", help="Clashing threshold in A", default=2.0, type=float)


def mindist(at1, at2):
    # print(at1.shape)
    # print(at2.shape)
    dif = at1[:, :, None] - at2[:, None, :]
    dif2 = np.sum(dif * dif, axis=3) ** 0.5
    mdif2 = dif2.min(axis=(1, 2))
    return mdif2


args = p.parse_args()

threshold = args.threshold
if threshold <= 0:
    err("Threshold must be positive")

fragments = np.load(getattr(args, "fragments.npy"))
tmpl = np.load(getattr(args, "template.npy"))
outp = getattr(args, "filter.npy")

if len(tmpl) != fragments.shape[1]:
    err(
        f"Template and fragments don't have the same number of atoms: {len(tmpl)} vs {fragments.shape[1]}"
    )
resids = np.unique(tmpl["resid"])

fragments_res = []
tmpl_res = []
for resid in resids:
    inds = tmpl["resid"] == resid
    tmpl_res.append(tmpl[inds])
    fragments_res.append(fragments[:, inds])

nres = len(resids)
keep = np.ones(len(fragments), bool)
for n1 in range(nres):
    for n2 in range(n1 + 1, nres):
        print(f"Interaction {n1} - {n2}", file=sys.stderr)
        fragments_res1 = fragments_res[n1]
        fragments_res2 = fragments_res[n2]
        if n2 - n1 == 1:
            tmpl_res1 = tmpl_res[n1]
            mask1 = tmpl_res1["name"] != b"O3'"
            fragments_res1 = fragments_res1[:, mask1]
            tmpl_res2 = tmpl_res[n2]
            mask2 = tmpl_res2["name"] != b"P"
            fragments_res2 = fragments_res2[:, mask2]
        d = mindist(fragments_res1, fragments_res2)
        curr_keep = d > threshold
        keep &= curr_keep
print(f"Keep {keep.sum()}/{len(fragments)} fragments", file=sys.stderr)
np.save(outp, keep)
