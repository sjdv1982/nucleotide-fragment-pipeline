#!/usr/bin/env python3
# Author: Sjoerd de Vries (CNRS)

import numpy as np
import sys, argparse


def err(*args):
    print(*args, file=sys.stderr)
    exit(1)


p = argparse.ArgumentParser(
    description="""Filter di/trinucleotides, removing fragments with disconnected bonds.
    Bonds are first detected in the template structure and then checked in the fragments."""
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
               Indicates if the fragment is kept (1) or rejected because of bad bonds (0).""",
)
p.add_argument(
    "--detection-threshold",
    help="Bond detection threshold in A",
    default=1.7,
    type=float,
)
p.add_argument(
    "--filter-threshold",
    help="Bond distance threshold for filtering the fragments",
    default=3.0,
    type=float,
)


args = p.parse_args()

detection_threshold = args.detection_threshold
if detection_threshold <= 0:
    err("Detection threshold must be positive")
filter_threshold = args.filter_threshold
if filter_threshold <= 0:
    err("Filter threshold must be positive")

fragments = np.load(getattr(args, "fragments.npy"))
tmpl = np.load(getattr(args, "template.npy"))
outp = getattr(args, "filter.npy")

if len(tmpl) != fragments.shape[1]:
    err(
        f"Template and fragments don't have the same number of atoms: {len(tmpl)} vs {fragments.shape[1]}"
    )

tmpl_coor = np.stack((tmpl["x"], tmpl["y"], tmpl["z"]), axis=1)
d = tmpl_coor[:, None, :] - tmpl_coor[None, :, :]
tmpl_dist = np.sqrt((d * d).sum(axis=2))
tmpl_bondmat = tmpl_dist < detection_threshold
p1, p2 = np.where(tmpl_bondmat)
bonds = []
for pp1, pp2 in zip(p1, p2):
    if pp1 >= pp2:
        continue

    a1 = tmpl[pp1]
    a2 = tmpl[pp2]

    """
    print(a1)
    print(a2)
    print(tmpl_dist[pp1, pp2])
    print()
    """

    if a1["resid"] != a2["resid"]:
        if a2["name"] != b"P":
            # print("REJECT")
            continue
    bonds.append((pp1, pp2))

p1 = [pp1 for pp1, _ in bonds]
p2 = [pp2 for _, pp2 in bonds]
frag_p1 = fragments[:, p1]
frag_p2 = fragments[:, p2]
d = frag_p1 - frag_p2
dis = np.sqrt((d * d).sum(axis=2))
maxdis = dis.max(axis=1)

keep = maxdis < filter_threshold
print(f"Keep {keep.sum()}/{len(fragments)} fragments", file=sys.stderr)
np.save(outp, keep)
