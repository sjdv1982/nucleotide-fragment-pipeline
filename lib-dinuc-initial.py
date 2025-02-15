print("TODO: # Remove dinucleotides with inter-nucleotide clashes")

import os
import sys
import itertools
import numpy as np
from seamless import Buffer
from nefertiti.functions.parse_mmcif import atomic_dtype
from crocodile.dinuc.from_ppdb import ppdb2coor


def err(msg):
    print(msg, file=sys.stderr)
    exit(1)


rna_struc_index, rna_strucs_data = Buffer.load("allpdb-rna-attract").deserialize(
    "mixed"
)

bases = ("A", "C", "G", "U")
dinuc_sequences = ["".join(s) for s in itertools.product(bases, repeat=2)]

lib_dinuc_coor = {}
lib_dinuc_origin = {}

template_pdbs = {}
for seq in dinuc_sequences:
    filename = f"templates/{seq}-template-ppdb.npy"
    template = np.load(filename)
    if template.dtype != atomic_dtype:
        err(f"Template '{filename}' does not contain a parsed PDB")
    template_pdbs[seq] = template

    lib_dinuc_coor[seq] = []
    lib_dinuc_origin[seq] = []

blacklist = [
    "2y8wB",
    "2y8yB",
    "3ercE",
    "3ercG",
    "3ercH",
    "4q5vB",
    "4v5gY",
    "4v5lY",
    "4v5pY",
    "4v5qY",
    "4v5rY",
    "4v5sY",
    "4v8qGB",
    "4w29X",
    "5afiY",
    "5bs3C",
    "7kqmB",
    "8agtWA",
]
from tqdm import tqdm

for code in tqdm(rna_struc_index):
    if code in blacklist:
        continue
    start, length = rna_struc_index[code]
    struc = rna_strucs_data[start : start + length]
    sequence, first_resids, coordinates = ppdb2coor(
        struc,
        template_pdbs=template_pdbs,
        rna=True,
        # do not ignore, since these structures come from aareduce and are therefore well-formatted
        ignore_unknown=False,
        ignore_missing=False,
        ignore_reordered=False,
    )
    for triseq, first_resid, coor in zip(sequence, first_resids, coordinates):
        lib_dinuc_coor[triseq].append(coor)
        lib_dinuc_origin[triseq].append((code, first_resid))

for seq in dinuc_sequences:
    coorf = f"lib-dinuc-initial-{seq}.npy"
    orif = f"lib-dinuc-initial-{seq}-origin.txt"
    if not lib_dinuc_coor[seq]:
        os.remove(coorf)
        os.remove(orif)
        continue
    coor = np.stack(lib_dinuc_coor[seq])
    np.save(coorf, coor)
    with open(orif, "w") as f:
        for code, first_resid in lib_dinuc_origin[seq]:
            print(code, first_resid, file=f)
