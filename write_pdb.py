"""Write binary format back to PDB file

Author: Sjoerd de Vries
License: GPLv3, https://www.gnu.org/licenses/gpl-3.0.en.html

This file is part of the Nefertiti project.
"""

import numpy as np
from typing import List
import importlib


def _import(name):
    name2 = name.lstrip(".")
    if name2 not in globals():
        if __name__.find(".") == -1:
            name = name2
        mod = importlib.import_module(name, __name__)
        globals()[name2] = mod


code = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "UNK": "X",
}
code_rev = {v: k for k, v in code.items()}

# adapted from Bio.PDB.PDBIO.py

_ATOM_FORMAT_STRING = (
    "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"
)


def write_pdb_atom(atom) -> str:
    name = atom["name"].decode()
    if not name.startswith("H"):
        name = " " + name
    occ = atom["occupancy"]
    if occ > 100:
        occ = 100
    if occ < 0:
        occ = 0
    args = (
        "ATOM  " if atom["hetero"].decode().strip() == "" else "HETATM",
        atom["index"] % 100000,
        name,
        atom["altloc"].decode(),
        atom["resname"].decode(),
        atom["chain"].decode()[0],
        atom["resid"],
        atom["icode"].decode(),
        atom["x"],
        atom["y"],
        atom["z"],
        occ,
        np.minimum(atom["bfactor"], 100),
        atom["segid"].decode(),
        atom["element"].decode(),
        "",
    )
    return _ATOM_FORMAT_STRING % args


# /adapted


def write_pdb(struc: np.ndarray, ter=False) -> str:
    _import("..parse_pdb")
    assert struc.dtype == parse_pdb.atomic_dtype
    assert struc.ndim == 1
    pdb = ""
    curr_id = None
    for atom in struc:
        line = write_pdb_atom(atom)
        new_id = atom["chain"], atom["resid"]
        if curr_id is not None and new_id != curr_id:
            if new_id[0] != curr_id[0] or new_id[1] != curr_id[1] + 1:
                pdb += "TER\n"
        pdb += line
        curr_id = new_id
    return pdb


def write_multi_pdb(struc: np.ndarray) -> str:
    _import("..parse_pdb")
    assert struc.dtype == parse_pdb.atomic_dtype
    assert struc.ndim == 2
    pdb = ""
    for n, substruc in enumerate(struc):
        pdb += "MODEL {}\n".format(n + 1)
        pdb += write_pdb(substruc)
        pdb += "ENDMDL\n"
    return pdb


def build_pdb_backbone(
    struc: np.ndarray,
    bb_atoms: List[str],
    sequence: str = None,
    *,
    atomindex_offset=0,
    resid_offset=0,
) -> np.ndarray:
    _import("..parse_pdb")
    assert struc.ndim == 3, struc.shape
    assert struc.shape[1] == len(bb_atoms), struc.shape
    assert struc.shape[2] in (3, 4), struc.shape
    assert len(struc) < 10000
    if sequence is not None:
        assert len(sequence) == len(struc)
    newstruc = np.zeros((len(struc), len(bb_atoms)), parse_pdb.atomic_dtype)
    for resnr in range(len(struc)):
        res = newstruc[resnr]
        for anr in range(len(bb_atoms)):
            s = struc[resnr][anr]
            a = newstruc[resnr][anr]
            a["model"] = 1
            a["hetero"] = ""
            a["name"] = bb_atoms[anr].encode()
            a["altloc"] = " "
            resname = code_rev[sequence[resnr]] if sequence is not None else "ALA"
            a["resname"] = resname.encode()
            a["chain"] = " "
            a["index"] = len(bb_atoms) * resnr + anr + atomindex_offset + 1
            a["icode"] = " "
            a["resid"] = resnr + resid_offset + 1
            a["x"] = s[0]
            a["y"] = s[1]
            a["z"] = s[2]
            a["occupancy"] = 1
            a["segid"] = ""
            a["element"] = bb_atoms[anr][0].encode()
    return newstruc.flatten()


def build_pdb_fragment_backbone(
    struc: np.ndarray, bb_atoms: List[str], sequence: str = None
):
    _import("..parse_pdb")
    assert struc.ndim == 4
    assert struc.shape[2] == len(bb_atoms)
    assert struc.shape[3] in (3, 4)
    struc = struc[:, :, :, :3]
    nfrag, fraglen = struc.shape[:2]
    if sequence is not None:
        assert len(sequence) == nfrag + fraglen - 1, (len(sequence), nfrag, fraglen)
    assert fraglen < 100
    offset = 10 if fraglen < 10 else 100
    assert len(sequence) * offset < 10000

    atoms = np.empty(nfrag * fraglen * len(bb_atoms), parse_pdb.atomic_dtype)
    for n in range(nfrag):
        stride = fraglen * len(bb_atoms)
        atoms[n * stride : (n + 1) * stride] = build_pdb_backbone(
            struc[n],
            bb_atoms,
            sequence[n : n + fraglen],
            atomindex_offset=n * len(bb_atoms) * fraglen,
            resid_offset=n * offset,
        )
    return atoms


if __name__ == "__main__":
    import sys
    import parse_pdb

    npy_file = sys.argv[1]
    outfile = sys.argv[2]
    struc = np.load(npy_file)
    pdb_data = write_pdb(struc, ter=True)
    with open(outfile, "w") as f:
        f.write(pdb_data)
