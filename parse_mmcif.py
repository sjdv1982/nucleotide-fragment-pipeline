"""Parse a mmCIF file into binary format

Author: Sjoerd de Vries
License: GPLv3, https://www.gnu.org/licenses/gpl-3.0.en.html

This file is part of the Nefertiti project.
"""

import warnings
import Bio.PDB
from Bio.PDB.StructureBuilder import PDBConstructionWarning

warnings.simplefilter("ignore", PDBConstructionWarning)
from io import StringIO
import numpy as np
from typing import List

atomic_dtype = [
    ("model", "uint16"),
    ("hetero", "S1"),
    ("name", "S4"),
    ("altloc", "S1"),
    ("resname", "S3"),
    ("chain", "S4"),
    ("index", "uint32"),
    ("icode", "S1"),
    ("resid", "uint16"),
    ("x", "float32"),
    ("y", "float32"),
    ("z", "float32"),
    ("occupancy", "float32"),
    ("bfactor", "float32"),
    ("segid", "S4"),
    ("element", "S2"),
]

atomic_dtype = np.dtype(atomic_dtype, align=True)


def parse_mmcif(mmcif_data: str, auth_chains=True, auth_residues=True) -> np.ndarray:

    mmcif_obj = StringIO(mmcif_data)

    try:
        p = Bio.PDB.MMCIFParser(auth_chains=auth_chains, auth_residues=auth_residues)
    except TypeError:
        if auth_residues is not True or auth_chains is not True:
            raise TypeError(
                "auth_residues / auth_chains is not supported in your version of Biopython"
            )
        p = Bio.PDB.MMCIFParser()
    struc = p.get_structure("PDB", mmcif_obj)
    natoms = len(list(struc.get_atoms()))
    atomstate = np.zeros(natoms, dtype=atomic_dtype)

    a = atomstate
    count = 0
    for modelnr, model in enumerate(struc.get_models()):
        atomlist = list(model.get_atoms())
        # atomlist.sort(key=lambda atom: atom.serial_number) # does not work
        index = 0  ###
        for atom in atomlist:
            residue = atom.get_parent()
            hetero, resid, icode = residue.get_id()
            segid = residue.segid
            resname = residue.resname
            chainid = residue.get_parent().id
            aa = a[count]
            aa["model"] = modelnr + 1
            aa["hetero"] = hetero
            aa["name"] = atom.name
            aa["altloc"] = atom.altloc
            aa["resname"] = resname
            aa["chain"] = chainid
            # aa["index"] = atom.serial_number
            aa["index"] = index + 1
            index += 1
            aa["icode"] = icode
            aa["resid"] = resid
            aa["x"] = atom.coord[0]
            aa["y"] = atom.coord[1]
            aa["z"] = atom.coord[2]
            occ = atom.occupancy
            if occ is None or occ < 0:
                occ = 0
            aa["occupancy"] = occ
            aa["segid"] = segid
            aa["element"] = atom.element
            count += 1
    return atomstate


def get_xyz(struc: np.ndarray) -> np.ndarray:
    struc2 = struc.flatten()
    result = np.zeros((len(struc2), 3))
    result[:, 0] = struc2["x"]
    result[:, 1] = struc2["y"]
    result[:, 2] = struc2["z"]
    return result.reshape(struc.shape + (3,))


def get_backbone(struc: np.ndarray, backbone_atoms: List[str]) -> np.ndarray:
    assert struc.ndim == 1
    selections = []
    for atomnr, atom in enumerate(backbone_atoms):
        mask = struc["name"] == atom.encode()
        selection = struc[mask]
        if atomnr > 0:
            assert len(selection) == len(selections[0])
        selections.append(selection)
    result = np.empty((len(selections[0]), len(backbone_atoms)), atomic_dtype)
    for snr, s in enumerate(selections):
        result[:, snr] = s
    return result


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
}


def get_sequence(struc: np.ndarray) -> str:
    result = ""
    assert struc.ndim == 1
    res = None
    for a in struc:
        if a["altloc"].decode() not in ("A", " "):
            continue
        curr_res = a["model"], a["chain"], a["icode"], a["resid"]
        if res == curr_res:
            continue
        res = curr_res
        c = a["resname"]
        aa = code.get(c.decode(), "X")
        result += aa
    return result


if __name__ == "__main__":
    import sys

    mmcif_file = sys.argv[1]
    outfile = sys.argv[2]
    auth_chains = bool(int(sys.argv[3]))
    auth_residues = bool(int(sys.argv[4]))
    print(f"auth_chains={auth_chains}, auth_residues={auth_residues}")
    data = parse_mmcif(
        open(mmcif_file).read(), auth_chains=auth_chains, auth_residues=auth_residues
    )
    np.save(outfile, data, allow_pickle=False)
