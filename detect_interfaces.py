import json
import numpy as np
from collections import namedtuple
import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.transform import Rotation

if not __name__.startswith("seamless_module"):
    # We are building this module from a Seamless Module object
    from nefertiti.functions.superimpose import superimpose


def get_coor(mol):
    mol_coor = np.stack((mol["x"], mol["y"], mol["z"]), axis=-1)
    mol_coor = mol_coor.astype(np.float64)
    # It is very important to do all computations in float64
    # Else you get rounding errors even above 0.1 A,
    #  and in a non-deterministic way, breaking reproducibility
    # Conversely, you must *store* all coordinates in float32
    return mol_coor


def mmult(curr_mat, curr_offset, new_mat, new_offset):
    result_mat = curr_mat.dot(new_mat)
    result_offset = curr_offset.dot(new_mat) + new_offset
    return result_mat, result_offset


def minv(mat, offset):
    result_mat = mat.T
    result_offset = (-offset).dot(mat.T)
    return result_mat, result_offset


def get_bb(coor):
    result = np.empty(6, dtype=coor.dtype)
    cmin = coor.min(axis=0)
    cmax = coor.max(axis=0)
    result[:3] = cmin
    result[3:] = cmax
    return result


def _arrayify_interfaces(interfaces):
    matrix_dtype = np.dtype(
        [("rotation", np.float32, (3, 3)), ("offset", np.float32, (3,))], align=True
    )
    interface_dtype = np.dtype(
        [
            ("chain1", "S4"),
            ("chain2", "S4"),
            ("m1", matrix_dtype),
            ("m2", matrix_dtype),
            ("m_relative", matrix_dtype),
        ],
        align=True,
    )
    result = np.zeros(len(interfaces), dtype=interface_dtype)
    for ifnr, interface in enumerate(interfaces):
        iface = result[ifnr]
        iface["chain1"], iface["chain2"] = interface[:2]
        iface["m1"]["rotation"], iface["m1"]["offset"] = interface[2]
        iface["m2"]["rotation"], iface["m2"]["offset"] = interface[3]
        iface["m_relative"]["rotation"], iface["m_relative"]["offset"] = interface[4]
    for f in "m1", "m2", "m_relative":
        result[f]["rotation"] = (
            np.round(result[f]["rotation"], decimals=3) + 0
        )  # to avoid negative zeroes
        result[f]["offset"] = (
            np.round(result[f]["offset"], decimals=2) + 0
        )  # to avoid negative zeroes
    return result


def detect_interfaces(struc: np.ndarray, header: dict) -> np.ndarray:

    if __name__.startswith("transformer"):
        # This module was built from a Seamless Module object,
        #  and is now being run inside a transformer
        from .superimpose import superimpose as superimpose_func
    else:
        superimpose_func = superimpose

    for k in "struct_ref", "entity_poly", "struct_asym", "pdbx_struct_assembly_gen":
        if k not in header:
            return _arrayify_interfaces([])

    operations = {}
    oper_list = header.get("pdbx_struct_oper_list", {})
    for oper_nr, oper_id in enumerate(oper_list["id"]):
        op = np.array(
            [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]], dtype=float
        )

        for i in range(3):
            op[i][3] = float(oper_list["vector[" + str(i + 1) + "]"][oper_nr])
            for j in range(3):
                op[i][j] = float(
                    oper_list["matrix[" + str(i + 1) + "][" + str(j + 1) + "]"][oper_nr]
                )
        operations[oper_id] = op

    db_codes = {
        k: v
        for k, v in zip(
            header["struct_ref"]["entity_id"], header["struct_ref"]["db_code"]
        )
    }
    poly_entity_types = {
        eid: typ
        for eid, typ in zip(
            header["entity_poly"]["entity_id"], header["entity_poly"]["type"]
        )
    }

    struct_asym = header["struct_asym"]

    ChainCopy = namedtuple("ChainCopy", ("source", "derivative", "rotmat", "offset"))
    chain_copies: list[ChainCopy] = []
    chain_copy_derivatives = {}

    """
    Detect monomeric chain copies

    A chain is defined at the level of the asymmetric unit (header.struct_asym.asym_id),
    this is also where the assembly matrices operate on.
    Two different chains are a copy of each other if:
    - They are instances of the same entity (header.struct_asym.entity_id) 
        OR: they are instances of different entities with the same db code
    - They have at least 80 % of the atoms in common, having the exact same resid, resname and name.
    - The common atoms can be fitted with an RMSD of at most 1 A. 
    Some common atoms can be eliminated to bring the RMSD under 1 A as long as they represent 80 % of the atoms.
    """

    chain_struc = {}
    chain_coors = {}
    chain_ball = {}
    chain_identifiers = {}
    struc_kd_trees = {}
    struc_resids = {}
    for chain_nr, chain in enumerate(struct_asym["id"]):
        assert chain not in chain_coors
        ent_id = struct_asym["entity_id"][chain_nr]
        if ent_id not in poly_entity_types:
            continue
        curr_struc = struc[struc["chain"] == chain.encode()]
        chain_struc[chain] = curr_struc
        coor = get_coor(curr_struc)
        chain_coors[chain] = coor
        com = coor.mean(axis=0)
        cen = coor - com
        dis = np.linalg.norm(cen, axis=1)
        radius = dis.max()
        chain_ball[chain] = com, radius
        struc_kd_tree = KDTree(coor)
        struc_kd_trees[chain] = struc_kd_tree
        struc_resids[chain] = curr_struc["resid"]

    def get_chain_identifiers(chain):
        if chain not in chain_identifiers:
            identifiers = {
                (a["resid"], a["resname"], a["name"]): anr
                for anr, a in enumerate(chain_struc[chain])
            }
            chain_identifiers[chain] = identifiers
        return chain_identifiers[chain]

    for chain1_nr, chain1 in enumerate(struct_asym["id"]):
        if chain1 in chain_copy_derivatives:
            continue
        ent1_id = struct_asym["entity_id"][chain1_nr]
        if ent1_id not in poly_entity_types:
            continue

        atoms1 = chain_struc[chain1]
        atom1_identifiers = get_chain_identifiers(chain1)

        for chain2_nr, chain2 in enumerate(struct_asym["id"]):
            if chain2 in chain_copy_derivatives:
                continue

            if chain2_nr <= chain1_nr:
                continue
            ent2_id = struct_asym["entity_id"][chain2_nr]
            if ent2_id not in poly_entity_types:
                continue
            if poly_entity_types[ent1_id] != poly_entity_types[ent2_id]:
                continue
            is_copy = False
            if ent1_id == ent2_id:
                is_copy = True
            else:
                if ent1_id in db_codes and ent2_id in db_codes:
                    if db_codes[ent1_id] == db_codes[ent2_id]:
                        is_copy = True
            if not is_copy:
                continue
            atoms2 = chain_struc[chain2]
            atom2_identifiers = get_chain_identifiers(chain2)
            common = set(atom1_identifiers.keys()).intersection(
                set(atom2_identifiers.keys())
            )
            biggest_natoms = max(len(atoms1), len(atoms2))
            if len(common) / biggest_natoms < 0.8:
                continue
            atom_pairs = np.array(
                [(atom1_identifiers[idf], atom2_identifiers[idf]) for idf in common]
            )
            atom1_fit0, atom2_fit0 = (
                chain_coors[chain1][atom_pairs[:, 0]],
                chain_coors[chain2][atom_pairs[:, 1]],
            )
            cm1 = atom1_fit0.mean(axis=0)
            cm2 = atom2_fit0.mean(axis=0)
            atom1_fit = atom1_fit0 - cm1
            atom2_fit = atom2_fit0 - cm2
            rotmat_initial, rmsd_initial = superimpose_func(atom2_fit, atom1_fit)
            if rmsd_initial > 0.5:
                atom2_fitted = atom2_fit.dot(rotmat_initial)
                delta = atom2_fitted - atom1_fit
                minimal_common = int(np.ceil(0.8 * biggest_natoms))
                delta_mask = np.argsort((delta * delta).sum(axis=1))[:minimal_common]
                rotmat, rmsd = superimpose_func(
                    atom2_fit[delta_mask], atom1_fit[delta_mask]
                )
            else:
                rotmat = rotmat_initial
                rmsd = rmsd_initial

            if rmsd > 2:
                continue

            # print(chain1, chain2, ent1_id, ent2_id, db_codes.get(ent1_id), db_codes.get(ent2_id), len(atoms1), len(atoms2), len(common), rmsd)
            cm2p = cm2.dot(rotmat)
            offset = -cm2p + cm1
            # atom2_fitted = atom2_fit0.dot(rotmat) + offset
            cop = ChainCopy(chain1, chain2, rotmat, offset)
            chain_copies.append(cop)
            chain_copy_derivatives[chain2] = cop

    """
    Eliminate redundant interfaces

    Consider an interface 1 between chain A and chain B, with M1 applied to A and M2 applied to B
    And an interface 2 between chain A* and chain B*, with M3 applied to A* and M4 applied to B*
    Interface 1 and 2 are redundant if all of the following is true:
    - chain A and A* are the same chain, giving a matrix M5=identity
    - chain A* is a copy of chain A, giving a matrix M5=the superposition of A* on A
    - chain B and B* are the same chain, giving a matrix M6=identity
    - chain B* is a copy of chain B, giving a matrix M6=the superposition of B* on B
    The receptor matrices X1 and X2 are computed as M1 and (M5.T).dot(M3)
    The ligand matrices Y1 and Y2 are computed as M2 and (M6.T).dot(M4)
    The relative ligand matrices R1 and R2 are computed as Y1.dot(X1.T) and Y2.dot(X2.T)
    The difference matrix D is then R1.dot(R2.T)
    If this gives an RMSD of less than 2 A on the ligand (this is a rigid body RMSD!), 
    then the interfaces are the same.
    If not, store X, Y and R for the new interface.
    """

    interfaces = []  # chain1, chain2, X, Y. Chain is original chain, i.e. not a copy.
    for assembly in range(len(header["pdbx_struct_assembly_gen"]["assembly_id"])):
        assembly_chains = [
            c.strip()
            for c in header["pdbx_struct_assembly_gen"]["asym_id_list"][assembly].split(
                ","
            )
        ]

        oper_expressions = header["pdbx_struct_assembly_gen"]["oper_expression"][
            assembly
        ]
        if len(oper_expressions) > 200:
            # skip too massive assemblies
            continue

        for chain1_nr, chain1 in enumerate(struct_asym["id"]):
            if chain1 not in assembly_chains:
                continue

            ent1_id = struct_asym["entity_id"][chain1_nr]
            if ent1_id not in poly_entity_types:
                continue

            coor_rec = chain_coors[chain1]
            bb_rec = get_bb(coor_rec)
            M3_44s = []
            bad_rotation = False
            for oper_expression in oper_expressions:
                M3_44 = np.eye(4)
                for oper in oper_expression:
                    M3_44 = M3_44.dot(operations[oper].T)
                    try:
                        M3_44[:3, :3] = Rotation.from_matrix(M3_44[:3, :3]).as_matrix()
                    except ValueError:
                        bad_rotation = True
                        break
                M3_44s.append(M3_44)

            if bad_rotation:
                continue

            print("HERE", assembly, chain1)
            curr_lig_coors = []
            for M3ind, M3_44 in enumerate(M3_44s):
                M3 = M3_44[:3, :3]
                M3O = M3_44[3, :3]
                M3T, M3TO = minv(M3, M3O)
                # print("ASSEM", chain1, assembly, M3ind, M3, M3O)

                for chain2_nr, chain2 in enumerate(struct_asym["id"]):
                    if chain2_nr < chain1_nr:
                        continue
                    ent2_id = struct_asym["entity_id"][chain2_nr]
                    if ent2_id not in poly_entity_types:
                        continue
                    if chain2 not in assembly_chains:
                        continue

                    # Detect if interface is relevant
                    # i.e. between chains with at least 3 residues per chain involved in contacts

                    lig_coor0 = chain_coors[chain2]

                    M4_44s = []
                    for oper_expression2 in oper_expressions:
                        M4_44 = np.eye(4)
                        for oper in oper_expression2:
                            M4_44 = M4_44.dot(operations[oper].T)
                            M4_44[:3, :3] = Rotation.from_matrix(
                                M4_44[:3, :3]
                            ).as_matrix()
                        M4_44s.append(M4_44)

                    for M4ind, M4_44 in enumerate(M4_44s):
                        if chain2_nr == chain1_nr and M4ind <= M3ind:
                            continue

                        M4 = M4_44[:3, :3]
                        M4O = M4_44[3, :3]

                        Mlig, MligO = mmult(M4, M4O, M3T, M3TO)

                        com1, radius1 = chain_ball[chain1]
                        com2_0, radius2 = chain_ball[chain2]
                        com2 = com2_0.dot(Mlig) + MligO
                        rdis = np.linalg.norm(com2 - com1)
                        if rdis > radius1 + radius2 + 5:
                            continue

                        lig_coor = lig_coor0.dot(Mlig) + MligO
                        lig_bb = get_bb(lig_coor)
                        curr_lig_coors.append((chain2, M3ind, M4ind, lig_bb))

            if not curr_lig_coors:
                continue

            rec_struc_tree = struc_kd_trees[chain1]
            curr_lig_bb = np.stack([coor[3] for coor in curr_lig_coors], axis=0)
            curr_lig_bb_min = curr_lig_bb[:, :3]
            curr_lig_bb_max = curr_lig_bb[:, 3:]
            bb_rec_min = bb_rec[:3]
            bb_rec_max = bb_rec[3:]
            filt1 = (curr_lig_bb_min[:, 0] < bb_rec_max[0]) & (
                bb_rec_min[0] < curr_lig_bb_max[:, 0]
            )
            filt2 = (curr_lig_bb_min[:, 1] < bb_rec_max[1]) & (
                bb_rec_min[1] < curr_lig_bb_max[:, 1]
            )
            filt3 = (curr_lig_bb_min[:, 2] < bb_rec_max[2]) & (
                bb_rec_min[2] < curr_lig_bb_max[:, 2]
            )
            bb_filter = filt1 & filt2 & filt3
            for ind, (chain2, M3ind, M4ind, _) in enumerate(curr_lig_coors):
                if not bb_filter[ind]:
                    continue

                M3_44 = M3_44s[M3ind]
                M3 = M3_44[:3, :3]
                M3O = M3_44[3, :3]
                M3T, M3TO = minv(M3, M3O)

                M4_44 = M4_44s[M4ind]
                M4 = M4_44[:3, :3]
                M4O = M4_44[3, :3]

                lig_coor0 = chain_coors[chain2]
                Mlig, MligO = mmult(M4, M4O, M3T, M3TO)
                lig_coor = lig_coor0.dot(Mlig) + MligO

                lig_struc_tree = KDTree(lig_coor)
                ok = False
                for tree1, tree2, rec_chain in (
                    (rec_struc_tree, lig_struc_tree, chain1),
                    (lig_struc_tree, rec_struc_tree, chain2),
                ):
                    rec_resid = struc_resids[rec_chain]
                    has_contact = np.zeros(len(rec_resid), dtype=bool)
                    contacts = tree1.query_ball_tree(tree2, r=5)
                    has_contact[:] = np.array([bool(l) for l in contacts])
                    ###has_contact[:] = 1 ###
                    if not has_contact.sum():
                        # print("NO CONTACTS")
                        break
                    n_rec_res = len(np.unique(rec_resid[has_contact]))
                    if n_rec_res < 3:
                        # print("TOO FEW RES", n_rec_res, rec_chain)
                        break
                    # print("ENOUGH RES", n_rec_res, rec_chain)
                else:
                    ok = True

                if not ok:
                    # print("SKIP", chain1, chain2)
                    # print()
                    continue

                # / detect if interface is relevant

                orichain1, orichain2 = chain1, chain2
                M5, M5O = np.eye(3), np.zeros(3)
                M6, M6O = np.eye(3), np.zeros(3)

                atoms1 = chain_coors[chain1]
                atoms2 = chain_coors[chain2]
                if chain1 in chain_copy_derivatives:
                    cop1 = chain_copy_derivatives[chain1]
                    M5 = cop1.rotmat
                    M5O = cop1.offset
                    orichain1 = cop1.source
                if chain2 in chain_copy_derivatives:
                    cop2 = chain_copy_derivatives[chain2]
                    M6 = cop2.rotmat
                    M6O = cop2.offset
                    orichain2 = cop2.source
                if len(atoms1) < len(atoms2):
                    orichain1, orichain2 = orichain2, orichain1
                    M5, M5O, M6, M6O = M6, M6O, M5, M5O

                # Compute matrix X2
                # (M5.T).dot(M3)
                M5T, M5TO = minv(M5, M5O)
                X2, X2O = mmult(M5T, M5TO, M3, M3O)

                # Compute matrix Y2
                # (M6.T).dot(M4)
                M6T, M6TO = minv(M6, M6O)
                Y2, Y2O = mmult(M6T, M6TO, M4, M4O)
                # print(X2, X2O)
                # print(Y2, Y2O)

                # Compute matrix R2
                X2T, X2TO = minv(X2, X2O)
                R2, R2O = mmult(Y2, Y2O, X2T, X2TO)
                # print(R2, R2O)

                # print("ADD?", chain1, chain2)

                R2T, R2TO = minv(R2, R2O)
                if orichain1 == orichain2:
                    direc1 = (np.sign(R2O) * R2O * R2O).sum()
                    direc2 = (np.sign(R2TO) * R2TO * R2TO).sum()
                    # print("HOMO", direc1, direc2, R2O, R2TO)
                    MIN_HOMO = 10
                    if direc2 > direc1 + MIN_HOMO:
                        # print("HOMO", direc2-direc1)
                        X2, X2O, Y2, Y2O = Y2, Y2O, X2, X2O
                        X2T, X2TO = minv(X2, X2O)
                        R2, R2O = mmult(Y2, Y2O, X2T, X2TO)
                        R2T, R2TO = minv(R2, R2O)

                # Compare with existing interfaces
                coor_chain2 = chain_coors[orichain2]
                min_delta = None
                for ifnr, interface in enumerate(interfaces):
                    orichain1_0, orichain2_0, (X1, X1O), (Y1, Y1O), (R1, R1O) = (
                        interface
                    )
                    if orichain1_0 != orichain1 or orichain2_0 != orichain2:
                        continue
                    D, DO = mmult(R1, R1O, R2T, R2TO)
                    coor_chain2_D = coor_chain2.dot(D) + DO
                    delta = coor_chain2_D - coor_chain2
                    delta_rmsd = np.sqrt((delta * delta).sum(axis=1).mean())
                    # print("DELTA", ifnr, delta_rmsd)
                    if delta_rmsd < 5:
                        break
                    if min_delta is None or min_delta > delta_rmsd:
                        min_delta = delta_rmsd
                else:
                    # print(len(interfaces)+1, "ADD", chain1, chain2, assembly, M3ind, M4ind, orichain1, orichain2, chain1_nr, chain2_nr, min_delta, R2O, R2)
                    interfaces.append(
                        (orichain1, orichain2, (X2, X2O), (Y2, Y2O), (R2, R2O))
                    )
    return _arrayify_interfaces(interfaces)


if __name__ == "__main__":
    import sys
    from nefertiti.functions.write_pdb import write_pdb

    # example usage: detect-interfaces.py data/2vyc.npy data/2vyc.json data/2vyc-interface.npy --write
    struc = np.load(sys.argv[1])
    with open(sys.argv[2]) as fp:
        header = json.load(fp)
    outfile = sys.argv[3]
    interfaces = detect_interfaces(struc, header)
    print(f"{len(interfaces)} interfaces detected", file=sys.stderr)
    np.save(outfile, interfaces)
    if len(sys.argv) > 4 and sys.argv[4].find("--write") > -1:

        # Write out (relative) interfaces in PDB format
        for ifnr, interface in enumerate(interfaces):
            chain1, chain2, _, _, R = interface
            # print(chain1, chain2)
            rec_struc = struc[struc["chain"] == chain1]
            lig_struc = struc[struc["chain"] == chain2]
            lig_coor00 = get_coor(lig_struc)
            lig_coor = lig_coor00.dot(R[0]) + R[1]
            lig_struc["x"] = lig_coor[:, 0]
            lig_struc["y"] = lig_coor[:, 1]
            lig_struc["z"] = lig_coor[:, 2]
            if chain1 == chain2:
                lig_struc["chain"] = b"Z"
            if_struc = np.concatenate((rec_struc, lig_struc), dtype=rec_struc.dtype)
            pdbdata = write_pdb(if_struc)
            with open(f"interface{ifnr+1}.pdb", "w") as f:
                f.write(pdbdata)

        # Write out absolute interfaces in PDB format
        for ifnr, interface in enumerate(interfaces):
            chain1, chain2, X, Y, R = interface
            rec_struc = struc[struc["chain"] == chain1]
            lig_struc = struc[struc["chain"] == chain2]

            lig_coor00 = get_coor(lig_struc)
            lig_coor = lig_coor00.dot(Y[0]) + Y[1]
            lig_struc["x"] = lig_coor[:, 0]
            lig_struc["y"] = lig_coor[:, 1]
            lig_struc["z"] = lig_coor[:, 2]

            rec_coor0 = get_coor(rec_struc)
            rec_coor = rec_coor0.dot(X[0]) + X[1]
            rec_struc["x"] = rec_coor[:, 0]
            rec_struc["y"] = rec_coor[:, 1]
            rec_struc["z"] = rec_coor[:, 2]

            if chain1 == chain2:
                lig_struc["chain"] = b"Z"
            if_struc = np.concatenate((rec_struc, lig_struc), dtype=rec_struc.dtype)
            pdbdata = write_pdb(if_struc)
            with open(f"abs-interface{ifnr+1}.pdb", "w") as f:
                f.write(pdbdata)
