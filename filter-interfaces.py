from seamless import Buffer
from tqdm import tqdm
import numpy as np


def mmult(curr_mat, curr_offset, new_mat, new_offset):
    result_mat = curr_mat.dot(new_mat)
    result_offset = curr_offset.dot(new_mat) + new_offset
    return result_mat, result_offset


def minv(mat, offset):
    result_mat = mat.T
    result_offset = (-offset).dot(mat.T)
    return result_mat, result_offset


def filter_interfaces(interfaces, header, na_type):
    h = header
    chains_to_entities = {
        k: v for k, v in zip(h["struct_asym"]["id"], h["struct_asym"]["entity_id"])
    }
    moltype = {
        eid: typ
        for eid, typ in zip(h["entity_poly"]["entity_id"], h["entity_poly"]["type"])
    }
    molweight = {
        eid: w for eid, w in zip(h["entity"]["id"], h["entity"]["formula_weight"])
    }

    result = []
    for interface in interfaces:
        chain1, chain2 = interface["chain1"].decode(), interface["chain2"].decode()
        if chain1 not in chains_to_entities or chain2 not in chains_to_entities:
            continue
        ent1 = chains_to_entities[chain1]
        ent2 = chains_to_entities[chain2]

        if moltype[ent1] == "polypeptide(L)" and moltype[ent2] == na_type:
            ent_prot, ent_na = ent1, ent2
            chain_prot, chain_na = chain1, chain2
            mat_prot, mat_na = interface["m1"], interface["m2"]
        elif moltype[ent2] == "polypeptide(L)" and moltype[ent1] == na_type:
            ent_prot, ent_na = ent2, ent1
            chain_prot, chain_na = chain2, chain1
            mat_prot, mat_na = interface["m2"], interface["m1"]
        else:
            continue

        if float(molweight[ent_na]) > 50000:
            continue

        iface = interface.copy()
        iface["chain1"] = chain_prot.encode()
        iface["chain2"] = chain_na.encode()
        iface["m1"]["rotation"] = mat_prot[0]
        iface["m1"]["offset"] = mat_prot[1]
        iface["m2"]["rotation"] = mat_na[0]
        iface["m2"]["offset"] = mat_na[1]
        X2, X2O = mat_prot
        Y2, Y2O = mat_na
        X2T, X2TO = minv(X2, X2O)
        R2, R2O = mmult(Y2, Y2O, X2T, X2TO)
        iface["m_relative"]["rotation"] = R2
        iface["m_relative"]["offset"] = R2O
        result.append(iface)
    return result


print("Load all interfaces")
interfaces_index, interfaces_data = Buffer.load(
    "intermediate/allpdb-interfaces.mixed"
).deserialize("mixed")
print("Load all headers...")
headers = Buffer.load("intermediate/allpdb-header-summarized.json").deserialize("plain")
print("..loaded")

na_type = "polyribonucleotide"  # for DNA: 'polydeoxyribonucleotide'

filtered_interfaces = []
filtered_interfaces_index = {}
offset = 0
for code in tqdm(interfaces_index):
    start, size = interfaces_index[code]
    if size == 0:
        continue
    curr_interfaces = interfaces_data[start : start + size]
    if code not in headers:
        continue
    curr_header = headers[code]
    ifaces = filter_interfaces(curr_interfaces, curr_header, na_type=na_type)
    if not len(ifaces):
        continue
    filtered_interfaces_index[code] = (offset, len(ifaces))
    filtered_interfaces.extend(ifaces)
    offset += len(ifaces)

if filtered_interfaces:
    filtered_interfaces = np.array(
        filtered_interfaces, dtype=filtered_interfaces[0].dtype
    )

allpdb_filtered_interfaces = (filtered_interfaces_index, filtered_interfaces)

buf = Buffer(allpdb_filtered_interfaces, celltype="mixed")
buf.save("allpdb-filtered-interfaces.mixed")
buf.get_checksum().save("allpdb-filtered-interfaces.mixed.CHECKSUM")
