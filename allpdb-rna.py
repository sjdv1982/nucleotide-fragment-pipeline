from seamless import Buffer
import numpy as np
from tqdm import tqdm

interface_index, interface_data = Buffer.load(
    "allpdb-filtered-interfaces.mixed"
).deserialize("mixed")
strucs = Buffer.load("intermediate/allpdb-interface-struc.mixed").deserialize("mixed")

rna_strucs = []
rna_struc_index = {}
offset = 0
for code, struc0 in tqdm(strucs.items()):
    if_start, if_size = interface_index[code]
    if not if_size:
        continue
    interfaces = interface_data[if_start : if_start + if_size]
    struc = struc0[struc0["model"] == 1]
    rna_chains = np.unique(interfaces["chain2"])
    for rna_chain in rna_chains:
        code2 = code.split(".")[0] + rna_chain.decode()
        rna_struc = struc[struc["chain"] == rna_chain]
        rna_strucs.append(rna_struc)
        rna_struc_index[code2] = (offset, len(rna_struc))
        offset += len(rna_struc)
rna_strucs = np.concatenate(rna_strucs)
allpdb_rna = rna_struc_index, rna_strucs
buf = Buffer(allpdb_rna, celltype="mixed")
buf.save("intermediate/allpdb-rna.mixed")
buf.get_checksum().save("intermediate/allpdb-rna.mixed.CHECKSUM")
