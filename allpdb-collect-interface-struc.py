import numpy as np
from seamless import Buffer, Checksum
from tqdm import tqdm

import seamless
seamless.delegate(level=2)
interface_index, interface_data = Buffer.load("allpdb-filtered-interfaces").deserialize("mixed")
struc_index = Buffer.load("allpdb-struc-index.json").deserialize("plain")
allpdb_keyorder = Checksum.load("allpdb-keyorder.CHECKSUM")
allpdb_keyorder = allpdb_keyorder.resolve("plain")

all_keys = []
for key in allpdb_keyorder:
    if key in interface_index:
        all_keys.append(key)

chunksize = 100
key_chunks = [all_keys[n:n+chunksize] for n in range(0, len(all_keys), chunksize)]

result = {}
for key_chunk in tqdm(key_chunks):
    checksum_chunk = [struc_index[key] for key in key_chunk]
    strucs = seamless.multi.resolve(checksum_chunk, nparallel=10, celltype="binary")
    for key, struc in zip(key_chunk, strucs):
        ifsize, nif = interface_index[key]
        interfaces = interface_data[ifsize:ifsize+nif]
        chains = set()
        for interface in interfaces:
            chains.add(interface["chain1"])
            chains.add(interface["chain2"])
        filtered_struc_mask = np.isin(struc["chain"], list(chains))
        filtered_struc = struc[filtered_struc_mask]
        #print(key, len(struc), len(filtered_struc))
        result[key] = filtered_struc

buf=Buffer(result, celltype="mixed")
buf.save("allpdb-interface-struc")
buf.checksum.save("allpdb-interface-struc.CHECKSUM")