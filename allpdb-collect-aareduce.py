import os
import seamless
from seamless import Buffer, Checksum
from parse_pdb import parse_pdb
from tqdm import tqdm
import numpy as np

outdir = "intermediate/rna-pdbs"
os.makedirs(outdir, exist_ok=True)

seamless.delegate(level=1)

struc_index, _ = Buffer.load("intermediate/allpdb-rna.mixed").deserialize(
    "mixed"
)

print("Load result checksums from file")
result = {}
for code in struc_index:
    outfile = os.path.join(outdir, code + "-aa.pdb.CHECKSUM")
    if os.path.exists(outfile):
        cs = Checksum.load(outfile)
        result[code] = cs

print(f"Result checksums available for {len(result)}/{len(struc_index)} structures")

keys = sorted(result.keys())
chunksize = 100
key_chunks = [keys[n : n + chunksize] for n in range(0, len(keys), chunksize)]
strucs = {}
print("Resolve result PDBs from checksum")
for key_chunk in tqdm(key_chunks):
    checksum_chunk = [result[key] for key in key_chunk]
    chunk_strucs = seamless.multi.resolve(checksum_chunk, nparallel=10, celltype="bytes")
    for key, struc in zip(key_chunk, chunk_strucs):
        strucs[key] = struc

print("Parse result PDBs")
offset = 0
new_strucs = []
new_struc_index = {}
for code in tqdm(strucs):
    pdbtxt2 = strucs[code].decode()
    try:
        new_struc = parse_pdb(pdbtxt2)
    except ValueError:
        continue
    new_strucs.append(new_struc)
    new_struc_index[code] = (offset, len(new_struc))
    offset += len(new_struc)
print(f"Parsing succeeded for {len(new_strucs)}/{len(struc_index)} structures")

new_strucs = np.concatenate(new_strucs)
allpdb_rna_aareduce = new_struc_index, new_strucs
buf = Buffer(allpdb_rna_aareduce, celltype="mixed")
buf.save("intermediate/allpdb-rna-aareduce.mixed")
buf.get_checksum().save("intermediate/allpdb-rna-aareduce.mixed.CHECKSUM")
