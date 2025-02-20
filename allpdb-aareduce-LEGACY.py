import os
import random
from tqdm import tqdm
import numpy as np
from seamless import Buffer
from write_pdb import write_pdb
from parse_pdb import parse_pdb
import subprocess

rna_struc_index, rna_strucs_data = Buffer.load(
    "intermediate/allpdb-rna.mixed"
).deserialize("mixed")

new_rna_strucs = []
new_rna_struc_index = {}
offset = 0
keys = list(rna_struc_index.keys())
random.shuffle(keys)
for code in tqdm(keys):
    start, length = rna_struc_index[code]
    rna_struc = rna_strucs_data[start : start + length]

    pdbtxt = write_pdb(rna_struc)
    tmpf = f"/tmp/tmp-{code}.pdb"
    with open(tmpf, "w") as f:
        f.write(pdbtxt)
    try:
        subprocess.check_call(
            f"python $ATTRACTDIR/../allatom/aareduce.py {tmpf} --rna --nalib --heavy",
            shell=True,
        )
    except subprocess.CalledProcessError:
        print("MISSING", code)
        continue
    finally:
        os.remove(tmpf)
    tmpf2 = f"/tmp/tmp-{code}-aa.pdb"
    try:
        with open(tmpf2) as f:
            pdbtxt2 = f.read()
    except Exception:
        print("MISSING", code)
        continue
    if not len(pdbtxt2):
        print("MISSING", code)
        continue

    new_rna_struc = parse_pdb(pdbtxt2)

    new_rna_strucs.append(new_rna_struc)
    new_rna_struc_index[code] = (offset, len(new_rna_struc))
    offset += len(new_rna_struc)

new_rna_strucs = np.concatenate(new_rna_strucs)
allpdb_rna_attract = new_rna_struc_index, new_rna_strucs
buf = Buffer(allpdb_rna_attract, celltype="mixed")
buf.save("intermediate/allpdb-rna-attract")
buf.checksum.save("intermediate/allpdb-rna-attract.CHECKSUM")
