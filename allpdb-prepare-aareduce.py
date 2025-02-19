import os
import seamless
from seamless import Buffer, Checksum, transformer
import write_pdb, parse_pdb
from tqdm import tqdm

outdir = "intermediate/rna-pdbs"
os.makedirs(outdir, exist_ok=True)

seamless.delegate(level=3)

struc_index, struc_data = Buffer.load("intermediate/allpdb-rna.mixed").deserialize(
    "mixed"
)


@transformer(direct_print=True)
def allpdb_write(struc_index, struc_data):
    from .write_pdb import write_pdb

    try:
        from tqdm import tqdm
    except ImportError:
        tqdm = lambda arg: arg

    result = {}
    for code in tqdm(struc_index):
        start, length = struc_index[code]
        struc = struc_data[start : start + length]
        pdbtxt = write_pdb(struc)
        result[code] = pdbtxt
    return result


allpdb_write.celltypes.result = "folder"
allpdb_write.modules.write_pdb = write_pdb
allpdb_write.modules.parse_pdb = parse_pdb

allpdb_pdb = allpdb_write(struc_index, struc_data)

with open(outdir + "/file.list", "w") as f:
    for k, cs in allpdb_pdb.items():
        outfile = os.path.join(outdir, k + ".pdb.CHECKSUM")
        Checksum(cs).save(outfile)
        print(k, file=f)
