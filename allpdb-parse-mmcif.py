"""
NOTE: essential to set up the fair dirs as both 
hashserver extra dirs (facing the user) 
  and as read buffer folders (facing the cluster workers):

export HASHSERVER_EXTRA_DIRS='/data3/sdevries/fairpdb-mmcif-buffers;/data3/sdevries/fairpdb-mmcif-deepbuffers;/data3/sdevries/fairpdb-mmcif/keyorder;/data3/sdevries/fairpdb-mmcif/access_index'
export SEAMLESS_READ_BUFFER_FOLDERS=$HASHSERVER_EXTRA_DIRS

NOTE: tricky to get poolsize working, poolsize 200 freezes all on the MBI cluster

NOTE: give workers to the hashserver using uvicorn
"""

import os

SEAMLESS_DELEGATION_LEVEL = os.environ["SEAMLESS_DELEGATION_LEVEL"]
# Must be defined in a config file

import seamless
from tqdm import tqdm

from seamless.workflow import Buffer

allpdb = Buffer.load("allpdb-index.json").deserialize("plain")
allpdb_keyorder = Buffer.load("allpdb-keyorder.json").deserialize("plain")

from seamless import transformer
import parse_mmcif


@transformer(return_transformation=True)
def parse_mmcifs(mmcifs):
    from .parse_mmcif import parse_mmcif

    result = {}
    for cifname, cifbuffer in mmcifs.items():
        cifdata = cifbuffer.decode()
        struc = parse_mmcif(cifdata, auth_chains=False, auth_residues=False)
        result[cifname] = struc
    return result


parse_mmcifs.celltypes.mmcifs = "folder"
parse_mmcifs.celltypes.result = "deepcell"
parse_mmcifs.modules.parse_mmcif = parse_mmcif

seamless.delegate(level=SEAMLESS_DELEGATION_LEVEL)

chunksize = 50
key_chunks = [
    allpdb_keyorder[n : n + chunksize]
    for n in range(0, len(allpdb_keyorder), chunksize)
]
nchunks = len(key_chunks)
forbidden = {"1ejg.cif", "4udf.cif"}
for key_chunk in key_chunks:
    key_chunk[:] = [k for k in key_chunk if k not in forbidden]

with tqdm(total=nchunks, desc="Parse mmCIF") as progress_bar:

    def parse_mmcifs_chunk(chunk_index):
        cif_chunk = {k: allpdb[k] for k in key_chunks[chunk_index]}
        return parse_mmcifs(cif_chunk)

    def callback(n, parsed_chunk):
        progress_bar.update(1)
        if parsed_chunk.checksum.value is None:
            print(
                f"""Failure for chunk {n}:
    status: {parsed_chunk.status}
    exception: {parsed_chunk.exception}
    logs: {parsed_chunk.logs}"""
            )

    with seamless.multi.TransformationPool(40) as pool:
        parsed_chunks = pool.apply(parse_mmcifs_chunk, nchunks, callback=callback)

if not any([parsed_chunk.checksum.value is None for parsed_chunk in parsed_chunks]):
    results = [parsed_chunk.checksum for parsed_chunk in parsed_chunks]
    result_index = {}

    """
    # Naive way
    print("Collect parsed mmCIF results...")
    for cs in tqdm(results):
        chunk_dict = cs.resolve("plain")
        result_index.update(chunk_dict)
    print("...done")
    """
    with tqdm(
        total=len(results), desc="Collect parsed mmCIF results..."
    ) as progress_bar:
        for chunk_dict in seamless.multi.resolve(
            results,
            nparallel=50,
            celltype="plain",
            callback=lambda *_: progress_bar.update(1),
        ):
            result_index.update(chunk_dict)

    buf = Buffer(result_index, celltype="plain")
    buf.save("allpdb-struc-index.json")
    buf.checksum.save("allpdb-struc-index.json.CHECKSUM")
