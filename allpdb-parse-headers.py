"""See allpdb-parse-mmcif.py
Cannot run simultaneously with that script on the MBI cluster:
40 jobs loading data is apparently too much.
"""

import os

SEAMLESS_DELEGATION_LEVEL = int(os.environ["SEAMLESS_DELEGATION_LEVEL"])
# Must be defined in a config file

import seamless
from tqdm import tqdm

import parse_mmcif_header

from seamless import Buffer

allpdb = Buffer.load("allpdb-index.json").deserialize("plain")
allpdb_keyorder = Buffer.load("allpdb-keyorder.json").deserialize("plain")

from seamless import transformer


@transformer(return_transformation=True)
def parse_mmcif_headers(mmcifs):
    from .parse_mmcif_header import parse_mmcif_header

    result = {}
    for cifname, cifbuffer in mmcifs.items():
        cifdata = cifbuffer.decode()
        header = parse_mmcif_header(cifdata)
        result[cifname] = header
    return result


parse_mmcif_headers.celltypes.mmcifs = "folder"
parse_mmcif_headers.celltypes.result = "deepcell"
parse_mmcif_headers.modules.parse_mmcif_header = parse_mmcif_header

seamless.delegate(level=SEAMLESS_DELEGATION_LEVEL)

chunksize = 50
key_chunks = [
    allpdb_keyorder[n : n + chunksize]
    for n in range(0, len(allpdb_keyorder), chunksize)
]
nchunks = len(key_chunks)

with tqdm(total=nchunks, desc="Parse mmCIF headers") as progress_bar:

    def parse_mmcif_headers_chunk(chunk_index):
        cif_chunk = {k: allpdb[k] for k in key_chunks[chunk_index]}
        return parse_mmcif_headers(cif_chunk)

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
        parsed_chunks = pool.apply(
            parse_mmcif_headers_chunk, nchunks, callback=callback
        )

if not any([parsed_chunk.checksum.value is None for parsed_chunk in parsed_chunks]):
    results = [parsed_chunk.checksum for parsed_chunk in parsed_chunks]
    result_index = {}

    """
    #Naive way
    print("Collect parsed mmCIF header results...")
    for cs in results:
        chunk_dict = cs.resolve("plain")
        result_index.update(chunk_dict)
    """

    with tqdm(
        total=len(results), desc="Collect parsed mmCIF header results..."
    ) as progress_bar:
        for chunk_dict in seamless.multi.resolve(
            results,
            nparallel=50,
            celltype="plain",
            callback=lambda *_: progress_bar.update(1),
        ):
            result_index.update(chunk_dict)

    print("...done")

    buf = Buffer(result_index, celltype="plain")
    buf.save("allpdb-header-index.json")
    buf.checksum.save("allpdb-header-index.json.CHECKSUM")
