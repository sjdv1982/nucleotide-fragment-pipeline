"""See allpdb-parse-mmcif.py

Reading input is the limiting factor.
Poolsize=20 is reasonable if the input files are available via a folder
(SEAMLESS_READ_BUFFER_FOLDERS). Hashserver is usually too slow.`

Runs in about ?? mins on the MBI cluster.
"""

import os

POOLSIZE = int(os.environ["POOLSIZE"])
SEAMLESS_DELEGATION_LEVEL = int(os.environ["SEAMLESS_DELEGATION_LEVEL"])
# Must be defined in a config file

import json
import os
import sys
import numpy as np
import seamless
from tqdm import tqdm

from seamless import Checksum, Buffer
from seamless import transformer
import detect_interfaces
from nefertiti.functions import superimpose

allpdb_keyorder = Buffer.load("allpdb-keyorder.json").deserialize("plain")
allpdb_headers = Buffer.load(
    "intermediate/allpdb-header-summarized-asym.json"
).deserialize("plain")
allpdb_struc = Buffer.load("allpdb-struc-index.json").deserialize("plain")


@transformer(return_transformation=True)
def detect_interfaces_chunk(strucs, headers):
    from .detect_interfaces import detect_interfaces

    result = {}
    for cifname, struc in strucs.items():
        if cifname not in headers:
            continue
        header = headers[cifname]
        interfaces = detect_interfaces(struc, header)
        result[cifname] = interfaces
        print(cifname, len(interfaces))
    return result


detect_interfaces_chunk.celltypes.strucs = "deepcell"
detect_interfaces_chunk.celltypes.headers = "plain"
detect_interfaces_chunk.celltypes.result = "mixed"
detect_interfaces_chunk.modules.detect_interfaces = detect_interfaces
detect_interfaces_chunk.modules.superimpose = superimpose

seamless.delegate(level=SEAMLESS_DELEGATION_LEVEL, raise_exceptions=True)

chunksize = 30
key_chunks = [
    allpdb_keyorder[n : n + chunksize]
    for n in range(0, len(allpdb_keyorder), chunksize)
]
for key_chunk in key_chunks:
    key_chunk[:] = [k for k in key_chunk if k in allpdb_headers and k in allpdb_struc]
nchunks = len(key_chunks)

with tqdm(total=nchunks, desc="Detect interfaces") as progress_bar:

    def process_chunk(chunk_index):
        key_chunk = [
            k
            for k in key_chunks[chunk_index]
            if k in allpdb_headers and k in allpdb_struc
        ]
        struc_chunk = {k: allpdb_struc[k] for k in key_chunk}
        header_chunk = {k: allpdb_headers[k] for k in key_chunk}
        return detect_interfaces_chunk(struc_chunk, header_chunk)

    def callback(n, processed_chunk):
        progress_bar.update(1)
        if processed_chunk.checksum.value is None:
            print(
                f"""Failure for chunk {n}:
    status: {processed_chunk.status}
    exception: {processed_chunk.exception}
    logs: {processed_chunk.logs}"""
            )

    with seamless.multi.TransformationPool(POOLSIZE) as pool:
        processed_chunks = pool.apply(process_chunk, nchunks, callback=callback)

if not any(
    [processed_chunk.checksum.value is None for processed_chunk in processed_chunks]
):
    results_cs = [processed_chunk.checksum for processed_chunk in processed_chunks]
    results = {}

    """
    # Naive way
    for cs in tqdm(results_cs, desc="Collect detected interface results..."):
        chunk_dict = cs.resolve("mixed")
        results.update(chunk_dict)
    """

    with tqdm(
        total=len(processed_chunks), desc="Collect detected interface results..."
    ) as progress_bar:
        for chunk_dict in seamless.multi.resolve(
            results_cs,
            nparallel=50,
            celltype="mixed",
            callback=lambda *_: progress_bar.update(1),
        ):
            results.update(chunk_dict)

    print("...done")

    offsets = {}
    offset = 0
    interfaces_concat = []
    for k in allpdb_keyorder:
        if k not in results:
            continue
        ifaces = results[k]
        if ifaces.size == 0:
            offsets[k] = (0, 0)
            continue
        offsets[k] = (offset, len(ifaces))
        interfaces_concat.append(ifaces)
        offset += len(ifaces)

    if interfaces_concat:
        interfaces_concat = np.concatenate(
            interfaces_concat, dtype=interfaces_concat[0].dtype
        )
    allpdb_interfaces = (offsets, interfaces_concat)

    buf = Buffer(allpdb_interfaces, celltype="mixed")
    buf.save("intermediate/allpdb-interfaces.mixed")
    buf.get_checksum().save("allpdb-interfaces.mixed.CHECKSUM")
