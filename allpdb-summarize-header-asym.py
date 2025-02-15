"""See allpdb-parse-mmcif.py
Maximum jobs on MBI cluster: 40
"""

POOLSIZE = 30

import json
import os
import sys
import numpy as np
import seamless
from tqdm import tqdm

###seamless.delegate()
seamless.delegate(level=3)  ###

from seamless.highlevel import Checksum
from seamless import transformer
import summarize_header_asym

from seamless.highlevel import Checksum, Buffer

allpdb_keyorder = Checksum.load("allpdb-keyorder.CHECKSUM")
allpdb_keyorder = allpdb_keyorder.resolve("plain")

with open("allpdb-header-index.json") as f:
    allpdb_headers = json.load(f)
# or:
# allpdb_headers = Buffer.load("allpdb-header-index.json").deserialize("plain")


@transformer(return_transformation=True)
def summarize_header_asym_chunk(headers):
    from .summarize_header_asym import summarize_header_asym

    result = {}
    for cifname, header in headers.items():
        header_asym = summarize_header_asym(header)
        result[cifname] = header_asym
        print(cifname)
    return result


summarize_header_asym_chunk.celltypes.headers = "deepcell"
summarize_header_asym_chunk.celltypes.result = "mixed"
summarize_header_asym_chunk.modules.summarize_header_asym = summarize_header_asym

chunksize = 500
key_chunks = [
    allpdb_keyorder[n : n + chunksize]
    for n in range(0, len(allpdb_keyorder), chunksize)
]
for key_chunk in key_chunks:
    key_chunk[:] = [k for k in key_chunk if k in allpdb_headers]
nchunks = len(key_chunks)

with tqdm(total=nchunks, desc="Summarize asym header") as progress_bar:

    def process_chunk(chunk_index):
        key_chunk = [k for k in key_chunks[chunk_index] if k in allpdb_headers]
        header_chunk = {k: allpdb_headers[k] for k in key_chunk}
        return summarize_header_asym_chunk(header_chunk)

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
    print("Collect summarized asym header results...")
    results_cs = [processed_chunk.checksum for processed_chunk in processed_chunks]
    results = {}
    for cs in results_cs:
        chunk_dict = cs.resolve("mixed")
        results.update(chunk_dict)
    print("...done")

    buf = Buffer(results, celltype="mixed")
    buf.save("allpdb-header-summarized-asym")
    buf.checksum.save("allpdb-header-summarized-asym.CHECKSUM")
