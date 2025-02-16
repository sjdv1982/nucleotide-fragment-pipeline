"""See allpdb-parse-mmcif.py
"""

POOLSIZE = 30

import os

SEAMLESS_DELEGATION_LEVEL = int(os.environ["SEAMLESS_DELEGATION_LEVEL"])
# Must be defined in a config file

import os
import seamless
from tqdm import tqdm

from seamless import Buffer
from seamless import transformer
import summarize_header_asym

allpdb_keyorder = Buffer.load("allpdb-keyorder.json").deserialize("plain")
allpdb_headers = Buffer.load("allpdb-header-index.json").deserialize("plain")


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

seamless.delegate(level=SEAMLESS_DELEGATION_LEVEL, raise_exceptions=True)

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
        chunk_dict = cs.resolve("plain")
        results.update(chunk_dict)
    print("...done")

    buf = Buffer(results, celltype="plain")
    os.makedirs("intermediate", exist_ok=True)
    buf.save("intermediate/allpdb-header-summarized-asym.json")
    buf.checksum.save("allpdb-header-summarized-asym.json.CHECKSUM")
