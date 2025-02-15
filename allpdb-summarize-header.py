"""See allpdb-parse-mmcif.py
Maximum jobs on Newton: 10
Good use case to test hashserver
"""

POOLSIZE = 3
# POOLSIZE = 30

import json
import seamless
from tqdm import tqdm

###seamless.delegate()
seamless.delegate(level=3)  ###

from seamless.highlevel import Checksum
from seamless import transformer
import summarize_header

from seamless.highlevel import Checksum, Buffer

allpdb_keyorder = Checksum.load("allpdb-keyorder.CHECKSUM")
allpdb_keyorder = allpdb_keyorder.resolve("plain")

with open("allpdb-header-index.json") as f:
    allpdb_headers = json.load(f)
# or:
# allpdb_headers = Buffer.load("allpdb-header-index.json").deserialize("plain")

header_keys_file = "summarize_header_keys.txt"
header_keys = summarize_header.load_header_keys(header_keys_file)


@transformer(return_transformation=True)
def summarize_header_chunk(headers, header_keys):
    from .summarize_header import summarize_header

    result = {}
    for cifname, header in headers.items():
        summarized_header = summarize_header(header, keys=header_keys)
        result[cifname] = summarized_header
        print(cifname)
    return result


summarize_header_chunk.celltypes.headers = "deepcell"
summarize_header_chunk.celltypes.result = "mixed"
summarize_header_chunk.modules.summarize_header = summarize_header

chunksize = 500
key_chunks = [
    allpdb_keyorder[n : n + chunksize]
    for n in range(0, len(allpdb_keyorder), chunksize)
]
for key_chunk in key_chunks:
    key_chunk[:] = [k for k in key_chunk if k in allpdb_headers]
nchunks = len(key_chunks)

with tqdm(total=nchunks, desc="Summarize header") as progress_bar:

    def process_chunk(chunk_index):
        key_chunk = [k for k in key_chunks[chunk_index] if k in allpdb_headers]
        header_chunk = {k: allpdb_headers[k] for k in key_chunk}
        return summarize_header_chunk(header_chunk, header_keys=header_keys)

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
    print("Collect summarized header results...")
    results_cs = [processed_chunk.checksum for processed_chunk in processed_chunks]
    results = {}
    for cs in results_cs:
        chunk_dict = cs.resolve("mixed")
        results.update(chunk_dict)
    print("...done")

    buf = Buffer(results, celltype="mixed")
    buf.save("allpdb-header-summarized")
    buf.checksum.save("allpdb-header-summarized.CHECKSUM")
