import numpy as np

from seamless import transformer, Buffer
from seamless.multi import TransformationPool

from nefertiti.functions import superimpose as superimpose_module

import deredundant

detect_close = transformer(deredundant.detect_close)
detect_close.modules.superimpose_module = superimpose_module

detect_close_async = transformer(deredundant.detect_close, return_transformation=True)
detect_close_async.modules.superimpose_module = superimpose_module

verify_close = transformer(deredundant.verify_close)
verify_close_async = transformer(deredundant.verify_close, return_transformation=True)


@transformer(return_transformation=True)
def detect_and_verify(
    struc1, struc2, offset1, offset2, residuals1, residuals2, representatives
):
    p1, p2 = deredundant.detect_close(struc1, struc2, representatives, 0.4)
    pp1, pp2 = deredundant.verify_close(
        p1, p2, struc1, struc2, residuals1, residuals2, 0.2
    )
    pp1 += offset1
    pp2 += offset2
    mask = pp1 < pp2
    if mask.sum():
        return pp1[mask], pp2[mask]
    else:
        return None, None


detect_and_verify.modules.deredundant = deredundant


def detect_pairs(struc, residuals, representatives, *, POOLSIZE):
    from tqdm import tqdm

    CHUNKSIZE = 1000
    result1 = []
    result2 = []

    offsets = []
    chunks = []
    chunk_checksums = []
    # It is very important to provide the structure chunks as checksum arguments.
    # All of the transformations are constructed immediately, even though only POOLSIZE are executed at any time.
    # If all chunks are serialized and held in memory at the same time, this will explode the memory usage
    chunks_residuals = []

    for offset in range(0, len(struc), CHUNKSIZE):
        # TODO: this gets out-of-memory errors if stuff hasn't been uploaded yet
        offsets.append(offset)
        chunk = struc[offset : offset + CHUNKSIZE]
        chunks.append(chunk)
        buf = Buffer(chunk, "binary")
        buf.upload()
        cs = buf.get_checksum()
        chunk_checksums.append(cs)
        chunk_residuals = residuals[offset : offset + CHUNKSIZE]
        chunks_residuals.append(chunk_residuals)

    representatives = np.stack(representatives)

    n1n2 = []
    for n1 in range(len(chunks)):
        for n2 in range(n1, len(chunks)):
            n1n2.append((n1, n2))

    result1 = []
    result2 = []

    with tqdm(total=len(n1n2)) as progress_bar:

        def detect_close_chunk(n):
            n1, n2 = n1n2[n]
            offset1 = offsets[n1]
            struc_checksum1 = chunk_checksums[n1]
            residuals1 = chunks_residuals[n1]

            offset2 = offsets[n2]
            struc_checksum2 = chunk_checksums[n2]
            residuals2 = chunks_residuals[n2]

            return detect_and_verify(
                struc_checksum1,
                struc_checksum2,
                offset1,
                offset2,
                residuals1,
                residuals2,
                representatives,
            )

        def callback(n, processed_chunk):
            progress_bar.update(1)
            if processed_chunk.checksum.value is None:
                print(
                    f"""Failure for chunk {n}:
        status: {processed_chunk.status}
        exception: {processed_chunk.exception}
        logs: {processed_chunk.logs}"""
                )
                return
            p1, p2 = processed_chunk.value
            if p1 is not None and len(p1):
                result1.append(p1)
                result2.append(p2)

        with TransformationPool(POOLSIZE) as pool:
            pool.apply(detect_close_chunk, len(n1n2), callback=callback)

    if len(result1):
        result1 = np.concatenate(result1)
        result2 = np.concatenate(result2)
    return result1, result2
