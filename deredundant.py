import numpy as np
from nefertiti.functions import superimpose as superimpose_module
from numpy.linalg import det, svd
from tqdm import tqdm


def detect_close(struc1, struc2, representatives, threshold):
    """Detect struc1-struc2 pairs of close structures.
    Close structures are detected by superimposing struc1 and struc2
     onto representative structures, and testing if the coordinates are close.
    This gives an upper estimation of their RMSD, but a good one for low RMSDs.

    struc1, struc2 and representatives must have been centered on the origin
    returns: s1, s2.
    s1 and s2 are arrays of equal length such that (s1[X], s2[X])
     is a closeness pair.

    """
    import opt_einsum
    import numpy as np

    assert struc1.ndim == struc2.ndim == 3
    assert struc1[0].shape == struc2[0].shape == representatives[0].shape
    natoms = struc1.shape[1]
    thresh = (threshold**2) * natoms
    close_mat = np.zeros((len(struc1), len(struc2)), bool)
    for rep in representatives:
        fstruc1_rotmat = superimpose_module.superimpose_array(struc1, rep)[0]
        fstruc1 = np.einsum("ijk,ikl->ijl", struc1, fstruc1_rotmat)
        fstruc2_rotmat = superimpose_module.superimpose_array(struc2, rep)[0]
        fstruc2 = np.einsum("ijk,ikl->ijl", struc2, fstruc2_rotmat)
        d = fstruc1[:, None, :, :] - fstruc2[None, :, :, :]
        sd = np.einsum("ijkl,ijkl->ij", d, d)
        close_mat0 = (sd < thresh).astype(bool)
        close_mat |= close_mat0
    p1, p2 = np.where(close_mat)
    return p1, p2


def verify_close(ind1, ind2, struc1, struc2, residuals1, residuals2, threshold):
    import numpy as np
    from numpy.linalg import det, svd

    assert struc1.ndim == struc2.ndim == 3
    assert struc1[0].shape == struc2[0].shape
    assert len(struc1) == len(residuals1)
    assert len(struc2) == len(residuals2)
    assert residuals1.ndim == residuals2.ndim == 1

    natoms = struc1.shape[1]

    r1 = residuals1[ind1]
    r2 = residuals2[ind2]

    s1 = struc1[ind1]
    s2 = struc2[ind2]
    covar = np.einsum("ijk,ijl->ikl", s1, s2)

    v, s, wt = svd(covar)
    reflect = det(v) * det(wt)
    s[:, -1] *= reflect
    v[:, :, -1] *= reflect[:, None]
    ss = (r1 + r2) - 2 * s.sum(axis=1)

    msd_thresh = (threshold**2) * natoms
    mask = ss < msd_thresh
    return ind1[mask], ind2[mask]


def detect_pairs(struc, residuals, representatives):
    CHUNKSIZE = 1000
    result1 = []
    result2 = []

    offsets = []
    chunks = []
    chunks_residuals = []

    for offset in range(0, len(struc), CHUNKSIZE):
        offsets.append(offset)
        chunk = struc[offset : offset + CHUNKSIZE]
        chunks.append(chunk)
        chunk_residuals = residuals[offset : offset + CHUNKSIZE]
        chunks_residuals.append(chunk_residuals)

    representatives = np.stack(representatives)

    n1n2 = []
    for n1 in range(len(chunks)):
        for n2 in range(n1, len(chunks)):
            n1n2.append((n1, n2))

    for n1, n2 in tqdm(n1n2):
        offset1 = offsets[n1]
        struc1 = chunks[n1]
        residuals1 = chunks_residuals[n1]

        offset2 = offsets[n2]
        struc2 = chunks[n2]
        residuals2 = chunks_residuals[n2]

        p1, p2 = detect_close(struc1, struc2, representatives, 0.4)
        pp1, pp2 = verify_close(p1, p2, struc1, struc2, residuals1, residuals2, 0.2)
        pp1 += offset1
        pp2 += offset2
        mask = pp1 < pp2
        if mask.sum():
            result1.append(pp1[mask])
            result2.append(pp2[mask])

    if len(result1):
        result1 = np.concatenate(result1)
        result2 = np.concatenate(result2)
    return result1, result2
