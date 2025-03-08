import os
import sys
from nefertiti.functions import superimpose as superimpose_module
from nefertiti.functions.superimpose import superimpose_array

from tqdm import tqdm
import numpy as np

import seamless
from seamless import transformer

SEAMLESS_DELEGATION_LEVEL = int(os.environ["SEAMLESS_DELEGATION_LEVEL"])
POOLSIZE = int(os.environ["POOLSIZE"])

# Must be defined in a config file

import deredundant
from deredundant_seamless_api import detect_pairs
from clusterlib import cluster_from_pairs, write_clustering


seamless.delegate(
    level=SEAMLESS_DELEGATION_LEVEL, raise_exceptions=True
)  ### to run calculations

"""
seamless.delegate(
    level=1, force_database=True, raise_exceptions=True
)  ### to verify calculation. Will not upload input buffer
"""


NREPRESENTATIVES = 3
NSAMPLE = 2000
RUN_TEST = False

diseq = sys.argv[1]
coorf = f"lib-dinuc-mutated-{diseq}.npy"

coor = np.load(coorf)
# coor = coor[:1000]  ###
# NSAMPLE = 500  ###
centered_coor = coor - coor.mean(axis=1)[:, None, :]
coor_residuals0 = np.einsum("ijk,ijk->i", centered_coor, centered_coor)
coor_residuals = np.round(coor_residuals0, decimals=4).astype(float)

np.random.seed(0)
sample_ind = np.random.choice(min(len(coor), 10000), NSAMPLE, replace=False)
sample = centered_coor[sample_ind]


def calc_rmsdmat(sample):
    import numpy as np
    from tqdm import tqdm

    superimpose_array = superimpose_module.superimpose_array

    rmsdmat = np.zeros(
        (len(sample), len(sample)),
    )
    for n in tqdm(range(len(sample))):
        _, rmsd = superimpose_array(sample, sample[n])
        rmsdmat[n] = rmsd
    return rmsdmat


if SEAMLESS_DELEGATION_LEVEL >= 3:
    calc_rmsdmat = transformer(calc_rmsdmat, direct_print=True)
    calc_rmsdmat.modules.superimpose_module = superimpose_module

rmsdmat = calc_rmsdmat(sample)

threshold = 5
while 1:
    pairmat = rmsdmat < threshold
    representatives = []
    # We want 3 representatives with at least 3 neighbours each
    for n in range(NREPRESENTATIVES):
        nb = pairmat.sum(axis=1)
        rep = nb.argmax()
        if nb[rep] <= 3:  # nb[rep][rep] = 1 ...
            break
        representatives.append(rep)
        repmask = pairmat[rep].copy()
        pairmat[repmask] = 0
        pairmat[:, repmask] = 0
    if len(representatives) == NREPRESENTATIVES:
        break
    threshold -= 0.5
    if threshold <= 0:
        raise Exception("FAIL")
repstrucs = sample[representatives]

if RUN_TEST:
    nonrep = np.ones(NSAMPLE, bool)
    nonrep[representatives] = 0
    sample = sample[nonrep]

    true_rmsdmat = rmsdmat[nonrep][:, nonrep]

    # Test matrix reproduction
    fitted_sample = []
    for repstruc in repstrucs:
        fsample_rotmat, _ = superimpose_array(sample, repstruc)
        fsample = np.einsum("ijk,ikl->ijl", sample, fsample_rotmat)
        fitted_sample.append(fsample)

    natoms = coor.shape[1]
    approx_rmsdmat = np.zeros(true_rmsdmat.shape)
    armsd = np.zeros(true_rmsdmat.shape[1])
    for n, struc in tqdm(enumerate(sample), total=len(sample)):
        armsd[:] = 9999
        for fsample in fitted_sample:
            refe = fsample[n]
            d = fsample - refe
            aarmsd = np.sqrt(np.einsum("ijk,ijk->i", d, d) / natoms)
            armsd = np.minimum(armsd, aarmsd)
        approx_rmsdmat[n] = armsd

    approx_rmsdmat = approx_rmsdmat.reshape(-1)
    true_rmsdmat = true_rmsdmat.reshape(-1)

    true_close = (true_rmsdmat < 0.2) & (true_rmsdmat > 0.01)
    print(
        np.histogram(
            approx_rmsdmat[true_close],
            bins=[0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 1.0, 10.0],
        )[0],
        "{:3f}".format(np.max(approx_rmsdmat[true_close])),
        (true_rmsdmat < 0.2).sum(),
        (approx_rmsdmat < 0.4).sum(),
    )
    # /test

    # Test closeness function
    true_rmsdmat = true_rmsdmat.reshape(len(sample), len(sample))
    true_closemat = true_rmsdmat < 0.2
    true_closemat[np.diag_indices(len(true_closemat))] = 0
    sample_residuals = coor_residuals[sample_ind][nonrep]
    pp1, pp2 = deredundant.detect_pairs(sample, sample_residuals, repstrucs)
    print("Positives", true_closemat.sum())
    print("False positives", (true_rmsdmat[pp1, pp2] > 0.2).sum())
    false_negatives = true_closemat
    false_negatives[pp1, pp2] = 0
    false_negatives[pp2, pp1] = 0
    print("False negatives", false_negatives.sum())

    # / test

pp1, pp2 = detect_pairs(centered_coor, coor_residuals, repstrucs, POOLSIZE=POOLSIZE)
pp = np.stack((pp1, pp2), axis=1)

clustering = cluster_from_pairs(pp, len(centered_coor))
write_clustering(f"lib-dinuc-nonredundant-{diseq}.clust", clustering)
