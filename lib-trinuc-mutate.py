import numpy as np
from mutate import mutate_to_AC
import itertools


bases = ("A", "C", "G", "U")
trinuc_sequences = ["".join(s) for s in itertools.product(bases, repeat=3)]

mut_coor = {}
mut_ori = {}

for seq in trinuc_sequences:
    coorf = f"lib-trinuc-initial-{seq}.npy"
    orif = f"lib-trinuc-initial-{seq}-origin.txt"
    coor = np.load(coorf)
    mut, mutseq = mutate_to_AC(coor, seq)
    if mutseq not in mut_coor:
        mut_coor[mutseq] = []
        mut_ori[mutseq] = []
    mut_coor[mutseq].append(mut)
    with open(orif) as f:
        for l in f.readlines():
            code, first_resid = l.split()
            mut_ori[mutseq].append((code, first_resid, seq))
    print(seq)

print()
for mutseq in mut_coor:
    coorf = f"lib-trinuc-mutated-{mutseq}.npy"
    orif = f"lib-trinuc-mutated-{mutseq}-origin.txt"
    coor = np.concatenate(mut_coor[mutseq])
    ori = mut_ori[mutseq]
    assert len(ori) == len(coor)
    orikeys = [
        (code, first_resid, seq, n) for n, (code, first_resid, seq) in enumerate(ori)
    ]
    orikeys.sort(key=lambda k: int(k[1]))
    orikeys.sort(key=lambda k: k[0])
    inds = [k[3] for k in orikeys]
    np.save(coorf, coor[inds])
    with open(orif, "w") as f:
        for ind in inds:
            code, first_resid, seq = ori[ind]
            print(code, first_resid, seq, file=f)
    print(mutseq)
