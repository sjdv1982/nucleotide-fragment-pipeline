import numpy as np
from mutate import mutate, mutate_to_AC

seq = "GCU"  # => ACC
tmpl = np.load(f"templates/{seq}-template-ppdb.npy")
coorf = f"lib-trinuc-initial-{seq}.npy"

coor = np.load(coorf)
coor = coor[:1000]

coor_mut, acc = mutate_to_AC(coor, seq)
assert acc == "ACC", acc

recoor = mutate(coor_mut, acc, seq)

recoor_mut, acc = mutate_to_AC(recoor, seq)
assert acc == "ACC", acc

print(np.allclose(coor_mut, recoor_mut))
