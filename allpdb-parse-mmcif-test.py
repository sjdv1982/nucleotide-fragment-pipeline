# First run ./test-seamless-startup.sh

import sys
import seamless

seamless.delegate(level=1)

from seamless.highlevel import Checksum

data_cs = Checksum.load("data.CHECKSUM")
data_index = data_cs.resolve("plain")
for k in list(data_index.keys()):
    if not k.endswith(".cif"):
        data_index.pop(k)

for k, cs0 in data_index.items():
    cs = Checksum(cs0)
    print(k, cs.resolve("text")[:50])

from seamless import transformer
from nefertiti.functions import parse_mmcif


@transformer(return_transformation=True)
def parse_mmcifs(mmcifs):
    from . import parse_mmcif

    result = {}
    for cifname, cifbuffer in mmcifs.items():
        print(cifname)
        cifdata = cifbuffer.decode()
        struc = parse_mmcif.parse_mmcif(cifdata, auth_chains=False, auth_residues=False)
        result[cifname] = struc
        print(cifname)
    return result


parse_mmcifs.celltypes.mmcifs = "folder"
parse_mmcifs.celltypes.result = "deepcell"
parse_mmcifs.modules.parse_mmcif = parse_mmcif
parse_mmcifs.direct_print = True

tf = parse_mmcifs(data_index)
tf.compute()
print(tf.logs)
print("Exception:", tf.exception)
if tf.checksum.value is None:
    sys.exit(1)
result_cs = Checksum(tf.checksum)
result_index = result_cs.resolve("plain")
for k, v in result_index.items():
    struc = Checksum(v).resolve("binary")
    print(k, struc[0])
