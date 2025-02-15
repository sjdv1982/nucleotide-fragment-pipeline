# First run ./test-seamless-startup.sh

import os
import sys
import seamless
seamless.delegate(level=1)

from seamless.highlevel import Checksum
from seamless import transformer
import detect_interfaces
from nefertiti.functions import superimpose

data_cs = Checksum.load("data.CHECKSUM")
data_index = data_cs.resolve("plain")

test_strucs = {}
test_headers = {}
for k in list(data_index.keys()):
    k0 = os.path.splitext(k)[0]
    if len(k0) > 4:
        continue
    k2 = k0 + ".cif"
    if k.endswith(".json"):
        test_headers[k2] = data_index[k]
    elif k.endswith(".npy"):
        test_strucs[k2] = data_index[k]

for k, cs0 in test_headers.items():
    cs = Checksum(cs0)
    print(k, cs.resolve("text")[:50])
for k, cs0 in test_strucs.items():
    cs = Checksum(cs0)
    print(k, cs.resolve("binary")[0])

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
        print(cifname, interfaces)
    return result
detect_interfaces_chunk.celltypes.strucs = "deepcell"
detect_interfaces_chunk.celltypes.headers = "deepcell"
detect_interfaces_chunk.celltypes.result = "deepcell"
detect_interfaces_chunk.modules.detect_interfaces = detect_interfaces
detect_interfaces_chunk.modules.superimpose = superimpose
detect_interfaces_chunk.direct_print = True

tf = detect_interfaces_chunk(test_strucs, test_headers)
tf.compute()
print(tf.logs)
print("Exception:", tf.exception)
if tf.checksum.value is None:
    sys.exit(1)
result_cs = Checksum(tf.checksum)
result_index = result_cs.resolve("plain")
print(result_index)
for k, v in result_index.items():
    struc = Checksum(v).resolve("mixed")
    print(k, struc[0])