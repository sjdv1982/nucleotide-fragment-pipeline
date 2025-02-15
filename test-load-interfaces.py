import seamless

seamless.delegate(level=1)
from seamless import Checksum

interfaces_cs = Checksum.load("allpdb-interfaces.CHECKSUM")
interfaces_tuple = interfaces_cs.resolve("mixed")
offsets, interfaces = interfaces_tuple

start, length = offsets["1b7f.cif"]
ifaces = interfaces[start : start + length]
print(len(ifaces))
print(ifaces)
