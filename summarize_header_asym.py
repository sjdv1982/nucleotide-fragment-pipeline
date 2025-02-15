from seamless.core.protocol.json import json_dumps

def summarize_header_asym(header):
    keys = (
        "pdbx_struct_assembly_gen",
        "pdbx_struct_oper_list",
        ("struct_ref", "entity_id"),
        ("struct_ref", "db_code"),
        ("entity_poly","entity_id"),
        ("entity_poly","type"),
        ("struct_asym","id"),
        ("struct_asym","entity_id"),
    )
    result = {}
    for knr, k in enumerate(keys):
        if isinstance(k, str):
            k = (k,)
        sub = header
        target = result
        for kk in k[:-1]:
            if kk not in sub:
                continue
            if kk not in target:
                target[kk] = {}
            target = target[kk]
            sub = sub[kk]
        kk = k[-1]
        if kk in sub:
            target[kk] = sub[kk]
    return result

if __name__ == "__main__":
    import sys
    import json
    full_header_file = sys.argv[1]
    asym_header_file = sys.argv[2]
    with open(full_header_file) as f:
        full_header = json.load(f)
    asym_header = summarize_header_asym(full_header)
    with open(asym_header_file, "w") as f:
        f.write(json_dumps(asym_header) + "\n")