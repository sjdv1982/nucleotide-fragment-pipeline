from seamless import Buffer

def summarize_header(header, keys):
    result = {}
    for k in keys:
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

def load_header_keys(filename):
    keys = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not len(line):
                continue
            fields = line.split(".")
            keys.append(fields)
    return keys

if __name__ == "__main__":
    import sys
    full_header_file = sys.argv[1]
    header_keys_file = sys.argv[2]
    summarized_header_file = sys.argv[3]
    header_keys = load_header_keys(header_keys_file)
    full_header = Buffer.load(full_header_file).deserialize("plain")
    summarized_header = summarize_header(full_header, keys=header_keys)
    summarized_header_buf = Buffer(summarized_header, "plain")
    summarized_header_buf.save(summarized_header_file)
