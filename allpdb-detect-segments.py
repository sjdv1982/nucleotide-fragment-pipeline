from crocodile.dinuc.from_ppdb import ppdb2nucseq
from seamless import Buffer

rna_struc_index, rna_strucs_data = Buffer.load(
    "intermediate/allpdb-rna-aareduce.mixed"
).deserialize("mixed")

all_segments = []
for code in rna_struc_index:
    start, length = rna_struc_index[code]
    struc = rna_strucs_data[start : start + length]

    monoseq, indices = ppdb2nucseq(
        struc, rna=True, ignore_unknown=True, return_index=True
    )

    segments = []
    seg_first_resid = None
    prev_resid = None
    prev_index = None
    for pos, (ind_first, ind_last) in enumerate(indices):
        curr_resid = struc[ind_first]["resid"]
        if prev_index is None:
            seg_first_resid = curr_resid

        if pos == len(indices) - 1:
            seg_length = curr_resid - seg_first_resid + 1
            segments.append((seg_first_resid, seg_length))
            break

        detect_break = False
        if prev_index is not None:
            if ind_first != prev_index:  # discontiguous atom selection
                detect_break = True
            elif curr_resid != prev_resid + 1:  # discontiguous numbering
                detect_break = True

        if detect_break:
            seg_length = prev_resid - seg_first_resid + 1
            segments.append((seg_first_resid, seg_length))
            seg_first_resid = curr_resid

        prev_index = ind_last
        prev_resid = curr_resid

    for first_resid, length in segments:
        all_segments.append((code, first_resid, length))

with open("allpdb-segments.txt", "w") as f:
    for code, first_resid, length in all_segments:
        print(code, first_resid, length, file=f)
