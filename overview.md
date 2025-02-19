Pipeline:

TODO: allpdb-detect-interfaces.mixed:
move from main to intermediate/  

TODO: add ./config subrepo => nucleotide-fragment-pipeline-config. Will also contain:
    -  .job files
    -  ~/seamless-deployment config (remove Dask!) . Corresponds to step 1 and 1a => adapt seamless-config-step1.sh
    - job .out files (add in .gitignore)
TODO: execution environment: seamless-exact + biopython + opt_einsum

Step 0. Starts with the entire PDB in mmcif format
    Size: about 300 GB
    TODO: input dir that contains all mmcifs, rsynced
    TODO: Seamless script that generates Seamless fairdir from this
    TODO: make a seamless-fair conda package (with the abilities of the current seamless-development conda environment). Remove it from the seamless Docker image? /home/sjoerd/seamless/tests/fairserver/README.txt is messy anyway... and slightly out-of-date (key order is random now)
    TODO: publish fairdir indexes
    TODO: non-Seamless checksum indexing (softlinking) script, that can replace Seamless
    TODO: make the keyorder incremental
    TODO: Seamless config file.
    Result: "allpdb" deepcell and "allpdb-keyorder" (TODO: store the files and add to the repo, not just the checksums)

Step 1: parse all PDBs using Biopython
    Duration:
        MBI cluster: about 1h on one compute node, mmcifs read via NFS, ppdbs already present
            If ppdbs need to be written via hashserver: about 2h
    Purpose: store all PDBs in ppdb ("parsed PDB") format, a Numpy structured dtype.
    Library code: parse_mmcif.py
    Pipeline code: allpdb-parse-mmcif.py
    Pipeline test: allpdb-parse-mmcif-test.py
    TODO: make a blacklist (now hardcoded), Seamless supports that.
    TODO: store ppdb dtype in a code file (now in parse_mmcif.py)
    TODO: make this incremental
    Result: allpdb-struc-index.json deepcell (18 MB, add to the repo) and checksum
        Deepcell size (i.e. all ppdbs): about 120 GB

Step 1a: parse all PDB headers using Biopython
    Purpose: create a dict that can be queried for useful properties, notably:
        stoichiometry and molecule type
        symmetry info
    Code: same as above, but all filenames are [.*]parse_mmcif_header[.*]
    Result: allpdb-header-index.json  deepcell (18 MB, add to the repo) and checksum

Step 1b: summarize all PDB headers
    allpdb-summarize-header.py
    allpdb-summarize-header-asym.py
    Result:
        intermediate/allpdb-header-summarized.json (1.1 GB)
        intermediate/allpdb-header-summarized-asym.json (0.5 GB)

Step 2: Detect interfaces
    allpdb-detect-interfaces.py

Step 3: filter interfaces
    filter-interfaces
    Result: "allpdb-filtered-interfaces.mixed" (4 MB, add to repo) Seamless deepcell containing protein-RNA interfaces

Step 4: Collect and extract RNA
    allpdb-collect-interface-struc.py
    allpdb-rna.py
    Result: allpdb-rna (.mixed) (470 MB, don't add to the repo, only checksum)

Step 5: Add missing atoms using RNA topology
    allpdb-rna-attract.py
    Use ATTRACT aareduce
    TODO: add blacklist
    Result: allpdb-rna-attract (.mixed) (439 MB, don't add to the repo, only checksum)

Step 6-9

For each library (lib = dinuc, trinuc):
For each sequence (trinuc: seq=A/C/G/U \* 3, e.g. ACU. dinuc: seq=A/C/G/U \* 2, e.g. UA)
Or: for each motif
A motif is a mutated sequence (U=>C, G=>A).
    (trinuc: motif=A/C \* 3, e.g. ACC. dinuc: seq=A/C \* 2, e.g. CA)

Step 5: Initial fragment extraction:
    Code: lib-$lib-initial.py
    Result:
        lib-$lib-initial-$seq.npy
        lib-$lib-initial-$seq-origin.txt

Step 6: Mutation
    Code: lib-$lib-initial.py
    Result:
        lib-$lib-mutated-$motif.npy
        lib-$lib-mutated-$motif-origin.txt

Step 7:
    Deredundant (clustering at 0.2 A)
    Code: lib-$lib-deredundant.py
    Library code: clusterlib/, deredundant.py
    Result:
        lib-$lib-nonredundant-$motif.clust

Step 8:
    Filtering

    Applies clustering to fragments
    Then, filtering based on:
    - internal clashes
    - disconnected bonds

    Code: filter-fragments.sh
    Tools code: ...
    Result:
        lib-$lib-nonredundant-filtered-$motif.npy
        lib-$lib-nonredundant-filtered-$motif-origin.txt
