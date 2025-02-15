#!/bin/bash
set -u -e

lib=$1 # dinuc or trinuc
motif=$2 # e.g. AA or ACA
python3 apply_clustering.py lib-$lib-mutated-$motif.npy lib-$lib-nonredundant-$motif.clust lib-$lib-nonredundant-$motif.npy
python3 apply_clustering.py lib-$lib-mutated-$motif-origin.txt lib-$lib-nonredundant-$motif.clust lib-$lib-nonredundant-$motif-origin.txt --concat

echo 'Discard clashing fragments' > /dev/stderr
python discard_clashing_fragments.py lib-$lib-nonredundant-$motif.npy templates/$motif-template-ppdb.npy lib-$lib-nonredundant-$motif.nonclash-filter.npy

echo 'Discard disconnected fragments' > /dev/stderr
python discard_disconnected_fragments.py lib-$lib-nonredundant-$motif.npy templates/$motif-template-ppdb.npy lib-$lib-nonredundant-$motif.connected-filter.npy
python combine-filters.py \
    lib-$lib-nonredundant-$motif.nonclash-filter.npy \
    lib-$lib-nonredundant-$motif.connected-filter.npy \
    lib-$lib-nonredundant-$motif.filter.npy
python apply_filter.py lib-$lib-nonredundant-$motif.npy lib-$lib-nonredundant-$motif.filter.npy lib-$lib-nonredundant-filtered-$motif.npy
python apply_filter.py lib-$lib-nonredundant-$motif-origin.txt lib-$lib-nonredundant-$motif.filter.npy lib-$lib-nonredundant-filtered-$motif-origin.txt
