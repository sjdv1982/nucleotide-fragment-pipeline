mkdir -p result
cd data
echo parse mmcif
for i in *.cif; do  
    python3 ../../parse_mmcif.py $i  ../result/${i%%.*}.npy 0 0
    echo $i
done
echo parse header
for i in *.cif; do  
    python3 ../../parse_mmcif_header.py $i  ../result/${i%%.*}-header.json
    echo $i
done
echo summarize header asym
for i in *.cif; do  
    python3 ../../summarize_header_asym.py ../result/${i%%.*}-header.json  ../result/${i%%.*}-header-asym.json
    echo $i
done
