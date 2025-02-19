mkdir -p .result
cd data
for i in *.cif; do  
    code=${i%%.*}
    python ../../detect_interfaces.py ../result/$code.npy ../result/$code-header-asym.json ../result/$code-interface.npy
    echo $i
done
