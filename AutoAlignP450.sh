#!/bin/bash

export CUDA_VISIBLE_DEVICES="1"

cd L

### load fasta
list=$(ls *.pdb)

###
cd ../


### run command
for f in ${list}; do
    /home/afonoda/seaborn/colabfold/colabfold-conda/bin/python3.7 ../../../../PDB_P450align.py -i L/${f}
done

