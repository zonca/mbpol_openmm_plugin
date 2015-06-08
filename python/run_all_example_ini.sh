#!/bin/bash
for f in example_ini/*.ini
do
    bin/mbpol_builder $f $f.py
done
for f in example_ini/*.ini
do
    echo $f.py
    python $f.py
done
