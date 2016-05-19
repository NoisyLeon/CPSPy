#!/bin/bash
# preparing data
#sprep96 -HR 0 -HS 1. -M tak135sph.mod -d dfile -NMOD 2 -R
sprep96 -HR 0 -HS 1. -M test.mod -d dfile -NMOD 2 -R
# compute dispersion curve
sdisp96
#Get ASCII disersion curve, not required for computing normal mode synthetics
sdpsrf96 -R -PER -XMIN 0.1 -XMAX 100 -ASC -TXT
# compute eigenfunctions
sregn96 -NOQ
# compute normal mode synthetics -M 0 fundamental mode only
spulse96 -d dfile -D -p -l 4 > myf96
f96tosac -B myf96
rm myf96
mkdir sac_dir_test_mod
mv *sac sac_dir_test_mod
