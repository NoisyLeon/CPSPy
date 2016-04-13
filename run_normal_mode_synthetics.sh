#!/bin/bash
sprep96 -HR 0 -HS 10 -M tak135sph.mod -d dfile -NMOD 1 -R
sdisp96
sdpsrf96 -R -PER -XMIN 0.1 -XMAX 100 -ASC
sregn96 
spulse96 -d dfile -D -M 0 -p -l 4 > myf96
f96tosac -B myf96
rm myf96
