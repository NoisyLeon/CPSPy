#!/bin/bash
# preparing data
#sprep96 -HR 0 -HS 0 -M tak135sph.mod -d dfile -NMOD 1 -R
sprep96 -HR 0 -HS 0 -M staircase_1km.mod -d dfile -NMOD 1 -R
## compute dispersion curve
sdisp96 
##Get ASCII disersion curve, not required for computing normal mode synthetics
sdpsrf96 -R -PER -XMIN 0.1 -XMAX 100 -ASC
