#/bin/bash
tprep96 -M my.mod -DT 0.1 -NPTS 16384 -NMOD 1 -HS 0 -HR 0 -L -R 
tdisp96 
tdpsrf96 -U -C -G -PER -XMIN 1 -XMAX 100 -TXT -R
tregn96 -HS 10 -HR 0 -NOQ
tpulse96 -d dfile -EX -D -i -FUND > tempf96 
f96tosac -B tempf96
#rm tdisp*
#mv TDISPR.TXT TDISPR.TXT_1
#rm TDISP*PLT
#rm *egn
#rm tempf96
#mv B00109ZEX.sac m1.sac
#rm B*sac
