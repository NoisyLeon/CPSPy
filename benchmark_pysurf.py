import modesum
import vmodel
import numpy as np

dset   = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=1)
model.ak135()
dset.getmodel(inmodel=model)
dset.run_disp(workingdir='/home/lili/code/CPSPy/benchmark_pysurf', outfname='/home/lili/code/CPSPy/benchmark_pysurf.asdf', nmodes=2,
            period=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
dset.run_eigen(hs=0., infname='/home/lili/code/CPSPy/benchmark_pysurf.asdf', outfname='/home/lili/code/CPSPy/benchmark_pysurf.asdf')

dset.write_disp(outfname='/home/lili/code/pysurf/benchmark_cps/ak135_cps_ray_ph.txt')
dset.write_disp(outfname='/home/lili/code/pysurf/benchmark_cps/ak135_cps_ray_gr.txt', dtype='U')
dset.write_disp(outfname='/home/lili/code/pysurf/benchmark_cps/ak135_cps_love_ph.txt', wavetype='love')
dset.write_disp(outfname='/home/lili/code/pysurf/benchmark_cps/ak135_cps_love_gr.txt',wavetype='love', dtype='U')
# dset.write_disp(outfname='../pysurf/benckmark_cps/ak135_cps_ray_ph.txt')