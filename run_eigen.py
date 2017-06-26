import modesum
import vmodel
import numpy as np

dset1   = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=1)
model.ak135()
dset1.getmodel(inmodel=model)
dset1.run_disp(workingdir='/home/lili/code/CPSPy/test_spherical', outfname='/home/lili/code/CPSPy/test_spherical.asdf', nmodes=2,
            period=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

# dset1.run_disp(workingdir='/home/leon/code/CPSPy/test_spherical', outfname='/home/leon/code/CPSPy/test_spherical.asdf', nmodes=2,
#             period=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
# dset1.run_eigen(hs=10., infname='/home/leon/code/CPSPy/test_spherical.asdf', outfname='/home/leon/code/CPSPy/test_spherical.asdf')
# dset1.write_disp('test.txt')

# dset2 = modesum.modesumASDF()
# model=vmodel.Model1d(modelindex=1, earthindex=1)
# model.ak135()
# dset2.getmodel(inmodel=model)
# dset2.run_disp(workingdir='/home/leon/code/CPSPy/test_flat', outfname='/home/leon/code/CPSPy/test_flat.asdf', nmodes=2,
#             period=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
# dset2.run_eigen(hs=0., infname='/home/leon/code/CPSPy/test_flat.asdf', outfname='/home/leon/code/CPSPy/test_flat.asdf')
# 


# dset.run_disp(workingdir='/home/leon/code/CPSPy/test_modesum_iso', outfname='/home/leon/code/CPSPy/eigen002.asdf', nmodes=2, period=[10, 20, 30, 40])
# dset.run_eigen(hs=0., infname='/home/leon/code/CPSPy/eigen002.asdf', outfname='/home/leon/code/CPSPy/eigen002.asdf')
# dset.load('./kerneltest_ti_0.asdf')
# dset.plot_eigen(wavetype='ray', period=10., dtype='dcdav', zmax=120, showfig=False)
# dset.plot_eigen(wavetype='ray', period=20., dtype='dcdav', zmax=120, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=30., dtype='dcdav', zmax=120, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=40., dtype='dcdav', zmax=120, showfig=True, newfig=False)
# # dset.plot_eigen(wavetype='ray', period=10., dtype='dcdah', zmax=120, showfig=False)
# dset.plot_eigen(wavetype='ray', period=20., dtype='dcdah', zmax=120, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=30., dtype='dcdah', zmax=120, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=40., dtype='dcdah', zmax=120, showfig=True, newfig=False)
# dset.plot_eigen(wavetype='ray', period=20., dtype='ur', zmax=80, newfig=False)

