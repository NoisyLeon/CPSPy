import modesum

dset = modesum.modesumASDF()
dset.getmodel(h=2. , zmax=400.)
# dset.run_disp(workingdir='/home/leon/code/CPSPy/ak135_eigen', outfname='/home/leon/code/CPSPy/eigen_ak135.asdf', nmodes=1,
#             period=[8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50])
# dset.run_disp(workingdir='/home/leon/code/CPSPy/ak135_eigen', outfname='/home/leon/code/CPSPy/eigen_ak135.asdf', nmodes=1,
#             period=[25, 40, 67, 100 ,125])
# dset.run_eigen(infname='/home/leon/code/CPSPy/eigen_ak135.asdf', outfname='/home/leon/code/CPSPy/eigen_ak135.asdf')
dset.load('/home/leon/code/CPSPy/eigen_ak135.asdf')
# dset.plot_eigen(wavetype='ray', period=8., dtype='dcdb', zmax=80, showfig=False)
# dset.plot_eigen(wavetype='ray', period=10., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=12., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=14., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=16., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=18., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=20., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=25., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=30., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=35., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=40., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=45., dtype='dcdb', zmax=80, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=50., dtype='dcdb', zmax=80, newfig=False)

# Yun Wang
# dset.plot_eigen(wavetype='ray', period=25., dtype='dcdb', zmax=380, showfig=False)
# dset.plot_eigen(wavetype='ray', period=40., dtype='dcdb', zmax=380, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=67., dtype='dcdb', zmax=380, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=100., dtype='dcdb', zmax=380, showfig=False, newfig=False)
# dset.plot_eigen(wavetype='ray', period=125., dtype='dcdb', zmax=380, newfig=False)