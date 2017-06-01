import modesum

# dset = modesum.modesumASDF('test_modesum.asdf')
dset = modesum.modesumASDF()
# dset.getmodel()
# dset.run_disp(workingdir='/home/leon/code/CPSPy/test_modesum_ti', outfname='/home/leon/code/CPSPy/eigen001.asdf', nmodes=2, period=[10, 20, 30, 40])
# dset.run_eigen(hs=0., infname='/home/leon/code/CPSPy/eigen001.asdf', outfname='/home/leon/code/CPSPy/eigen001.asdf')

dset.load('/home/leon/code/CPSPy/eigen001.asdf')
c=dset.plot_eigen(wavetype='ray', period=10., dtype='ur', zmax=80)