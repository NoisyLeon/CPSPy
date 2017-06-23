import modesum
import vmodel

dset    = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=2)
model.ak135()
dset.getmodel(inmodel=model)
dset.run_disp(workingdir='/home/leon/code/CPSPy/test_synthetic_spherical', outfname='/home/leon/code/CPSPy/test_synthetic_spherical.asdf', nmodes=1)
dset.run_eigen(hs=10., runtype='SYN', infname='/home/leon/code/CPSPy/test_synthetic_spherical.asdf',
            outfname='/home/leon/code/CPSPy/test_synthetic_spherical.asdf')
dset.run_pulse(infname='/home/leon/code/CPSPy/test_synthetic_spherical.asdf', run=True)



