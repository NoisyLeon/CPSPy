import modesum
import vmodel

dset    = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=2)
model.read_axisem_bm('/home/leon/code/axisem/SOLVER/MESHES/prem_aniso_10s/1dmodel_axisem.bm')
# model.trim(zmax=1000.)
dset.getmodel(inmodel=model)
dset.run_disp(workingdir='/home/leon/code/CPSPy/axisem_benchmark', outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf', nmodes=1)
dset.run_eigen(hs=10., runtype='SYN', infname='/home/leon/code/CPSPy/axisem_benchmark.asdf',
            outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf', run=True)
dset.run_pulse(infname='/home/leon/code/CPSPy/axisem_benchmark.asdf', run=False, dist0=1000, Nd=1)



