
import modesum
import vmodel
import numpy as np

dset    = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=2)
# model.read_axisem_bm('../pysurf/1dmodel_ak135f.txt')
# model.ak135()
# model.trim(zmax=2000.)

model.addlayer(H=200., vsv=7.850000, vsh=7.850000, vpv=13.500000, vph=13.500000, vpf=None, rho=5.75,
                Qp=310., Qs=150., etap=0.0, etas=0.0, frefp=1.0, frefs=1.0, zmin=9999.)

model=model.relayerize(h=2.)


dset.getmodel(inmodel=model)
dset.run_disp(workingdir='/home/leon/code/CPSPy/axisem_benchmark', outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf', nmodes=1, hs=10., 
            period=(np.arange(10, 52, 2)).tolist())
dset.run_eigen(hs=10., infname='/home/leon/code/CPSPy/axisem_benchmark.asdf', outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf')
dset.load('/home/leon/code/CPSPy/axisem_benchmark.asdf')
dset.write_disp('/home/leon/code/CPSPy/axisem_benchmark/homo_ph.ray.txt', wavetype='ray', dtype='ph')
