
import modesum
import vmodel
import numpy as np

dset    = modesum.modesumASDF()
model   = vmodel.Model1d(modelindex=2, earthindex=2)
model.read_axisem_bm('/home/leon/code/axisem/SOLVER/MESHES/prem_aniso_10s/1dmodel_axisem.bm')
# dset.getmodel(inmodel=model)
# dset.run_disp(workingdir='/home/leon/code/CPSPy/axisem_benchmark', outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf', nmodes=1, hs=10., 
#             period=(np.arange(10, 105, 5)).tolist())
# dset.run_eigen(hs=10., infname='/home/leon/code/CPSPy/axisem_benchmark.asdf', outfname='/home/leon/code/CPSPy/axisem_benchmark.asdf')
dset.load('/home/leon/code/CPSPy/axisem_benchmark.asdf')
dset.write_disp('/home/leon/code/CPSPy/axisem_benchmark/gr.disp.txt', dtype='ph')
