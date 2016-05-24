import vmodel
import copy 
model=vmodel.Model1d()
model.ak135()

model.perturb(dm=0.2, zmax=21., datatype='vs')
model.perturb(dm=0.2, zmax=21., datatype='vp')
# model.perturb(dm=0.1, datatype='rho')
model.write('./disp_dir/ak135_p0.2_20.mod')

