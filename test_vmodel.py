import vmodel
import copy 
model=vmodel.Model1d()
model.ak135()
model1=copy.deepcopy(model)


model.perturb(dm=0.5, datatype='vs')
model.perturb(dm=0.5, datatype='vp')
# model.perturb(dm=0.1, datatype='rho')
model.write('./disp_dir/ak135_vsvp_+0.5.mod')

