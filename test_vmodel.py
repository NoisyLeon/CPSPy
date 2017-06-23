import vmodel
import copy
import matplotlib.pyplot as plt
import numpy as np
model=vmodel.Model1d(modelindex=2)
# model.read('ak135.mod')
model.ak135()
model.write_axisem_bm('tt.bm')
# model2=model.relayerize(h=1.)

# model.read_layer_txt('../SW4Py/cpsin_staircase_1km.txt')
# model.write('./test.mod')
# model1=model.copy()
# model2=model.copy()
# zmin=13
# zmax=128
# model.perturb(dm=.3, zmin=zmin, zmax= zmax)
# model.perturb(dm=.3, zmin=zmin, zmax= zmax, datatype='vp')
# model.perturb(dm=0.3, zmin=zmin, zmax= zmax, datatype='rho')
# model.plotvsak135(zmax=200., datatype='vs')
# model.plotvsak135(zmax=200., datatype='vp')
# model.plotvsak135(zmax=200., datatype='rho')
# model1.perturb(dm=0.1, zmax=21., datatype='vs')
# model2.perturb(dm=-0.1, zmax=21., datatype='vs')
# model.perturb(dm=0.1, zmax=21., datatype='vp')
# model.perturb(dm=0.1, datatype='rho')
# model.write('./disp_dir/ak135_p0.2_20.mod')

#############
# dtype='rho'
# dataArr, depthArr=model.getArr4plot(zmax=400, datatype=dtype)
# dataArr1, depthArr1=model1.getArr4plot(zmax=400., datatype=dtype)
# dataArr2, depthArr2=model2.getArr4plot(zmax=400., datatype=dtype)
# fig, ax= plt.subplots(figsize=(8,12))
# ax, = plt.plot(dataArr, depthArr, 'k-', lw=3 )
# plt.plot(dataArr1, depthArr1, 'b-', lw=10, label='High Vs')
# plt.plot(dataArr2, depthArr2, 'r-', lw=10, label='Low Vs')
# plt.plot(dataArr, depthArr, 'k-', lw=10, label='Reference' )
# plt.ylabel('Depth (km)', fontsize=30)
# plt.xlabel('Vs (km/s)', fontsize=30)
# ax.tick_params(axis='x', labelsize=20)
# ax.tick_params(axis='y', labelsize=20)
# plt.yticks(np.array([0, 20., 35., 50., 100., 150., 200.]))
# plt.legend(loc='lower left', fontsize=25)
# plt.xlim((0,))
# plt.ylim((0, 200))
# plt.gca().invert_yaxis()
# 
# model.plot(zmax=200., datatype='vs')
# model1.plot(zmax=200., datatype='vs')
# model2.plot(zmax=200., datatype='vs')
# plt.show()