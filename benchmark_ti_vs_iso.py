import modesum, vmodel
import matplotlib.pyplot as plt
import numpy as np
import os
run =1
per=[10., 15., 20., 25., 30., 35., 40]
dm=-0.05
zmax=20.
wavetype='love'
wavetype='ray'

dti0 = modesum.modesumASDF()
mti  = vmodel.Model1d(modelindex=2)
mti.ak135()
mti.trim(zmax=400.)

dti0.getmodel(h=1, zmax=200., modelindex=2)
# dti0.getmodel(inmodel=mti)

if run:
    dti0.run_disp(workingdir='/home/leon/code/CPSPy/benchmark_ti_iso', outfname='/home/leon/code/CPSPy/benchmark_ti_0.asdf', nmodes=1,
                period=per)
    dti0.run_eigen(infname='/home/leon/code/CPSPy/benchmark_ti_0.asdf', outfname='/home/leon/code/CPSPy/benchmark_ti_0.asdf')
    dti0.del_dir()

    
os.chdir('/home/leon/code/CPSPy')
diso0 = modesum.modesumASDF()
# miso  = vmodel.Model1d(modelindex=1)
# miso.ak135()
# miso.trim(zmax=400.)
diso0.getmodel(h=1, zmax=200., modelindex=1)
# diso0.getmodel(inmodel=miso)
if run:
    diso0.run_disp(workingdir='/home/leon/code/CPSPy/benchmark_ti_iso', outfname='/home/leon/code/CPSPy/benchmark_iso_0.asdf', nmodes=1,
                period=per)
    diso0.run_eigen(infname='/home/leon/code/CPSPy/benchmark_iso_0.asdf', outfname='/home/leon/code/CPSPy/benchmark_iso_0.asdf')
    diso0.del_dir()


mti=dti0.model1d
mti.VsvArr[:10]  =mti.VsvArr[:10]*(1+dm)
mti.VshArr[:10]  = mti.VshArr[:10]*(1+dm)
mti.VpfArr     = np.sqrt ( (mti.VphArr**2 - 2.*((mti.VshArr)**2) ) )
# mti.perturb(zmax=zmax, dm=dm, datatype='vsv')
# mti.perturb(zmax=zmax, dm=dm, datatype='vsh')
miso=diso0.model1d
# miso.perturb(zmax=zmax, dm=dm, datatype='vs')
miso.VsArr[:10]  = miso.VsArr[:10]*(1+dm)


dti1 = modesum.modesumASDF()
dti1.getmodel(inmodel=mti)
if run:
    dti1.run_disp(workingdir='/home/leon/code/CPSPy/benchmark_ti_iso', outfname='/home/leon/code/CPSPy/benchmark_ti_1.asdf', nmodes=1,
                period=per)
    dti1.run_eigen(infname='/home/leon/code/CPSPy/benchmark_ti_1.asdf', outfname='/home/leon/code/CPSPy/benchmark_ti_1.asdf')
    dti1.del_dir()
    
diso1 = modesum.modesumASDF()
diso1.getmodel(inmodel=miso)
if run:
    diso1.run_disp(workingdir='/home/leon/code/CPSPy/benchmark_ti_iso', outfname='/home/leon/code/CPSPy/benchmark_iso_1.asdf', nmodes=1,
                period=per)
    diso1.run_eigen(infname='/home/leon/code/CPSPy/benchmark_iso_1.asdf', outfname='/home/leon/code/CPSPy/benchmark_iso_1.asdf')
    # diso1.del_dir()

dti0.load('/home/leon/code/CPSPy/benchmark_ti_0.asdf')
diso0.load('/home/leon/code/CPSPy/benchmark_iso_0.asdf')
diso1.load('/home/leon/code/CPSPy/benchmark_iso_1.asdf')
dti1.load('/home/leon/code/CPSPy/benchmark_ti_1.asdf')


cti0 = dti0.get_disp_egn(wavetype=wavetype, mode=0, perlst=per)
ciso0 = diso0.get_disp_egn(wavetype=wavetype, mode=0, perlst=per)
ciso1 = diso1.get_disp_egn(wavetype=wavetype, mode=0, perlst=per)
cti1 = dti1.get_disp_egn(wavetype=wavetype, mode=0, perlst=per)


fig, ax=plt.subplots()
plt.plot(per, cti0 , 'o-', lw=3, ms=10, label='ti 0')
plt.plot(per, ciso0, '^--', lw=3, ms=10, label='iso 0')
plt.plot(per, cti1, 'v-', lw=3, ms=10, label='ti 1')
plt.plot(per, ciso1, 'ko--', lw=3, ms=10, label='iso 1')
# plt.legend(numpoints=1, fontsize=20, loc=0)
# ax.tick_params(axis='x', labelsize=20)
# ax.tick_params(axis='y', labelsize=20)
# plt.ylabel('C (km/s)', fontsize=30)
# plt.xlabel('Period (sec)', fontsize=30)

# fig, ax=plt.subplots()
# plt.plot(per, cti1, 'v', ms=10, label='m0 love kernel predicted, montagner')
# plt.plot(per, ciso1, 'ko', ms=10, label='m1')

plt.legend(numpoints=1, fontsize=20, loc=0)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.ylabel('C (km/s)', fontsize=30)
plt.xlabel('Period (sec)', fontsize=30)
plt.show()