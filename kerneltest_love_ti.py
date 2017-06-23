import modesum, vmodel
import matplotlib.pyplot as plt
import numpy as np

run =1
per=[10., 15., 20., 25., 30., 35., 40]
datatype='vs'
dm=-0.05
zmax=10.
wavetype='Love'
wavetype='Rayleigh'

dset0 = modesum.modesumASDF()
dset0.getmodel(h=1, zmax=200., modelindex=2)
if run:
    dset0.run_disp(workingdir='/home/leon/code/CPSPy/kerneltest_ti', outfname='/home/leon/code/CPSPy/kerneltest_ti_0.asdf', nmodes=1,
                period=per)
    dset0.run_eigen(infname='/home/leon/code/CPSPy/kerneltest_ti_0.asdf', outfname='/home/leon/code/CPSPy/kerneltest_ti_0.asdf')
dset0.load('/home/leon/code/CPSPy/kerneltest_ti_0.asdf')

dset3 = dset0.copy()

model1 = dset0.model1d.copy()
model1.VsvArr[:10]  =model1.VsvArr[:10]*(1+dm)
model1.VshArr[:10]  = model1.VshArr[:10]*(1+dm)
model1.VpfArr     = np.sqrt ( (model1.VphArr**2 - 2.*((model1.VshArr)**2) ) )
# model1.perturb(zmax=zmax, dm=dm, datatype=datatype)
# model1.perturb(zmax=zmax, dm=dm, datatype='vsh')
# dset0.perturb(inmodel=dset0.model1d.copy(), wavetype='ray', perlst=[10., 20., 30., 40])

dset1 = modesum.modesumASDF()
dset1.getmodel(inmodel=model1)
if run:
    dset1.run_disp(workingdir='/home/leon/code/CPSPy/kerneltest_ti', outfname='/home/leon/code/CPSPy/kerneltest_ti_1.asdf', nmodes=1,
               period=per)
    dset1.run_eigen(infname='/home/leon/code/CPSPy/kerneltest_ti_1.asdf', outfname='/home/leon/code/CPSPy/kerneltest_ti_1.asdf')

dset1.load('/home/leon/code/CPSPy/kerneltest_ti_1.asdf')

c0, cpre1, c = dset0.compare_disp(inmodel=model1, indbase=dset1, wavetype=wavetype, perlst=per)
dset0.compute_love_kernel()
c0, cpre2= dset0.perturb_love(inmodel=model1, wavetype=wavetype, perlst=per)
dset3.compute_love_kernel_montagner()
c0, cpre3= dset3.perturb_love(inmodel=model1, wavetype=wavetype, perlst=per)

fig, ax=plt.subplots()
plt.title(wavetype+' '+datatype+' perturb %g percent' %(dm*100)+' at top %g' %zmax +' km', fontsize=40)
# plt.plot(per, np.array(c0), 'x--', lw=3, ms=10, label='m0')
plt.plot(per, np.array(cpre1) - np.array(c0), 'o', ms=10, label='m0 velocity kernel predicted')
plt.plot(per, np.array(cpre2) - np.array(c0), '^', ms=10, label='m0 love parameter kernel predicted')
plt.plot(per, np.array(cpre3) - np.array(c0), 'v', ms=10, label='m0 love parameter kernel predicted, montagner')
plt.plot(per, np.array(c) - np.array(c0), 'ko', ms=10, label='m1')

# plt.plot(per, (np.array(cpre1) - np.array(c0))/(np.array(c) - np.array(c0)), 'o', ms=10, label='m0 kernel predicted')
# plt.plot(per, (np.array(cpre2) - np.array(c0))/(np.array(c) - np.array(c0)), 'o', ms=10, label='m0 love kernel predicted')

plt.legend(numpoints=1, fontsize=20, loc=0)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
plt.ylabel('C (km/s)', fontsize=30)
plt.xlabel('Period (sec)', fontsize=30)
plt.show()