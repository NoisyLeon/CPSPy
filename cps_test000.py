import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from obspy.geodetics import kilometer2degrees


infname = './syndata_dir_000/Tgr_10.0.txt'
infname = './syndata_dir_000_DD/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
plt.figure();
T=inArr[:,2]
DistArr=inArr[:,3]
DeltaArr=kilometer2degrees(DistArr)
VgrArr=DistArr/T
mindist=DistArr.min()
indexmin=DistArr.argmin()
plt.plot(DistArr, VgrArr,'o' );
# plt.plot(DistArr, (VgrArr-VgrArr[indexmin])/VgrArr[indexmin]*100.,'o' );
# plt.ylabel('Relative Difference in Vgr (%)');
plt.ylabel('Vgr(km/s)');
plt.xlabel('Distance(km)');
infname = './syndata_dir_000/Amp_10.0.txt'
infname = './syndata_dir_000_DD/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]
DistArr=inArr2[:,3]
plt.figure();
mindist=DistArr.min()
indexmin=DistArr.argmin()
maxamp=AmpArr[indexmin]
# plt.plot(DistArr, AmpArr*1e9,'o' );
plt.ylabel('Amplitude');
plt.xlabel('Distance(km)');
# 
# plt.figure();
# # plt.plot(DistArr, VgrArr, 'x');
# # plt.plot(DistArr, AmpArr*np.sqrt(np.sin(DeltaArr*np.pi/180.) ));
# CampArr=AmpArr*np.sqrt(np.sin(DeltaArr*np.pi/180.) )/np.sqrt(np.sin(DeltaArr[indexmin]*np.pi/180.) )
CampArr=AmpArr*np.sqrt(DistArr/mindist )  
# CampArr=AmpArr*DistArr/ DistArr[0] 
# plt.plot(DistArr, (CampArr-CampArr[indexmin])/CampArr[indexmin]*100.,'o' );
plt.plot(DistArr, CampArr,'o' );
# # plt.
# plt.ylabel('Relative Difference in Corrected Amp (%)');
# plt.xlabel('Distance(km)');
# # plt.axis([ DistArr.min(), DistArr.max(), CampArr.min(), CampArr.max()])
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err

