import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from obspy.geodetics import kilometer2degrees


infname = '../cps_data_dir/syndata_dir_000/Tgr_10.0.txt'
# infname = '/media/lili/BCD29BBBD29B787A/cps_data_dir/syndata_dir_all_DD/Tgr_10.0.txt'
# infname = './syndata_dir_000_DD/Tgr_10.0.txt'
inArr=np.loadtxt(infname)
fig1, ax1 = plt.subplots()
T=inArr[:,2]
DistArr=inArr[:,3]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
T=T[ind]
DeltaArr=kilometer2degrees(DistArr)
VgrArr=DistArr/T
ax1.plot(DistArr, VgrArr,'--o' , markersize=10);
plt.ylabel('Vgr (km/s) ', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);
plt.title('Group Velocity', fontsize=30);
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)


fig2, ax2 = plt.subplots()
infname = '../cps_data_dir/syndata_dir_000/Amp_10.0.txt'
# infname = '/media/lili/BCD29BBBD29B787A/cps_data_dir/syndata_dir_all_DD/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]*1000000.
DistArr=inArr2[:,3]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
AmpArr=AmpArr[ind]
mindist=DistArr.min()
indexmin=DistArr.argmin()
maxamp=AmpArr[indexmin]
plt.ylabel('Amplitude', fontsize=30);
plt.xlabel('Distance (km)', fontsize=30);

ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)

CampArr=AmpArr*np.sqrt(DistArr/mindist )  / maxamp
# plt.plot(DistArr,AmpArr,'g--o', markersize=10 );
plt.plot(DistArr,CampArr,'r--o', markersize=10 );
# plt.title('Amplitude ', fontsize=30);
plt.title('Corrected Amplitude (normalized to the closest point)', fontsize=30);
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err

