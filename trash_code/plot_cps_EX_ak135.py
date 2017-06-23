import obspy
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from obspy.geodetics import kilometer2degrees


infname = '../cps_data_dir/syndata_dir_all_EX/Tgr_10.0.txt'
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
plt.ylabel('Vgr (km/s) ', fontsize=40);
plt.xlabel('Distance (km)', fontsize=40);
plt.title('Group Velocity', fontsize=40);
ax1.tick_params(axis='x', labelsize=30)
ax1.tick_params(axis='y', labelsize=30)
plt.ylim([2.94, 3.1])

fig2, ax2 = plt.subplots()
infname = '../cps_data_dir/syndata_dir_all_EX/Amp_10.0.txt'
# infname = '/media/lili/BCD29BBBD29B787A/cps_data_dir/syndata_dir_all_DD/Amp_10.0.txt'
inArr2=np.loadtxt(infname)
AmpArr=inArr2[:,2]*1000000.
DistArr=inArr2[:,3]
ind=np.argsort(DistArr)
DistArr=DistArr[ind]
AmpArr=AmpArr[ind]
mindist=DistArr.min()
indexmin=DistArr.argmin()
indexmin=np.where(DistArr==1000.)[0][0]
maxamp=AmpArr[indexmin]
plt.ylabel('Amplitude', fontsize=40);
plt.xlabel('Distance (km)', fontsize=40);

ax2.tick_params(axis='x', labelsize=30)
ax2.tick_params(axis='y', labelsize=30)

CampArr=AmpArr*np.sqrt(DistArr/1000. )  / maxamp
# plt.plot(DistArr,AmpArr,'g--o', markersize=10 );
plt.plot(DistArr,CampArr,'r--o', markersize=10 );
# plt.title('Amplitude ', fontsize=30);
plt.title('Corrected Normalized Amplitude', fontsize=40);
plt.ylim([0.9, 1.1])
plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(DistArr, DistArr/VgrArr);
print slope, intercept, r_value, p_value, std_err

