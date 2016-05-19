import symdata
# import noisepy as npy
import obspy 
# dbase=symdata.cpsASDF('cpssynthetics_all_EX.h5') -0.00663607176745
# dbase=symdata.cpsASDF('cpssynthetics_all_DD.h5')
# dbase=symdata.cpsASDF('cpssynthetics_all_EX.h5')
# dbase=symdata.cpsASDF('cpssynthetics_all.h5')
dbase=symdata.cpsASDF('cpssynthetics_all_DD.h5')
# dbase.Readsac('station.lst', './syndata_dir_all/sac_dir/', comptype=['ZDD'], verbose=True)

inftan=symdata.InputFtanParam()
inftan.ffact=5.
inftan.pmf=True
try:
    del dbase.auxiliary_data.DISPbasic1
    del dbase.auxiliary_data.DISPbasic2
    del dbase.auxiliary_data.DISPpmf1
    del dbase.auxiliary_data.DISPpmf2
except:
    print 'aftan data already deleleted'
dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/cps_working_dir', inftan=inftan, tb=-0.752628782077, basic2=True,
            pmf1=True, pmf2=True)
# dbase.aftan(inftan=inftan, tb=-0.729324905982, basic2=True,
#             pmf1=True, pmf2=True)
# # # del dbase.events
del dbase.auxiliary_data.DISPbasic1interp
dbase.InterpDisp(data_type='DISPbasic1')
# # dbase.InterpDisp(data_type='DISPpmf2')
# # ndbase=dbase.SelectData(outfname='sw4synthetics001.h5', stafile='station_001.lst')
# 
# # inv=dbase.Readsac('ak135_station.lst',datadir='/home/lili/sw4_working_dir/ak135_001', comptype='u')
# # dbase.AddEvent(x=1000, y=1000, z=0)
# 
# # del dbase.auxiliary_data.FieldDISPbasic1interp
# # dbase.GetField(outdir='./homo', fieldtype='amp', data_type='DISPpmf2')
dbase.GetField(outdir='./syndata_dir_all_DD', fieldtype='amp',  data_type='DISPbasic1')
dbase.GetField(outdir='./syndata_dir_all_DD', fieldtype='Vgr',  data_type='DISPbasic1')
# 
# dbase.aftanMP(outdir='/lustre/janus_scratch/life9360/sw4_working_dir/DISP', inftan=inftan, basic2=True,
#             pmf1=True, pmf2=True)