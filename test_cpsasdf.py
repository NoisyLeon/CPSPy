
import symdata

dbase1=symdata.cpsASDF('cpssynthetics_000_DD.h5')

# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
# dbase.PlotStreamsDistance()
dbase1.write2sac('../sac_DD_000')

dbase2=symdata.cpsASDF('cpssynthetics_001_DD.h5')

# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
# dbase.PlotStreamsDistance()
dbase2.write2sac('../sac_DD_001')

dbase3=symdata.cpsASDF('cpssynthetics_all_DD.h5')

# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
# dbase.PlotStreamsDistance()
dbase3.write2sac('../sac_DD_all')