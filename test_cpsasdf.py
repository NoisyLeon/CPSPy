import symdata

dbase=symdata.cpsASDF('cpssynthetics_000.h5')

# dbase.Readsac('station.lst', '../syndata_dir_000/sac_dir/', verbose=True)
dbase.PlotStreamsDistance()