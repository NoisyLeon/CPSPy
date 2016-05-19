import symdata

dbase=symdata.cpsASDF('cpssynthetics.h5')

dbase.Readsac('station.lst', './syndata_dir_all/sac_dir/', verbose=True)