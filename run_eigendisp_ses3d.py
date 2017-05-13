import eigendisp


# minlat=24.
# maxlat=50.
# minlon=-125.0
# maxlon=-66.
# dbase = eigendisp.eigendispHDF5('../ses3d_US_model.eigenh5')
# # dbase.readh5model('../USmodel.h5', minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlon=0.5, dlat=0.5,
# #             groupname= 'US_model_410km_ses3d', maxdepth = 210., vsmin=2.0)
# # dbase.getdisp(workingdir='./cps_working_ses3d')
# dbase.get_2D_map(outdir='../Pygm/US_data')

minlat=22.
maxlat=52.
minlon=85.
maxlon=133.
dbase = eigendisp.eigendispHDF5('../ses3d_EA_model.eigenh5')
# dbase.readh5model('../EAmodel.h5', minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlon=0.5, dlat=0.5,
#             groupname= 'EA_model_210_ses3d', maxdepth = 210., vsmin=2.0)
# dbase.getdisp(workingdir='./cps_working_ses3d')
dbase.get_2D_map(outdir='../Pygm/EA_data')