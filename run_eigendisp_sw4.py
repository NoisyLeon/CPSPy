import eigendisp

# dbase=eigendisp.eigendispASDF('/lustre/janus_scratch/life9360/cps_working_dir/eigendisp_homobasin.h5')
# dbase.sw4Vprofile('/lustre/janus_scratch/life9360/cps_working_dir/cpsinput_3km_2000.txt')
# dbase=eigendisp.eigendispASDF('/lustre/janus_scratch/life9360/cps_working_dir/eigendisp_ringbasin.h5')
# dbase.sw4Vprofile('/lustre/janus_scratch/life9360/cps_working_dir/cpsin.txt')
# dbase.getdisp(workingdir='/lustre/janus_scratch/life9360/cps_working_dir')


dbase=eigendisp.eigendispASDF('./eigendisp_ringbasin.h5')
dbase.sw4Vprofile('./cpsin.txt')
dbase.getdisp(workingdir='./cps_working_dir')