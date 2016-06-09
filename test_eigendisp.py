import eigendisp

dbase=eigendisp.eigendispASDF('./eigendisp.h5')
dbase.sw4Vprofile('../SW4Py/cpsinput_002.txt')
dbase.getdisp()