import cpsfile

dispfile=cpsfile.DispFile('./disp_dir/SREGN_vsvp_+0.1.TXT')
dispfile.write('./disp_dir/ak135_interp_vsvp_+0.1.disp', datatype='phase')
dispfile=cpsfile.DispFile('./disp_dir/SREGN_vsvp_+0.2.TXT')
dispfile.write('./disp_dir/ak135_interp_vsvp_+0.2.disp', datatype='phase')
dispfile=cpsfile.DispFile('./disp_dir/SREGN_vsvp_+0.5.TXT')
dispfile.write('./disp_dir/ak135_interp_vsvp_+0.5.disp', datatype='phase')
# dispfile=cpsfile.DispFile('./disp_dir/SREGN_vs_-0.1.TXT')
# dispfile.write('./disp_dir/ak135_interp_-.disp')
# dispfile=cpsfile.DispFile('./disp_dir/SREGN_ref.TXT')
# dispfile.write('./disp_dir/ak135_interp.disp', datatype='phase')

