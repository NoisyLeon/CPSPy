import subprocess, os

# os.makedirs('/home/leon/code/CPSPy/test_subprocess')
os.chdir('/home/leon/code/CPSPy/test_subprocess')
# subprocess.call('cp ../ak135.mod .')
subprocess.call(['sprep96', '-DT', '0.1', '-NPTS', '1024', '-M', 'ak135.mod', '-HS', '3.0', '-HR', '4.4', '-NMOD', '10', '-R', '-L'])