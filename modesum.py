
import numpy as np
import asdf
import vmodel
import cpsfile
import subprocess, os
import copy
import matplotlib.pyplot as plt


class modesumASDF(asdf.AsdfFile):
    """
    ASDF database for mode summation computation
    Note that the ASDF here is Advanced Scientific Data Format, NOT Adaptable Seismic Data Format !
    """
    def getmodel(self, inmodel=None):
        if isinstance(inmodel, vmodel.Model1d):
            self.model1d    = inmodel
        else:
            self.model1d    = vmodel.Model1d(modelindex=2)
            self.model1d.ak135()
            self.model1d.trim(zmax=400.)
            self.model1d=self.model1d.relayerize(h=2.)
        return
    
    def getdfile(self, distfile=None):
        if isinstance(distfile, cpsfile.DistFile()):
            self.distfile   = distfile
        elif os.path.isfile(distfile):
            self.distfile   = cpsfile.DistFile(distfname=distfile)
        return
    
    def load(self, infname): self.tree.update((asdf.AsdfFile.open(infname)).tree)
    
    
    def run_disp(self, workingdir, outfname=None, nmodes=1, hr=0., hs=0., love=True, rayleigh=True, dt=1., N2=14,
                freq=[], period=[], freqper='PER', run=True, xmin=1, xmax=100):
        ###
        # prepare for sprep96/tprep96
        ###
        if len(freq) != 0 and len(period) != 0:
            raise ValueError('Specify only one of the freq or period array!')
        if len(freq) != 0: outdata  = np.asarray(freq)
        else: outdata       = np.asarray(period)
        input_parameters    = {'inparam':{'workingdir': workingdir, 'nmodes': nmodes, 'hr': hr, 'hs': hs, 'love': love,
                'rayleigh': rayleigh, 'dt': dt, 'N2': N2, 'freq': freq, 'period': period, 'freqper': freqper, 'xmin': xmin, 'xmax': xmax,
                'redo': False}}
        try: self.model1d
        except:
            print 'WARNING: No input model! Use ak135 as default!'
            self.getmodel()
        if not os.path.isdir(workingdir): os.makedirs(workingdir)
        if self.model1d.modeltype == 'ISOTROPIC': command = ['sprep96', '-M', 'out.mod']
        elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC': command = ['tprep96', '-M', 'out.mod']
        self.model1d.write(outfname=workingdir+'/out.mod')
        try:
            self.distfile.write(distfname=workingdir+'/dfile')
            command.append('-d'); command.append('dfile')
        except: 
            command.append('-DT'); command.append('%g' %dt)
            npts = 2**N2
            command.append('-NPTS'); command.append('%g' %npts)
        command.append('-NMOD'); command.append('%g' %nmodes)
        command.append('-HS'); command.append('%g' %hs)
        command.append('-HR'); command.append('%g' %hr)
        if love: command.append('-L')
        if rayleigh: command.append('-R')
        if len(freq) == 1: command.append('-FREQ'); command.append('%g' %freq[0])
        elif len(freq) > 1:
            command.append('-FARR'); command.append('freq_lst')
            np.savetxt(workingdir+'/freq_lst', freq, fmt='%g')
        if len(period) == 1: command.append('-PER'); command.append('%g' %period[0])
        elif len(period) > 1:
            command.append('-PARR'); command.append('per_lst')
            np.savetxt(workingdir+'/per_lst', period, fmt='%g')
        command_prep = copy.deepcopy(command)
        ###
        # prepare for sdisp96/tdisp96
        ###
        if self.model1d.modeltype == 'ISOTROPIC': command_disp = ['sdisp96']
        elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC': command_disp = ['tdisp96']
        command.append('\n'+command_disp[0])
        ###
        # prepare for sdpsrf96/tdpsrf96
        ###
        if self.model1d.modeltype == 'ISOTROPIC': command_dpsrf=['sdpsrf96', '-U', '-C', '-G']
        elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC': command_dpsrf=['tdpsrf96', '-U', '-C', '-G']
        command_dpsrf.append('-'+freqper)
        command_dpsrf.append('-XMIN'); command_dpsrf.append('%g' %xmin)
        command_dpsrf.append('-XMAX'); command_dpsrf.append('%g' %xmax)
        command_dpsrf.append('-TXT')
        ###
        # saving input parameters
        ###
        self.tree.update(input_parameters)
        disptree = {'ray':{}, 'love': {}}
        if run:
            os.chdir(workingdir)
            subprocess.call(command_prep)
            subprocess.call(command_disp)
            if love:
                command_dpsrf_love = copy.deepcopy(command_dpsrf)
                command_dpsrf_love.append('-L')
                subprocess.call(command_dpsrf_love)
                if self.model1d.modeltype == 'ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/SDISPL.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/TDISPL.TXT')
                disptree.update(dispfile.get_tree(intree=disptree))
            if rayleigh:
                command_dpsrf_ray = copy.deepcopy(command_dpsrf)
                command_dpsrf_ray.append('-R')
                subprocess.call(command_dpsrf_ray)
                if self.model1d.modeltype == 'ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/SDISPR.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/TDISPR.TXT')
                disptree.update(dispfile.get_tree(intree=disptree))
            # subprocess.call(command_dpsrf)
            ###
            # Save dispersion curves
            ###
            self.tree.update({'disp': disptree})
            if outfname!=None: self.write_to(outfname)
            return True
        else:
            if love:
                command_dpsrf_love = copy.deepcopy(command_dpsrf)
                command_dpsrf_love.append('-L')
                command.append('\n'+command_dpsrf_love[0])
                command += command_dpsrf_love[1:]
            if rayleigh:
                command_dpsrf_ray = copy.deepcopy(command_dpsrf)
                command_dpsrf_ray.append('-R')
                command.append('\n'+command_dpsrf_ray[0])
                command += command_dpsrf_ray[1:]
            outstr = ' '.join(command)
            print outstr
            return outstr
        return
    
    def run_eigen(self, workingdir=None, infname=None, outfname=None, love=True, rayleigh=True, runtype='ALL', noq=True, hr=None, hs=None, freqper='PER',
                xmin=1, xmax=100, verbose=True, redo=False, run=True):
        if runtype != 'ALL'  and runtype != 'SYN' and runtype != 'DER' and runtype != 'DE' and runtype != 'DR' and runtype != 'DH'\
                and runtype != 'DA' and runtype != 'DB':
            raise ValueError('Unrecognized runtype (should be ALL, SYN, DER, DE, DR, DH, DA, DB)!')
        if infname != None: self.load(infname=infname)
        try:
            inparam    = copy.deepcopy(self.tree['inparam'])
            if workingdir == None: workingdir = inparam['workingdir']
            if hr==None: hr     = inparam['hr']
            if hs==None: hs     = inparam['hs']
            if not redo: redo   = inparam['redo']
        except:
            inparam     = {}
            print 'WARNING: No input parameters for run_disp'
        if workingdir == None: raise ValueError('Working directory is not specified!')
        if hr==None: hr = 0.
        if hs==None: hs = 0.
        ###
        # Update input parameters
        ###
        inparam.update({'workingdir': workingdir, 'hr': hr, 'hs': hs, 'love': love, 'rayleigh': rayleigh, 'freqper': freqper, 'xmin': xmin, 'xmax': xmax,
                'redo': False})
        command_love        = []
        command_ray         = []
        command_dpegn       = []
        command_dpegn_love  = []
        command_dpegn_ray   = []
        command_dpder       = []
        command_dpder_love  = []
        command_dpder_ray   = []
        command             = []
        ###############
        # Love command
        ###############
        if love:
            if self.model1d.modeltype == 'ISOTROPIC':
                command_love = ['slegn96']
                if redo:
                    if not os.path.isfile(workingdir+'/tsdisp96.lov'): raise ValueError('No tsdisp96.lov file!')
                    command_love.append('-T')
                else:
                    if not os.path.isfile(workingdir+'/sdisp96.lov'): raise ValueError('No sdisp96.lov file!')
            elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                command_love = ['tlegn96']
                if redo:
                    if not os.path.isfile(workingdir+'/ttdisp96.lov'): raise ValueError('No ttdisp96.lov file!')
                    command_love.append('-T')
                else:
                    if not os.path.isfile(workingdir+'/tdisp96.lov'): raise ValueError('No tdisp96.lov file!')
            command_love.append('-HS'); command_love.append('%g' %hs)
            command_love.append('-HR'); command_love.append('%g' %hr)
            if noq: command_love.append('-NOQ')
        command += command_love
        ###############
        # Rayleigh command
        ###############
        if rayleigh:
            if self.model1d.modeltype == 'ISOTROPIC':
                command_ray = ['sregn96']
                if redo:
                    if not os.path.isfile(workingdir+'/tsdisp96.ray'): raise ValueError('No tsdisp96.ray file!')
                    command_ray.append('-T')
                else:
                    if not os.path.isfile(workingdir+'/sdisp96.ray'): raise ValueError('No sdisp96.ray file!')
            elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                command_ray = ['tregn96']
                if redo:
                    if not os.path.isfile(workingdir+'/ttdisp96.ray'): raise ValueError('No ttdisp96.ray file!')
                    command_ray.append('-T')
                else:
                    if not os.path.isfile(workingdir+'/tdisp96.ray'): raise ValueError('No tdisp96.ray file!')
            command_ray.append('-HS'); command_ray.append('%g' %hs)
            command_ray.append('-HR'); command_ray.append('%g' %hr)
            if noq: command_ray.append('-NOQ')
        if love and rayleigh:
            command.append('\n'+command_ray[0])
            command += command_ray[1:]
        else: command += command_ray
        ###############
        # convert dispersion data to TXT 
        ###############
        if runtype == 'ALL' or runtype == 'SYN':
            if self.model1d.modeltype == 'ISOTROPIC': command_dpegn = ['sdpegn96', '-TXT', '-C', '-U', '-G']
            elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC': command_dpegn = ['tdpegn96', '-TXT', '-C', '-U', '-G']
            command_dpegn.append('-'+freqper)
            command_dpegn.append('-XMIN'); command_dpegn.append('%g' %xmin)
            command_dpegn.append('-XMAX'); command_dpegn.append('%g' %xmax)
            if love:
                command_dpegn_love  = copy.deepcopy(command_dpegn)
                command_dpegn_love.append('-L')
            if rayleigh:
                command_dpegn_ray   = copy.deepcopy(command_dpegn)
                command_dpegn_ray.append('-R')
        ###############
        # convert eigenfunction data to TXT 
        ###############
        if runtype != 'SYN':
            try:
                if len(self.tree['inparam']['period'])==0:
                    raise ValueError ('period array needs to be specified for eigenfunction extraction!')
            except:
                raise ValueError ('period array needs to be specified for eigenfunction extraction!')
            if self.model1d.modeltype == 'ISOTROPIC': command_dpder = ['sdpder96', '-TXT']
            elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC': command_dpder = ['tdpder96', '-TXT']
            if love:
                command_dpder_love  = copy.deepcopy(command_dpder)
                command_dpder_love.append('-L')
            if rayleigh:
                command_dpder_ray   = copy.deepcopy(command_dpder)
                command_dpder_ray.append('-R')
        ###############
        # run command
        ###############
        try:
            disptree    = copy.deepcopy(self.tree['disp'])
        except:
            self.tree.update( {'disp':{'ray':{}, 'love': {}}} )
            disptree    = {'ray':{}, 'love': {}}
        dertree={'model': {}, 'egn': {'ray':{}, 'love': {} } } 
        if run:
            os.chdir(workingdir)
            if runtype == 'ALL':
                subprocess.call(command_love)
                subprocess.call(command_ray)
                command_love.append('-DER')
                command_ray.append('-DER')
                subprocess.call(command_love)
                subprocess.call(command_ray)
            elif runtype == 'DER' or runtype == 'DH' or runtype == 'DA' or runtype == 'DB' or runtype == 'DR':
                command_love.append('-'+runtype)
                command_ray.append('-'+runtype)
                subprocess.call(command_love)
                subprocess.call(command_ray)
            else:
                subprocess.call(command_love)
                subprocess.call(command_ray)
            ###
            # Read/Upate dispersion data
            ###
            if len(command_dpegn_love) !=0:
                subprocess.call(command_dpegn_love)
                if self.model1d.modeltype == 'ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/SLEGN.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/TLEGN.TXT')
                disptree.update( dispfile.get_tree(intree=disptree) )
            if len(command_dpegn_ray) !=0:
                subprocess.call(command_dpegn_ray)
                if self.model1d.modeltype == 'ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/SREGN.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    dispfile = cpsfile.DispFile(workingdir+'/TREGN.TXT')
                disptree.update( dispfile.get_tree(intree=disptree) )
            self.tree['disp'].update(disptree)
            ###
            # Read eigenfunction data
            ###
            if len(command_dpder_love) !=0:
                subprocess.call(command_dpder_love)
                if self.model1d.modeltype == 'ISOTROPIC':
                    derfile     = cpsfile.derFile(workingdir+'/SLDER.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    derfile     = cpsfile.derFile(workingdir+'/TLDER.TXT')
                dertree = derfile.get_tree(intree=dertree)
            if len(command_dpder_ray) !=0:
                subprocess.call(command_dpder_ray)
                if self.model1d.modeltype == 'ISOTROPIC':
                    derfile     = cpsfile.derFile(workingdir+'/SRDER.TXT')
                elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
                    derfile     = cpsfile.derFile(workingdir+'/TRDER.TXT')
                dertree = derfile.get_tree(intree=dertree)
            self.tree.update({'der': dertree})
            ###
            # Update input parameters
            ###
            self.tree.update({'inparam': inparam})
            if outfname!=None: self.write_to(outfname)
            return True
        if runtype == 'ALL':
            command_love.append('-DER')
            command_ray.append('-DER')
            command.append('\n'+command_ray[0])
            command += command_ray[1:]
            command.append('\n'+command_love[0])
            command += command_love[1:]
        elif runtype == '-DER' or runtype == '-DH' or runtype == '-DA' or runtype == '-DB' or runtype == '-DR':
            command_love.append('-'+runtype)
            command_ray.append('-'+runtype)
            command = []
            command += command_love
            command.append('\n'+command_ray[0])
            command += command_ray[1:]
        if len(command_dpegn_love) !=0:
            command.append('\n'+command_dpegn_love[0])
            command += command_dpegn_love[1:]
        if len(command_dpegn_ray) !=0:
            command.append('\n'+command_dpegn_ray[0])
            command += command_dpegn_ray[1:]
        if len(command_dpder_love) !=0:
            command.append('\n'+command_dpder_love[0])
            command += command_dpder_love[1:]
        if len(command_dpder_ray) !=0:
            command.append('\n'+command_dpder_ray[0])
            command += command_dpder_ray[1:]
        outstr = ' '.join(command)
        print outstr
        return outstr
    
    def plot_eigen(self, wavetype, period, dtype, mode=0, zmax=9999.):
        dataArr = self.tree['der']['egn'][wavetype][mode][period][dtype]
        HArr    = self.tree['der']['model']['H']
        zArr    = np.cumsum(HArr) - HArr
        dataArr = dataArr[zArr<zmax]
        zArr    = zArr[zArr<zmax]
        fig, ax=plt.subplots()
        plt.title('')
        plt.plot(dataArr, zArr, 'k-', lw=3)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        plt.xlabel(dtype.upper(), fontsize=30)
        plt.ylabel('Depth (km)', fontsize=30)
        plt.gca().invert_yaxis()
        plt.show()
            
        
    
    
        
        