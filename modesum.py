
import numpy as np
import asdf
import vmodel
import cpsfile
import subprocess, os
import copy
import matplotlib.pyplot as plt
import shutil

class modesumASDF(asdf.AsdfFile):
    """
    ASDF database for mode summation computation
    Note that the ASDF here is Advanced Scientific Data Format, NOT Adaptable Seismic Data Format !
    """
    
    def __str__(self):
        outstr  = '==================== ASDF Database for Mode Summation Computation ====================\n'
        try:
            inparam = self.tree['inparam']
            tempstr = ''
            tempstr +='---------------------------------- Input Parameters ----------------------------------\n'
            tempstr +='working directory: ' + inparam['workingdir'] + ' \n'
            tempstr +='number of modes = ' + str(inparam['nmodes']) + ' \n'
            tempstr +='Love wave : ' + str(inparam['love']) + ' Rayleigh wave : ' + str(inparam['rayleigh']) + ' \n'
            if len(inparam['period']) !=0:
                tempstr +='period = ' +  str(inparam['period']) + ' \n'
            if len(inparam['freq']) !=0:
                tempstr +='frequency = ' +  str(inparam['freq']) + ' \n'
            tempstr +='hr = ' +  str(inparam['hr']) + ' hr = ' +  str(inparam['hs']) + ' \n'
            # tempstr +='--------------------------------------------------------------------------------------\n'
            outstr  += tempstr
        except:
            outstr  +='----------------------------- No input parameters in database -------------------------\n'
        try:
            disp = self.tree['disp']
            tempstr = ''
            tempstr +='----------------------------------- Dispersion Data ----------------------------------\n'
            tempstr +='Love wave modes: '+ str(disp['love'].keys()) + ' \n'
            if inparam['love'] and len(disp['love'].keys()) != inparam['nmodes']:
                tempstr +='In compatible with input parameters\n' 
            tempstr +='Love wave period : '+ str(np.sort(disp['love'][0]['T'])) + ' \n'
            if not np.allclose(np.sort(inparam['period']), np.sort(disp['love'][0]['T']) ) : 
                tempstr +='In compatible with input parameters\n'
            tempstr +='Rayleigh wave modes : '+ str(disp['ray'].keys())    + ' \n'
            if inparam['rayleigh'] and len(disp['ray'].keys()) != inparam['nmodes']:
                tempstr +='In compatible with input parameters\n'
            tempstr +='Rayleigh wave period : '+ str(np.sort(disp['ray'][0]['T'])) + ' \n'
            if not np.allclose(np.sort(inparam['period']), np.sort(disp['ray'][0]['T']) ) : 
                tempstr +='In compatible with input parameters\n'
            outstr  += tempstr
        except:
            outstr  +='----------------------------- No dispersion data in database --------------------------\n'
        try:
            der     = self.tree['der']
            tempstr = ''
            tempstr +='--------------------------------- Eigenfunction Data ---------------------------------\n'
            model   = der['model']
            if model['isotropic']: tempstr +='Isotropic model \n'
            else: tempstr +='TI model \n'
            tempstr +='Love wave modes : '+ str(der['egn']['love'].keys()) + ' \n'
            if inparam['love'] and len(der['egn']['love'].keys()) != inparam['nmodes']:
                tempstr +='In compatible with input parameters\n'
            tempstr +='Love wave period : '+ str(np.sort(der['egn']['love'][0].keys())) + ' \n'
            if not np.allclose(np.sort(inparam['period']), np.sort(der['egn']['love'][0].keys()) ) : 
                tempstr +='In compatible with input parameters\n'
            tempstr +='Rayleigh wave modes : '+ str(der['egn']['ray'].keys())    + ' \n'
            if inparam['rayleigh'] and len(der['egn']['ray'].keys()) != inparam['nmodes']:
                tempstr +='In compatible with input parameters\n'
            tempstr +='Rayleigh wave period : '+ str(np.sort(der['egn']['love'][0].keys())) + ' \n'
            if not np.allclose(np.sort(inparam['period']), np.sort(der['egn']['ray'][0].keys()) ) : 
                tempstr +='In compatible with input parameters\n'
            outstr  += tempstr
        except:
            outstr  +='----------------------------- No dispersion data in database --------------------------\n'
        outstr +='======================================================================================\n'
        return outstr
    
    def __repr__(self): return self.__str__()
    
    def getmodel(self, inmodel=None, modelindex=1, h=2., zmax=400.):
        """
        Get velocity model
        =====================================================================
        Input parameters:
        inmodel     - input model
        modelindex  - model index (1 - ISOTROPIC, 2 - TRANSVERSE ISOTROPIC)
        h           - layer thickness for re-layerize
        zmax        - maximum depth for trim 
        =====================================================================
        """
        if isinstance(inmodel, vmodel.Model1d):
            self.model1d    = inmodel
        else:
            self.model1d    = vmodel.Model1d(modelindex=modelindex)
            self.model1d.ak135()
            self.model1d.trim(zmax=zmax)
            self.model1d=self.model1d.relayerize(h=h)
        return
    
    def getdfile(self, distfile=None, dist0=500, dD=100, Nd=20):
        """
        Get distance file
        """
        if distfile==None:
            self.distfile = cpsfile.DistFile()
            self.distfile.addEqualDist(dist0=dist0, dD=dD, Nd=Nd)
        else:
            if isinstance(distfile, cpsfile.DistFile):
                self.distfile   = distfile
            elif os.path.isfile(distfile):
                self.distfile   = cpsfile.DistFile(distfname=distfile)
        return
    
    def load(self, infname):
        """
        Load ASDF file
        """
        self.tree.update((asdf.AsdfFile.open(infname)).tree)
    
    
    def run_disp(self, workingdir, outfname=None, nmodes=1, hr=0., hs=0., love=True, rayleigh=True, dt=1., N2=14,
                freq=[], period=[], freqper='PER', run=True, xmin=1, xmax=100):
        """
        Compute dispersion curves
        =====================================================================
        Input parameters:
        workingdir  - working directory
        outfname    - output file name
        nmodes      - number of modes (default - 1)
        hr, hs      - depth of receiver/source
        love        - compute Love wave dispersion or not
        rayleigh    - compute Rayleigh wave dispersion or not
        dt          - time interval
        N2          - NPTS = 2**N2
        freq, period- frequency/period list
        freqper     - output x axis
        run         - run the code or not
        xmin, xmax  - output xmin/xmax
        =====================================================================
        """
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
            cdir    = os.getcwd()
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
            ###
            # Save dispersion curves
            ###
            self.tree.update({'disp': disptree})
            if outfname!=None: self.write_to(outfname)
            os.chdir(cdir)
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
        """
        Compute eigenfunctions
        ================================================================================
        Input parameters:
        workingdir  - working directory 
        infname     - input file name 
        outfname    - output file name 
        love        - compute Love wave dispersion or not
        rayleigh    - compute Rayleigh wave dispersion or not
        runtype     - type of run (ALL, SYN, DER, DE, DR, DH, DA, DB)
        noq         - computation with attenuation or not
        hr, hs      - depth of receiver/source
        freqper     - output x axis
        xmin, xmax  - output xmin/xmax
        redo        - redo the run with scomb96/tcomb96 or not (not implemented yet)
        run         - run the code or not
        ================================================================================
        """
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
        if runtype != 'SYN':
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
            cdir    = os.getcwd()
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
            os.chdir(cdir)
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
    
    def run_pulse(self, distfile=None, dist0=500, dD=100, Nd=10, workingdir=None, infname=None, outfname=None, stftype='-i',
            duraparam=2, sourcetype='-ALL', outtype='-D', fund=True, verbose=False, run=True):
        """
        Compute synthetics seismograms given eigenfunction data
        ================================================================================
        Input parameters:
        distfile    - distance file, object or file path
        dist0       - distance for origin point
        dD          - distance interval
        Nd          - number of distance point 
        workingdir  - working directory 
        infname     - input file name 
        outfname    - output file name 
        stftype     - type of source time function (-i, -o, -t, -p)
        duraparam   - duration controle parameters for source pulse
        sourcetype  - type of source (-EX, -EQ, -ALL)
        outtype     - output type (-D, -V, -A), displacement, velocity, acceleration
        fund        - only include fundamental mode or not
        run         - run the code or not
        ================================================================================
        """
        if infname != None: self.load(infname=infname)
        inparam    = copy.deepcopy(self.tree['inparam'])
        if workingdir == None: workingdir = inparam['workingdir']
        if outtype!='-D' and outtype!='-V' and outtype!='-A':
            raise ValueError('Unrecognized type of output seismogram '+ outtype +' should be -D, -V or -A')
        if sourcetype!='-EX' and sourcetype!='-EQ' and sourcetype!='-ALL':
            raise ValueError('Unrecognized type of souce '+ sourcetype +' should be -EX, -EQ or -ALL')
        if stftype!='-i' and stftype!='-t' and stftype!='-p' and stftype!='-o':
            raise ValueError('Unrecognized type of souce '+ stftype +' should be -i, -t, -p or -o')
        ###
        # distance file
        ###
        self.getdfile(distfile=distfile, dist0=dist0, dD=dD, Nd=Nd)
        try: self.distfile
        except: raise ValueError('Distance file not specified!')
        self.distfile.write(distfname=workingdir+'/dfile')
        ###
        #
        ###
        if self.model1d.modeltype == 'ISOTROPIC':
            command = ['spulse96']
        elif self.model1d.modeltype == 'TRANSVERSE ISOTROPIC':
            command = ['tpulse96']
        command.append('-d'); command.append('dfile')
        command.append(sourcetype); command.append(outtype)
        if verbose: command.append('-v')
        command.append(stftype)
        if stftype == '-t' or stftype == '-p':
            command.append('-l'); command.append('%d' %duraparam)
        if fund==True: command.append('-FUND')
        ###
        # Update input parameters
        ###
        inparam.update({'workingdir': workingdir, 'stf': stftype, 'sourcetype': sourcetype, 'outtype': outtype, 'fund': fund})
        self.tree.update({'inparam': inparam})
        command2sac=['f96tosac', '-B', 'tempf96']
        if outfname!=None: self.write_to(outfname)
        if run:
            cdir    = os.getcwd()
            os.chdir(workingdir)
            with open('tempf96', 'w') as f:
                subprocess.call(command, stdout=f)
            subprocess.call(command2sac)
            os.remove(workingdir+'/tempf96')
            os.chdir(cdir)
        else:
            command.append('>'); command.append('tempf96')
            command.append('\n'+command2sac[0])
            command+=command2sac[1:]
            outstr = ' '.join(command)
            print outstr
            return outstr
    
    def write_disp(self, outfname, wavetype='ray', mode=0, dtype='ph'):
        """
        Write dispersion curve to txt file
        ================================================================================
        Input parameters:
        outfname    - output file name 
        wavetype    - type of wave (Love or Rayleigh)
        mode        - mode id (0: fundamental mode, 1, 2,... overtones)
        dtype       - data type
        ================================================================================
        """
        wavetype = wavetype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        if wavetype != 'ray' and wavetype != 'love':
            raise ValueError('wavetype can only be ray or love')
        try:
            disptree= self.tree['disp'][wavetype][mode]
            if dtype == 'ph' or dtype == 'C': dtype = 'ph'
            elif dtype == 'gr' or dtype == 'U': dtype = 'gr'
            else: raise ValueError('Unrecognized data type!')
            vel     = disptree['V'+dtype]
            periods = disptree['T']
        except:
            egntree = self.tree['der']['egn'][wavetype][mode]
            periods = np.sort(egntree.keys())
            vel     = []
            if dtype == 'ph' or dtype == 'C': dtype = 'C'
            elif dtype == 'gr' or dtype == 'U': dtype = 'U'
            else: raise ValueError('Unrecognized data type!')
            for per in periods:
                vel.append( egntree[per][dtype] )
            vel = np.asarray(vel)
        outArr  = np.append(periods, vel)
        outArr  = (outArr.reshape(2, periods.size)).T
        np.savetxt(outfname, outArr, fmt='%g')
        return
        
    
    def del_dir(self, workingdir=None):
        """
        Delete working directory
        """
        try:
            inparam    = copy.deepcopy(self.tree['inparam'])
            if workingdir == None: workingdir = inparam['workingdir']
        except:
            pass
        if workingdir == None: raise ValueError('Working directory is not specified!')
        shutil.rmtree(workingdir)
        return
    
    def plot_eigen(self, wavetype, period, dtype, style=None, mode=0, zmax=9999., newfig=True, showfig=True):
        """
        Plot eigenfunction/sensitivity kernels
        ================================================================================
        Input parameters:
        wavetype    - type of wave (Love or Rayleigh)
        dtype       - data type
                    ------------------------- eigenfunctions ---------------------------
                    ur, tr, uz, tz  - Rayleigh wave eigenfunctions
                    ut, tt          - Love wave eigenfunctions
                    ----------------------- sensitivity kernels ------------------------
                    dcdh            - layer thickness
                    dcda, dcdb      - P/S wave velocity
                    dcdr            - density
                    dcdav/dcdah     - PV/PH wave velocity
                    dcdbv/dcdbh     - SV/SH wave velocity
                    dcdn            - eta (eta = F/(A-2L); A, L, F are Love parameters)
        style       - line style for plot 
        mode        - mode id (0: fundamental mode, 1, 2,... overtones)
        zmax        - maximum depth for trim
        ================================================================================
        """
        ###
        # Check input
        ###
        wavetype = wavetype.lower(); #dtype   = dtype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        if wavetype != 'ray' and wavetype != 'love':
            raise ValueError('wavetype can only be ray or love')
        if dtype !='ur' and dtype !='tr' and dtype !='uz' and dtype !='tz' and dtype !='ut' and dtype !='tt' and dtype !='dcdh' and dtype !='dcda' \
            and dtype !='dcdb' and dtype !='dcdr' and dtype !='dcdav' and dtype !='dcdah' and dtype !='dcdbv' and dtype !='dcdbh' and dtype !='dcdn'\
            and dtype !='dcdA' and dtype !='dcdC' and dtype !='dcdL' and dtype !='dcdF' and dtype !='dcdN':
            raise ValueError('dtype should be one of Ur, Tr, Uz, Tz, Ut, Tt, DCDH, DCDA, DCDB, DCDR, DCDAV, DCDAH, DCDBV, DCDBH, DCDN')
        if self.tree['der']['model']['isotropic']:
            if dtype == 'dcdav' or dtype == 'dcdah' or dtype == 'dcdbv' or dtype == 'dcdbh' or dtype == 'dcdn':
                raise ValueError('For isotropic model, dtype should NOT be one of DCDAV, DCDAH, DCDBV, DCDBH, DCDN')
        else:
            if dtype == 'dcda' or dtype == 'dcdb':
                raise ValueError('For TI model, dtype should NOT be one of DCDA, DCDB')
        ###
        # Get data for plot
        ###
        dataArr = self.tree['der']['egn'][wavetype][mode][period][dtype]
        HArr    = self.tree['der']['model']['H']
        zArr    = np.cumsum(HArr) - HArr
        dataArr = dataArr[zArr<zmax]
        zArr    = zArr[zArr<zmax]
        ###
        # Plot data
        ###
        if newfig: fig, ax=plt.subplots()
        plt.title('')
        if style !=None:
            plt.plot(dataArr, zArr, style, lw=3, label='mode: '+str(mode)+', T = '+str(period))
        else:
            plt.plot(dataArr, zArr, lw=3, label='mode: '+str(mode)+', T = '+str(period))
        if newfig: 
            ax.tick_params(axis='x', labelsize=20)
            ax.tick_params(axis='y', labelsize=20)
            plt.gca().invert_yaxis()
            plt.xlabel(dtype.upper(), fontsize=30)
            plt.ylabel('Depth (km)', fontsize=30)
        if showfig:
            plt.legend(numpoints=1, fontsize=20, loc=0)
            plt.show()
            
    
    def compute_love_kernel(self, outfname=None, modelst=None, perlst=None):
        refmodel    = self.tree['der']['model']
        if refmodel['isotropic']:  raise ValueError('Sensitivity kernels for Love parameters can not be computed for isotropic model')
        if modelst == None: modelst = self.tree['der']['egn']['ray'].keys()
        if perlst == None: perlst = self.tree['der']['egn']['ray'][modelst[0]].keys()
        ###
        # model parameters
        ###
        A       = refmodel['A']
        L       = refmodel['L']
        C       = refmodel['C']
        F       = refmodel['F']
        N       = refmodel['N']
        eta     = F/(A-2.*L)
        rho     = refmodel['rho']
        for mode in modelst:
            for period in perlst:
                VphR    = self.tree['der']['egn']['ray'][mode][period]['C']
                VgrR    = self.tree['der']['egn']['ray'][mode][period]['U']
                kR       = 2.*np.pi/VphR/period
                VphL    = self.tree['der']['egn']['love'][mode][period]['C']
                VgrL    = self.tree['der']['egn']['love'][mode][period]['U']
                kL       = 2.*np.pi/VphL/period
                # eigenfuntions
                Ur      = self.tree['der']['egn']['ray'][mode][period]['ur']
                Uz      = self.tree['der']['egn']['ray'][mode][period]['uz']
                Tr      = self.tree['der']['egn']['ray'][mode][period]['tr']
                Tz      = self.tree['der']['egn']['ray'][mode][period]['tz']
                Ut      = self.tree['der']['egn']['love'][mode][period]['ut']
                Tt      = self.tree['der']['egn']['love'][mode][period]['tt']
                # derivative of eigenfunctions, $5.8 of R.Herrmann
                durdz   = 1./L*Tr - kR*Uz
                duzdz   = kR*F/C*Ur + Tz/C
                dutdz   = Tt/L
                # Rayleigh wave sensitivity kernel
                dcRdav  = self.tree['der']['egn']['ray'][mode][period]['dcdav']
                dcRdah  = self.tree['der']['egn']['ray'][mode][period]['dcdah']
                dcRdbv  = self.tree['der']['egn']['ray'][mode][period]['dcdbv']
                dcRdbh  = self.tree['der']['egn']['ray'][mode][period]['dcdbh']
                dcRdeta = self.tree['der']['egn']['ray'][mode][period]['dcdn']
                dcRdr   = self.tree['der']['egn']['ray'][mode][period]['dcdr']
                # Love wave sensitivity kernel
                dcLdbv  = self.tree['der']['egn']['love'][mode][period]['dcdbv']
                dcLdbh  = self.tree['der']['egn']['love'][mode][period]['dcdbh']
                dcLdr   = self.tree['der']['egn']['love'][mode][period]['dcdr']
                # get trees for update
                rayTree = self.tree['der']['egn']['ray'][mode][period]
                loveTree= self.tree['der']['egn']['love'][mode][period]
                # compute partial derivatives using chain rule
                # Rayleigh wave
                dcRdA   = -dcRdeta * F/((A-2.*L)**2)  + dcRdah *0.5 /np.sqrt(rho*A)
                dcRdC   = dcRdav *0.5 /np.sqrt(rho*C)
                dcRdF   = dcRdeta /(A-2.*L)
                dcRdL   = dcRdeta * 2.*F/((A-2.*L)**2)  + dcRdbv *0.5 /np.sqrt(rho*L)
                dcRdN   = dcRdbh *0.5 /np.sqrt(rho*N)
                # Love wave
                dcLdA   = np.array([])
                dcLdC   = np.array([])
                dcLdF   = np.array([])
                dcLdL   = dcLdbv *0.5 /np.sqrt(rho*L)
                dcLdN   = dcLdbh *0.5 /np.sqrt(rho*N)
                # Update trees
                addRtree= {'dcdA': dcRdA, 'dcdC': dcRdC, 'dcdF': dcRdF, 'dcdL': dcRdL, 'dcdN': dcRdN}
                addLtree= {'dcdA': dcLdA, 'dcdC': dcLdC, 'dcdF': dcLdF, 'dcdL': dcLdL, 'dcdN': dcLdN}
                rayTree.update(addRtree)
                loveTree.update(addLtree)
        if outfname!=None: self.write_to(outfname)
        return
    
    def compute_love_kernel_montagner(self, outfname=None, modelst=None, perlst=None):
        refmodel    = self.tree['der']['model']
        if refmodel['isotropic']:  raise ValueError('Sensitivity kernels for Love parameters can not be computed for isotropic model')
        if modelst == None: modelst = self.tree['der']['egn']['ray'].keys()
        if perlst == None: perlst = self.tree['der']['egn']['ray'][modelst[0]].keys()
        ###
        # model parameters
        ###
        A       = refmodel['A']
        L       = refmodel['L']
        C       = refmodel['C']
        F       = refmodel['F']
        N       = refmodel['N']
        eta     = F/(A-2.*L)
        rho     = refmodel['rho']
        for mode in modelst:
            for period in perlst:
                VphR    = self.tree['der']['egn']['ray'][mode][period]['C']
                VgrR    = self.tree['der']['egn']['ray'][mode][period]['U']
                kR       = 2.*np.pi/VphR/period
                VphL    = self.tree['der']['egn']['love'][mode][period]['C']
                VgrL    = self.tree['der']['egn']['love'][mode][period]['U']
                kL       = 2.*np.pi/VphL/period
                # eigenfuntions
                Ur      = self.tree['der']['egn']['ray'][mode][period]['ur']
                Uz      = self.tree['der']['egn']['ray'][mode][period]['uz']
                Tr      = self.tree['der']['egn']['ray'][mode][period]['tr']
                Tz      = self.tree['der']['egn']['ray'][mode][period]['tz']
                Ut      = self.tree['der']['egn']['love'][mode][period]['ut']
                Tt      = self.tree['der']['egn']['love'][mode][period]['tt']
                # derivative of eigenfunctions, $5.8 of R.Herrmann
                durdz   = 1./L*Tr - kR*Uz
                duzdz   = kR*F/C*Ur + Tz/C
                dutdz   = Tt/L
                R0      = np.sum(rho*(Ur**2+Uz**2))
                L0      = np.sum(rho*(Ut**2))
                # Rayleigh wave sensitivity kernel
                dcRdav  = self.tree['der']['egn']['ray'][mode][period]['dcdav']
                dcRdah  = self.tree['der']['egn']['ray'][mode][period]['dcdah']
                dcRdbv  = self.tree['der']['egn']['ray'][mode][period]['dcdbv']
                dcRdbh  = self.tree['der']['egn']['ray'][mode][period]['dcdbh']
                dcRdeta = self.tree['der']['egn']['ray'][mode][period]['dcdn']
                dcRdr   = self.tree['der']['egn']['ray'][mode][period]['dcdr']
                # Love wave sensitivity kernel
                dcLdbv  = self.tree['der']['egn']['love'][mode][period]['dcdbv']
                dcLdbh  = self.tree['der']['egn']['love'][mode][period]['dcdbh']
                dcLdr   = self.tree['der']['egn']['love'][mode][period]['dcdr']
                # get trees for update
                rayTree = self.tree['der']['egn']['ray'][mode][period]
                loveTree= self.tree['der']['egn']['love'][mode][period]
                # compute partial derivatives using chain rule
                # Rayleigh wave
                dcRdA   = (Uz**2)/2./VgrR/R0
                dcRdC   = ((durdz/kR)**2)/2./VgrR/R0
                dcRdF   = (2.*durdz*Uz/kR)/2./VgrR/R0
                dcRdL   = ((duzdz/kR-Ur)**2)/2./VgrR/R0
                dcRdN   = dcRdbh *0.5 /np.sqrt(rho*N)
                # Love wave
                dcLdA   = np.array([])
                dcLdC   = np.array([])
                dcLdF   = np.array([])
                dcLdL   = ((dutdz/kL)**2)/2./VgrL/L0
                dcLdN   = (Ut**2)/2./VgrL/L0
                # Update trees
                addRtree= {'dcdA': dcRdA, 'dcdC': dcRdC, 'dcdF': dcRdF, 'dcdL': dcRdL, 'dcdN': dcRdN}
                addLtree= {'dcdA': dcLdA, 'dcdC': dcLdC, 'dcdF': dcLdF, 'dcdL': dcLdL, 'dcdN': dcLdN}
                rayTree.update(addRtree)
                loveTree.update(addLtree)
        if outfname!=None: self.write_to(outfname)
        return
    
    def perturb(self, inmodel, wavetype, mode=0, perlst=[10.]):
        """
        Compute phase velocity dispersion curve given a perturbed model
        ================================================================================
        Input parameters:
        inmodel     - input perturbed model
        wavetype    - type of wave (Love or Rayleigh)
        mode        - mode id (0: fundamental mode, 1, 2,... overtones)
        perlst      - period list
        ================================================================================
        """
        wavetype = wavetype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        refmodel    = self.tree['der']['model']
        if (refmodel['isotropic'] and inmodel.modeltype != 'ISOTROPIC') or ((not refmodel['isotropic']) and inmodel.modeltype == 'ISOTROPIC'):
            raise ValueError('Model need to be the same type!')
        if not np.allclose(refmodel['H'], inmodel.HArr):
            raise ValueError('Model layers are not the same!')
        HArr    = inmodel.HArr
        c0Arr   = []
        cArr    = []
        if refmodel['isotropic']:
            dvs = inmodel.VsArr - refmodel['vs']
            dvp = inmodel.VpArr - refmodel['vp']
            drho= inmodel.rhoArr - refmodel['rho']
            for period in perlst:
                dcda    = self.tree['der']['egn'][wavetype][mode][period]['dcda']
                dcdb    = self.tree['der']['egn'][wavetype][mode][period]['dcdb']
                dcdr    = self.tree['der']['egn'][wavetype][mode][period]['dcdr']
                # # # dc      = np.sum(dcda*dvp*HArr) + np.sum(dcdb*dvs*HArr) + np.sum(dcdr*drho*HArr)
                if wavetype=='ray':
                    dc      = np.sum(dcda*dvp) + np.sum(dcdb*dvs) + np.sum(dcdr*drho)
                else:
                    dc      = np.sum(dcdb*dvs) + np.sum(dcdr*drho)
                c0      = self.tree['der']['egn'][wavetype][mode][period]['C']
                c       = c0+dc
                c0Arr.append(c0)
                cArr.append(c)
        else:
            rho     = refmodel['rho']
            vsv0    = np.sqrt(refmodel['L'] / rho)
            vsh0    = np.sqrt(refmodel['N'] / rho)
            vpv0    = np.sqrt(refmodel['C'] / rho)
            vph0    = np.sqrt(refmodel['A'] / rho)
            eta0    = refmodel['F'] / (refmodel['A'] - 2.* refmodel['L'])
            eta     = (inmodel.VpfArr**2)/(inmodel.VphArr**2 - 2.*(inmodel.VsvArr**2))
            
            drho    = inmodel.rhoArr - rho
            dvsv    = inmodel.VsvArr - vsv0
            dvsh    = inmodel.VshArr - vsh0
            dvpv    = inmodel.VpvArr - vpv0
            dvph    = inmodel.VphArr - vph0
            deta    = eta - eta0
            # print 'HERE:', drho.min(), dvsv.max()
            for period in perlst:
                dcdav   = self.tree['der']['egn'][wavetype][mode][period]['dcdav']
                dcdah   = self.tree['der']['egn'][wavetype][mode][period]['dcdah']
                dcdbv   = self.tree['der']['egn'][wavetype][mode][period]['dcdbv']
                dcdbh   = self.tree['der']['egn'][wavetype][mode][period]['dcdbh']
                dcdeta  = self.tree['der']['egn'][wavetype][mode][period]['dcdn']
                dcdr    = self.tree['der']['egn'][wavetype][mode][period]['dcdr']
                # # # dc      = np.sum(dcda*dvp*HArr) + np.sum(dcdb*dvs*HArr) + np.sum(dcdr*drho*HArr)
                if wavetype=='ray':
                    dc      = np.sum(dcdav*dvpv) + np.sum(dcdah*dvph) + np.sum(dcdbv*dvsv) + np.sum(dcdbh*dvsh) \
                        + np.sum(dcdr*drho) + np.sum(dcdeta*deta)
                else:
                    dc      = np.sum(dcdbv*dvsv) + np.sum(dcdbh*dvsh) + np.sum(dcdr*drho)
                c0      = self.tree['der']['egn'][wavetype][mode][period]['C']
                c       = c0+dc
                c0Arr.append(c0)
                cArr.append(c)
            
        return c0Arr, cArr
    
    def perturb_love(self, inmodel, wavetype, mode=0, perlst=[10.]):
        """
        Compute phase velocity dispersion curve given a perturbed model,
        using partial derivatives for Love parameters.
        ================================================================================
        Input parameters:
        inmodel     - input perturbed model
        wavetype    - type of wave (Love or Rayleigh)
        mode        - mode id (0: fundamental mode, 1, 2,... overtones)
        perlst      - period list
        ================================================================================
        """
        wavetype = wavetype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        refmodel    = self.tree['der']['model']
        if refmodel['isotropic']:
            raise ValueError('This function only works for TI model')
        if not np.allclose(refmodel['H'], inmodel.HArr):
            raise ValueError('Model layers are not the same!')
        HArr    = inmodel.HArr
        c0Arr   = []
        cArr    = []
        A       = refmodel['A']
        L       = refmodel['L']
        C       = refmodel['C']
        F       = refmodel['F']
        N       = refmodel['N']
        rho0    = refmodel['rho']
        
        rho     = inmodel.rhoArr
        A1      = (inmodel.VphArr**2)*rho
        L1      = (inmodel.VsvArr**2)*rho
        C1      = (inmodel.VpvArr**2)*rho
        F1      = (inmodel.VpfArr**2)*rho
        N1      = (inmodel.VshArr**2)*rho
        
        drho    = inmodel.rhoArr - rho0
        dA      = A1 - A
        dC      = C1 - C
        dF      = F1 - F
        dL      = L1 - L
        dN      = N1 - N
        for period in perlst:
            dcdA    = self.tree['der']['egn'][wavetype][mode][period]['dcdA']
            dcdC    = self.tree['der']['egn'][wavetype][mode][period]['dcdC']
            dcdF    = self.tree['der']['egn'][wavetype][mode][period]['dcdF']
            dcdL    = self.tree['der']['egn'][wavetype][mode][period]['dcdL']
            dcdN    = self.tree['der']['egn'][wavetype][mode][period]['dcdN']
            dcdr    = self.tree['der']['egn'][wavetype][mode][period]['dcdr']
            # # # dc      = np.sum(dcda*dvp*HArr) + np.sum(dcdb*dvs*HArr) + np.sum(dcdr*drho*HArr)
            if wavetype=='ray':
                dc  = np.sum(dcdA*dA) + np.sum(dcdC*dC) + np.sum(dcdF*dF) + np.sum(dcdL*dL) \
                        + np.sum(dcdr*drho) + np.sum(dcdN*dN)
            else:
                dc  = np.sum(dcdL*dL) + np.sum(dcdr*drho) + np.sum(dcdN*dN)
            c0      = self.tree['der']['egn'][wavetype][mode][period]['C']
            c       = c0+dc
            c0Arr.append(c0)
            cArr.append(c)
            
        return c0Arr, cArr
    
    
    def compare_disp(self, inmodel, indbase, wavetype, mode=0, perlst=[10.]):
        """
        Compare dispersion curve from sensitivity kernel based computation with input database
        =======================================================================================
        Input parameters:
        inmodel     - input perturbed model
        indbase     - input database
        wavetype    - type of wave (Love or Rayleigh)
        mode        - mode id (0: fundamental mode, 1, 2,... overtones)
        perlst      - period list
        =======================================================================================
        """
        wavetype = wavetype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        c0Arr, cpreArr = self.perturb(inmodel=inmodel, wavetype=wavetype, mode=mode, perlst=perlst)
        cArr = []
        for period in perlst:
            c       = indbase.tree['der']['egn'][wavetype][mode][period]['C']
            cArr.append(c)
        return c0Arr, cpreArr, cArr
    
    def test_derivative(self, zmax, period, mode=0):
        refmodel    = self.tree['der']['model']
        ###
        # Get data for plot
        ###
        HArr    = self.tree['der']['model']['H']
        A       = refmodel['A']
        L       = refmodel['L']
        C       = refmodel['C']
        F       = refmodel['F']
        eta     = refmodel['F'] / (refmodel['A'] - 2.* refmodel['L'])
        rho     = refmodel['rho']
        VphR    = self.tree['der']['egn']['ray'][mode][period]['C']
        VgrR    = self.tree['der']['egn']['ray'][mode][period]['U']
        kR       = 2.*np.pi/VphR/period
        VphL    = self.tree['der']['egn']['love'][mode][period]['C']
        kL       = 2.*np.pi/VphL/period
        Ur      = self.tree['der']['egn']['ray'][mode][period]['ur']
        Uz      = self.tree['der']['egn']['ray'][mode][period]['uz']
        Tr      = self.tree['der']['egn']['ray'][mode][period]['tr']
        Tz      = self.tree['der']['egn']['ray'][mode][period]['tz']
        Ut      = self.tree['der']['egn']['love'][mode][period]['ut']
        Tt      = self.tree['der']['egn']['love'][mode][period]['tt']
        
        diffur1 = (Ur[1:] - Ur[:-1]) /1./HArr[1:]
        diffur2 = 1./L*Tr - kR*Uz
        
        diffuz1 = (Uz[1:] - Uz[:-1]) /1./HArr[1:]
        diffuz2 = kR*F/C*Ur + Tz/C
        
        diffut1 = (Ut[1:] - Ut[:-1]) /1./HArr[1:]
        diffut2 = Tt/L
        
        R0      = np.sum(rho*(Ur**2+Uz**2))
        dcdA1   = (Uz**2)/R0
        dcdah1  = self.tree['der']['egn']['ray'][mode][period]['dcdah']
        
        dcdah2  = 1./VgrR/R0*( eta*rho*np.sqrt(A/rho)*(Ur**2- 2*eta/kR*Ur*diffuz2) )
        
        # dcdeta  = 
        
        # return diffuz1, diffuz2
        # return diffut1, diffut2
        # return diffur1, diffur2
        return dcdah1, dcdah2
    
    def benchmark_montagner_nataf(self, zmax, period, mode=0):
        refmodel    = self.tree['der']['model']
        ###
        # Get data for plot
        ###
        # model parameters
        HArr    = self.tree['der']['model']['H']
        zArr    = np.cumsum(HArr)
        A       = refmodel['A']
        L       = refmodel['L']
        C       = refmodel['C']
        F       = refmodel['F']
        N       = refmodel['N']
        eta     = F/(A-2.*L)
        rho     = refmodel['rho']
        VphR    = self.tree['der']['egn']['ray'][mode][period]['C']
        VgrR    = self.tree['der']['egn']['ray'][mode][period]['U']
        kR       = 2.*np.pi/VphR/period
        VphL    = self.tree['der']['egn']['love'][mode][period]['C']
        VgrL    = self.tree['der']['egn']['love'][mode][period]['U']
        kL       = 2.*np.pi/VphL/period
        # eigenfuntions
        Ur      = self.tree['der']['egn']['ray'][mode][period]['ur']
        Uz      = self.tree['der']['egn']['ray'][mode][period]['uz']
        Tr      = self.tree['der']['egn']['ray'][mode][period]['tr']
        Tz      = self.tree['der']['egn']['ray'][mode][period]['tz']
        Ut      = self.tree['der']['egn']['love'][mode][period]['ut']
        Tt      = self.tree['der']['egn']['love'][mode][period]['tt']
        # derivative of eigenfunctions, $5.8 of R.Herrmann
        durdz   = 1./L*Tr - kR*Uz
        duzdz   = kR*F/C*Ur + Tz/C
        dutdz   = Tt/L
        # sensitivity kernel
        dcRdav  = self.tree['der']['egn']['ray'][mode][period]['dcdav']
        dcRdah  = self.tree['der']['egn']['ray'][mode][period]['dcdah']
        dcRdbv  = self.tree['der']['egn']['ray'][mode][period]['dcdbv']
        dcRdbh  = self.tree['der']['egn']['ray'][mode][period]['dcdbh']
        dcRdeta = self.tree['der']['egn']['ray'][mode][period]['dcdn']
        dcRdr   = self.tree['der']['egn']['ray'][mode][period]['dcdr']
        R0      = np.sum(rho*(Ur**2+Uz**2))
        # 1 : Montagner & Nataf, 2: derived from kernels computed by CPS
        dcdA1   = VphR/VgrR*( 1./2./VphR * (Uz**2/R0)) # if replace Uz with Ur, it matches...
        # dcdA1   = 1./VgrR/R0*(0.5*eta*Ur**2- eta/kR*Ur*duzdz*(1-eta))
        # because: eta = F/(A-2L), ah = sqrt(A/rho)
        # thus: dc/dA = dc/deta * deta/dA + dc/dah * dah/dA
        #             = -dc/deta * F/(A-2L)**2 + dc/dah * 1/2/sqrt(rho*A)
        dcdA2   = -dcRdeta * F/((A-2.*L)**2)  + dcRdah *0.5 /np.sqrt(rho*A) 
        
        dcdA3   = 1./VgrR/R0*(0.5*eta*Ur**2- eta/kR*Ur*duzdz*(1-eta))
        plt.plot(dcdA1, zArr-HArr[0]/2, 'b-', lw=3, label='dcR/dA Montagner')
        plt.plot(dcdA3, zArr-2., 'g.-', lw=3, label='dcR/dA derived')
        plt.plot(dcdA2, zArr-HArr[0]/2, 'r--', lw=3, label='dcR/dA CPS')
        plt.gca().invert_yaxis()
        plt.xlabel('dC/dA', fontsize=30)
        plt.ylabel('Depth (km)', fontsize=30)
        plt.ylim([50, 0])
        plt.legend(fontsize=20, loc=0)
        plt.show()
        
    def get_disp_egn(self, wavetype, mode=0, perlst=[10.]):
        cArr = []
        for period in perlst:
            c       = self.tree['der']['egn'][wavetype][mode][period]['C']
            cArr.append(c)
        return np.array(cArr)
        
        
        
        
                
    
    
        
        