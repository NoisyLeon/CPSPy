
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
    def getmodel(self, inmodel=None, modelindex=1, h=2., zmax=400.):
        """
        Get velocity model
        =====================================================================
        Input parameters:
        inmodel     - input model
        modelindex  - model index (1 - ISOTROPIC, 2 - TRANSVERSE ISOTROPIC)
        h           - layer thickness for re-layerize
        zmax        - maximum depth fro trim 
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
    
    def getdfile(self, distfile=None):
        """
        Get distance file
        """
        if isinstance(distfile, cpsfile.DistFile()):
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
    
    def del_dir(self, workingdir=None):
        try:
            inparam    = copy.deepcopy(self.tree['inparam'])
            if workingdir == None: workingdir = inparam['workingdir']
        except:
            pass
        if workingdir == None: raise ValueError('Working directory is not specified!')
        shutil.rmtrees(workingdir)
        return
    
    def plot_eigen(self, wavetype, period, dtype, style=None, mode=0, zmax=9999., newfig=True, showfig=True):
        ###
        # Check input
        ###
        wavetype = wavetype.lower(); dtype   = dtype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        if wavetype != 'ray' and wavetype != 'love':
            raise ValueError('wavetype can only be ray or love')
        if dtype !='ur' and dtype !='tr' and dtype !='uz' and dtype !='tz' and dtype !='ut' and dtype !='tt' and dtype !='dcdh' and dtype !='dcda' \
            and dtype !='dcdb' and dtype !='dcdr' and dtype !='dcdav' and dtype !='dcdah' and dtype !='dcdbv' and dtype !='dcdbh' and dtype !='dcdn':
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
    
    def perturb(self, inmodel, wavetype, mode=0, perlst=[10.]):
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
    
    def compare_disp(self, inmodel, indbase, wavetype, mode=0, perlst=[10.]):
        wavetype = wavetype.lower()
        if wavetype == 'rayleigh': wavetype = 'ray'
        c0Arr, cpreArr = self.perturb(inmodel=inmodel, wavetype=wavetype, mode=mode, perlst=perlst)
        cArr = []
        for period in perlst:
            c       = indbase.tree['der']['egn'][wavetype][mode][period]['C']
            cArr.append(c)
        return c0Arr, cpreArr, cArr
        
                
    
    
        
        