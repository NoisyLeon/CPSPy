import numpy as np
import pyasdf
import cpsfile
import vmodel
import os
import warnings
import subprocess 

class eigendispASDF(pyasdf.ASDFDataSet):
    
    def sw4Vprofile(self, infname):
        self.Vprofile=vmodel.vprofile(infname=infname)
        self.Vprofile.HArr=self.Vprofile.HArr/1000.
        index = self.Vprofile.dtypeArr==0
        self.Vprofile.vsArr[index]=self.Vprofile.vsArr[index]/1000.
        self.Vprofile.vpArr[index]=self.Vprofile.vpArr[index]/1000.
        self.Vprofile.rhoArr[index]=self.Vprofile.rhoArr[index]/1000.
        self.Vprofile.z0Arr[index]=self.Vprofile.z0Arr[index]/1000.
        return
    
    def cpsmodel1d(self):
        self.model1d=vmodel.Model1d()
        self.model1d.ak135()
        return
    
    def getdisp(self, dt=0.1, N2=12 , nmodes=1, mode=0, hr=0., hs=0., workingdir='./cpspy_working_dir', deletemod=True, verbose=True):
        """ Get theoretical dispersion curves for all vertical profiles and save them to ASDF file
        ==========================================================================================
        Input Parameters:
        dt                 - time interval
        N2               - npts = 2**N2
        nmodes        - number of modes
        mode           - mode index (0 for fundamental mode, 1 for 1st overtone ...)
        hr                 - receiver depth
        hs                - source depth
        workingdir  - working directory for Computer Program for Seismology
        deletemod   - delete model files or not
        ==========================================================================================
        """
        npts=2**N2
        FNULL = open(os.devnull, 'w')
        if mode+1 > nmodes:
            nmodes=mode+1
            warnings.warn('Mode required is not available, re-assign nmodes = mode+1 !', UserWarning, stacklevel=1)
        if not os.path.isdir(workingdir):
            os.makedirs(workingdir)
        ak135mod=workingdir+'/ak135.mod'
        tempCPS=workingdir+'/run_dispersion.sh'
        try:
            self.model1d.write(ak135mod)
        except AttributeError:
            self.cpsmodel1d()
            self.model1d.write(ak135mod)
        with open(tempCPS, 'wb') as f:
            f.writelines('sprep96 -DT %f -NPTS %d -HR %f -HS %f -M %s -NMOD %d -R \n'
                        %(dt, npts, hr, hs, ak135mod, nmodes) )
            f.writelines('sdisp96 \nsregn96 -NOQ \nsdpegn96 -R -U -ASC -TXT -PER \n')
        if verbose == False:
            subprocess.call(['bash', tempCPS], stdout=FNULL, stderr=subprocess.STDOUT)
        else:
            subprocess.call(['bash', tempCPS])
        os.remove('sdisp96.dat')
        os.remove('sdisp96.ray')
        os.remove('sregn96.egn')
        os.remove('SREGN.ASC')
        os.remove('SREGNU.PLT')
        dispfile=cpsfile.DispFile('SREGN.TXT')
        os.remove('SREGN.TXT')
        os.remove(tempCPS)
        if deletemod==True:
            os.remove(ak135mod)
        dispcurve=dispfile.DispLst[mode]
        dispcurve.InterpDisp()
        auxArr=np.append(dispcurve.period, dispcurve.Vph)
        auxArr=np.append(auxArr, dispcurve.Vgr)
        auxArr=auxArr.reshape(3, dispcurve.period.size )
        # vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
        parameters={'Rmax': 0., 'Rmin': 0., 'x': 0, 'y': 0, 'T': 0, 'Vph': 1, 'Vgr': 2}
        self.add_auxiliary_data(data=auxArr, data_type='Disp', path='VP000', parameters=parameters)
        Np=self.Vprofile.vsArr.size
        ############################################################################
        # Get theoretical dispersion curve for each vertical profile 
        ############################################################################
        for i in xrange(Np):
            model1d=self.model1d.copy()
            if self.Vprofile.dtypeArr[i]==0:
                model1d.addlayer(H=self.Vprofile.HArr[i], vs=self.Vprofile.vsArr[i], vp=self.Vprofile.vpArr[i],
                        rho=self.Vprofile.rhoArr[i], zmin=self.Vprofile.z0Arr[i])
                tempmod=workingdir+'/temp_%03d' %i +'.mod' 
                model1d.write(tempmod)
            else:
                z1=self.Vprofile.z0Arr[i]
                z2=self.Vprofile.HArr[i] + z1
                if self.Vprofile.vsArr[i] !=0.:
                    model1d.perturb(dm=self.Vprofile.vsArr[i], zmin=z1, zmax=z2, datatype='vs')
                if self.Vprofile.vpArr[i] !=0.:
                    model1d.perturb(dm=self.Vprofile.vpArr[i], zmin=z1, zmax=z2, datatype='vp')
                if self.Vprofile.rhoArr[i] !=0.:
                    model1d.perturb(dm=self.Vprofile.rhoArr[i], zmin=z1, zmax=z2, datatype='rho')
                tempmod=workingdir+'/temp_%03d' %i +'.mod' 
                model1d.write(tempmod)
            with open(tempCPS, 'wb') as f:
                f.writelines('sprep96 -DT %f -NPTS %d -HR %f -HS %f -M %s -NMOD %d -R \n'
                            %(dt, npts, hr, hs, tempmod, nmodes) )
                f.writelines('sdisp96 \nsregn96 -NOQ \nsdpegn96 -R -U -ASC -TXT -PER \n')
            if verbose == False:
                subprocess.call(['bash', tempCPS], stdout=FNULL, stderr=subprocess.STDOUT)
            else:
                subprocess.call(['bash', tempCPS])
            os.remove('sdisp96.dat')
            os.remove('sdisp96.ray')
            os.remove('sregn96.egn')
            os.remove('SREGN.ASC')
            os.remove('SREGNU.PLT')
            dispfile=cpsfile.DispFile('SREGN.TXT')
            os.remove('SREGN.TXT')
            os.remove(tempCPS)
            if deletemod==True:
                os.remove(tempmod)
            dispcurve=dispfile.DispLst[mode]
            dispcurve.InterpDisp()
            auxArr=np.append(dispcurve.period, dispcurve.Vph)
            auxArr=np.append(auxArr, dispcurve.Vgr)
            auxArr=auxArr.reshape(3, dispcurve.period.size )
            # vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
            parameters={'Rmax': self.Vprofile.RmaxArr[i], 'Rmin': self.Vprofile.RminArr[i], 'x': self.Vprofile.xArr[i],
                        'y': self.Vprofile.yArr[i], 'T': 0, 'Vph': 1, 'Vgr': 2}
            path='VP%03d' %i 
            self.add_auxiliary_data(data=auxArr, data_type='Disp', path=path, parameters=parameters)
        FNULL.close()
        return
        
        
    
            
        
            
                
                
                
        
        
        
        
        
    
    
    
