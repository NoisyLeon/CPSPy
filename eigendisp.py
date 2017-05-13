import numpy as np
import pyasdf, h5py
import cpsfile
import vmodel
import os
import warnings
import subprocess 


def _write_txt(fname, outlon, outlat, outZ):
    outArr=np.append(outlon, outlat)
    outArr=np.append(outArr, outZ)
    outArr=outArr.reshape((3,outZ.size))
    outArr=outArr.T
    np.savetxt(fname, outArr, fmt='%g')
    return

class eigendispASDF(pyasdf.ASDFDataSet):
    """
    An object for generating theoretical dispersion curve for input from SW4
    """
    def sw4Vprofile(self, infname):
        self.Vprofile               = vmodel.vprofile(infname=infname)
        self.Vprofile.HArr          = self.Vprofile.HArr/1000.
        index                       = self.Vprofile.dtypeArr==0
        self.Vprofile.vsArr[index]  = self.Vprofile.vsArr[index]/1000.
        self.Vprofile.vpArr[index]  = self.Vprofile.vpArr[index]/1000.
        self.Vprofile.rhoArr[index] = self.Vprofile.rhoArr[index]/1000.
        self.Vprofile.z0Arr[index]  = self.Vprofile.z0Arr[index]/1000.
        return
    
    def cpsmodel1d(self):
        self.model1d=vmodel.Model1d()
        self.model1d.ak135()
        return
    
    def getdisp(self, workingdir, dt=0.1, N2=12 , nmodes=1, mode=0, hr=0., hs=0., deletemod=False, verbose=True, noq=False):
        """ Get theoretical dispersion curves for all vertical profiles and save them to ASDF file
        ==========================================================================================
        Input Parameters:
        dt          - time interval
        N2          - npts = 2**N2
        nmodes      - number of modes
        mode        - mode index (0 for fundamental mode, 1 for 1st overtone ...)
        hr          - receiver depth
        hs          - source depth
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
            if noq: f.writelines('sdisp96 \nsregn96 -NOQ \nsdpegn96 -R -U -ASC -TXT -PER \n')
            else: f.writelines('sdisp96 \nsregn96 \nsdpegn96 -R -U -ASC -TXT -PER \n')        
        if verbose == False: subprocess.call(['bash', tempCPS], stdout=FNULL, stderr=subprocess.STDOUT)
        else: subprocess.call(['bash', tempCPS])
        os.remove('sdisp96.dat')
        os.remove('sdisp96.ray')
        os.remove('sregn96.egn')
        os.remove('SREGN.ASC')
        os.remove('SREGNU.PLT')
        dispfile=cpsfile.DispFile('SREGN.TXT')
        os.remove('SREGN.TXT')
        os.remove(tempCPS)
        if deletemod: os.remove(ak135mod)
        dispcurve=dispfile.DispLst[mode]
        dispcurve.InterpDisp()
        auxArr=np.append(dispcurve.period, dispcurve.Vph)
        auxArr=np.append(auxArr, dispcurve.Vgr)
        if noq:
            auxArr=auxArr.reshape(3, dispcurve.period.size )
            # vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
            parameters={'Rmax': 0., 'Rmin': 0., 'x': 0, 'y': 0, 'T': 0, 'Vph': 1, 'Vgr': 2}
        else:
            auxArr=np.append(auxArr, dispcurve.gamma)
            auxArr=auxArr.reshape(4, dispcurve.period.size )
            parameters={'Rmax': 0., 'Rmin': 0., 'x': 0, 'y': 0, 'T': 0, 'Vph': 1, 'Vgr': 2, 'gamma': 3}
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
                if noq: f.writelines('sdisp96 \nsregn96 -NOQ \nsdpegn96 -R -U -ASC -TXT -PER \n')
                else: f.writelines('sdisp96 \nsregn96 \nsdpegn96 -R -U -ASC -TXT -PER \n')        
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
            # os.remove(tempCPS)
            if deletemod: os.remove(tempmod)
            dispcurve=dispfile.DispLst[mode]
            dispcurve.InterpDisp()
            auxArr=np.append(dispcurve.period, dispcurve.Vph)
            auxArr=np.append(auxArr, dispcurve.Vgr)
            if noq:
                auxArr=auxArr.reshape(3, dispcurve.period.size )
                # vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
                parameters={'Rmax': self.Vprofile.RmaxArr[i], 'Rmin': self.Vprofile.RminArr[i], 'x': self.Vprofile.xArr[i],
                            'y': self.Vprofile.yArr[i], 'T': 0, 'Vph': 1, 'Vgr': 2}
            else:
                auxArr=np.append(auxArr, dispcurve.gamma)
                auxArr=auxArr.reshape(4, dispcurve.period.size )
                parameters={'Rmax': self.Vprofile.RmaxArr[i], 'Rmin': self.Vprofile.RminArr[i], 'x': self.Vprofile.xArr[i],
                            'y': self.Vprofile.yArr[i], 'T': 0, 'Vph': 1, 'Vgr': 2, 'gamma': 3}
            path='VP%03d' %(i+1) 
            self.add_auxiliary_data(data=auxArr, data_type='Disp', path=path, parameters=parameters)
        FNULL.close()
        return
    
class eigendispHDF5(h5py.File):
    """
    An object for generating theoretical dispersion curve for input from SES3D model
    """
    def readh5model(self, infname, groupname, minlat=-999, maxlat=999, minlon=-999, maxlon=999, dlat=None, dlon=None, maxdepth=None, vsmin=None):
        """
        Read hdf5 model
        ================================================================================================
        Input parameters:
        infname         - input filename
        groupname       - group name
        minlon, maxlon  - defines study region, default is to read corresponding data from hdf5 file 
        minlat, maxlat  -
        maxdepth        - maximum depth to be truncated
        ================================================================================================
        """
        header={'H': 0, 'vs':1, 'vp':2, 'rho':3, 'Qs':4}
        if vsmin!=None:
            vpmin=0.9409+2.0947*vsmin-0.8206*vsmin**2+0.2683*vsmin**3-0.0251*vsmin**4
            rhomin=1.6612*vpmin-0.4721*vpmin**2+0.0671*vpmin**3-0.0043*vpmin**4+0.000106*vpmin**5
        MDataset = h5py.File(infname)
        # get latitude/longitude information
        if minlat < MDataset[groupname].attrs['minlat']: minlat = MDataset[groupname].attrs['minlat']
        if maxlat > MDataset[groupname].attrs['maxlat']: maxlat = MDataset[groupname].attrs['maxlat']
        if minlon < MDataset[groupname].attrs['minlon']: minlon = MDataset[groupname].attrs['minlon']
        if maxlon > MDataset[groupname].attrs['maxlon']: maxlon = MDataset[groupname].attrs['maxlon']
        if dlon==None or dlat ==None:
            dlon = MDataset[groupname].attrs['dlon']; dlat = MDataset[groupname].attrs['dlat']
        self.attrs.create(name = 'dlon', data=dlon, dtype='f')
        self.attrs.create(name = 'dlat', data=dlat, dtype='f')
        self.attrs.create(name = 'minlat', data=minlat, dtype='f')
        self.attrs.create(name = 'maxlat', data=maxlat, dtype='f')
        self.attrs.create(name = 'minlon', data=minlon, dtype='f')
        self.attrs.create(name = 'maxlon', data=maxlon, dtype='f')
        # determine number of subvolumes, max depth and whether to interpolate or not 
        if maxdepth == None:
            dz          = MDataset[groupname].attrs['dz']
            depth       = MDataset[groupname].attrs['depth']
            depthArr    = MDataset[groupname].attrs['depthArr']
        else:
            dz          = MDataset[groupname].attrs['dz']
            depth       = MDataset[groupname].attrs['depth']
            depthArr    = MDataset[groupname].attrs['depthArr']
            if maxdepth > depth[-1]: raise ValueError('maximum depth is too large!')
            depth = depth [ np.where(depth<maxdepth)[0] ]
            bdz = dz[depth.size]
            if depth.size == 0 and maxdepth % bdz != 0:
                maxdepth = int( maxdepth / bdz ) * bdz
                print 'actual max depth:', maxdepth, 'km'
            elif ( maxdepth - depth[-1]) % bdz !=0:
                maxdepth = int( ( maxdepth - depth[-1]) / bdz ) * bdz + depth[-1]
                print 'actual max depth:', maxdepth, 'km'
            depth = np.append(depth, maxdepth)
            dz = dz[:depth.size]
        ##############################
        # Get verfical profile information
        ##############################
        # generate node array
        i=0
        for depb, dzb in zip(depth, dz): 
            if i==0:
                blockdep    = depthArr[depthArr<depb]
                nodeArr     = blockdep - dzb/2.
            else:
                dep0        = depth[i-1]
                blockdep    = depthArr[(depthArr<depb)*(depthArr>dep0)]
                nodeArr     = np.append(nodeArr, blockdep - dzb/2.)
            if i == depth.size-1: nodeArr     = np.append(nodeArr, depb)
            i+=1
        ## Getting data
        if MDataset[groupname].attrs['isblock']:
            group = MDataset[groupname]
        else:
            try: group = MDataset[groupname+'_block']
            except:
                group = MDataset[groupname]
                warnings.warn('Input model is NOT block model! ', UserWarning, stacklevel=1)
        latarr = minlat + np.arange( (maxlat-minlat)/dlat + 1)*dlat
        lonarr = minlon + np.arange( (maxlon-minlon)/dlon + 1)*dlon
        Nlayer  = nodeArr.size-1
        Harr    = nodeArr[1:] - nodeArr[:-1]
        outgroup = self.create_group( name = '3Dmodel' )
        for hkey in header.keys(): outgroup.attrs.create(name = hkey, data=header[hkey], dtype='i')
        for lat in latarr:
            for lon in lonarr:
                name='%g_%g' %(lon, lat)
                inProf  = group[name][...]
                HProf   = inProf[:Nlayer, :]
                HProf[:, 0] = Harr
                if vsmin != None:
                    vsArr = HProf[:, header['vs']]
                    vsArr[vsArr<vsmin] = vsmin
                    HProf[:, header['vs']] = vsArr
                    
                    vpArr = HProf[:, header['vp']]
                    vpArr[vpArr<vpmin] = vpmin
                    HProf[:, header['vp']] = vpArr
                    
                    rhoArr = HProf[:, header['rho']]
                    rhoArr[rhoArr<rhomin] = rhomin
                    HProf[:, header['rho']] = rhoArr
                dset = outgroup.create_dataset( name=name, shape=HProf.shape, data=HProf)
                dset.attrs.create(name = 'lon', data=lon, dtype='f')
                dset.attrs.create(name = 'lat', data=lat, dtype='f')
        return
    
    def getdisp(self, workingdir, dt=0.1, N2=12 , nmodes=1, mode=0, hr=0., hs=0., deletemod=False, verbose=True, noq=True):
        """ Get theoretical dispersion curves for all vertical profiles and save them to ASDF file
        ==========================================================================================
        Input Parameters:
        dt          - time interval
        N2          - npts = 2**N2
        nmodes      - number of modes
        mode        - mode index (0 for fundamental mode, 1 for 1st overtone ...)
        hr          - receiver depth
        hs          - source depth
        workingdir  - working directory for Computer Program for Seismology
        deletemod   - delete model files or not
        ==========================================================================================
        """
        ingroup = self['3Dmodel']
        header={'H': 0, 'vs':1, 'vp':2, 'rho':3, 'Qs':4}
        npts=2**N2
        FNULL = open(os.devnull, 'w')
        if mode+1 > nmodes:
            nmodes=mode+1
            warnings.warn('Mode required is not available, re-assign nmodes = mode+1 !', UserWarning, stacklevel=1)
        if not os.path.isdir(workingdir): os.makedirs(workingdir)
        minlat = self.attrs['minlat']
        maxlat = self.attrs['maxlat']
        minlon = self.attrs['minlon']
        maxlon = self.attrs['maxlon']
        dlon = self.attrs['dlon']; dlat = self.attrs['dlat']
        latarr = minlat + np.arange( (maxlat-minlat)/dlat + 1)*dlat
        lonarr = minlon + np.arange( (maxlon-minlon)/dlon + 1)*dlon
        tempCPS=workingdir+'/run_dispersion.sh'
        ############################################################################
        # Get theoretical dispersion curve for each vertical profile 
        ############################################################################
        outgroup= self.create_group( name = 'phase_vel' )
        for lat in latarr:
            for lon in lonarr:
                name='%g_%g' %(lon, lat)
                VProf  = ingroup[name][...]
                model1d=vmodel.Model1d()
                QpArr=3/4.*VProf[:,header['Qs']]*(VProf[:,header['vp']]/VProf[:,header['vs']])**2
                model1d.getmodel(modelname=name, HArr=VProf[:,header['H']], VpArr=VProf[:,header['vp']], VsArr=VProf[:,header['vs']],
                        rhoArr=VProf[:,header['rho']], QpArr=QpArr, QsArr=VProf[:,header['Qs']])
                model1d.check_model()
                tempmod=workingdir+'/temp_'+ name +'.mod' 
                model1d.write(tempmod)
                with open(tempCPS, 'wb') as f:
                    f.writelines('sprep96 -DT %f -NPTS %d -HR %f -HS %f -M %s -NMOD %d -R \n'
                                %(dt, npts, hr, hs, tempmod, nmodes) )
                    if noq: f.writelines('sdisp96 \nsregn96 -NOQ \nsdpegn96 -R -U -ASC -TXT -PER \n')
                    else: f.writelines('sdisp96 \nsregn96 \nsdpegn96 -R -U -ASC -TXT -PER \n')        
                if verbose == False: subprocess.call(['bash', tempCPS], stdout=FNULL, stderr=subprocess.STDOUT)
                else: subprocess.call(['bash', tempCPS])
                os.remove('sdisp96.dat')
                os.remove('sdisp96.ray')
                os.remove('sregn96.egn')
                os.remove('SREGN.ASC')
                os.remove('SREGNU.PLT')
                dispfile=cpsfile.DispFile('SREGN.TXT')
                os.remove('SREGN.TXT')
                os.remove(tempCPS)
                if deletemod: os.remove(tempmod)
                dispcurve=dispfile.DispLst[mode]
                dispcurve.InterpDisp()
                auxArr=np.append(dispcurve.period, dispcurve.Vph)
                auxArr=np.append(auxArr, dispcurve.Vgr)
                if noq:
                    auxArr=auxArr.reshape(3, dispcurve.period.size )
                    # # vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
                    # parameters={'Rmax': self.Vprofile.RmaxArr[i], 'Rmin': self.Vprofile.RminArr[i], 'x': self.Vprofile.xArr[i],
                    #             'y': self.Vprofile.yArr[i], 'T': 0, 'Vph': 1, 'Vgr': 2}
                else:
                    auxArr=np.append(auxArr, dispcurve.gamma)
                    auxArr=auxArr.reshape(4, dispcurve.period.size )
                    # parameters={'Rmax': self.Vprofile.RmaxArr[i], 'Rmin': self.Vprofile.RminArr[i], 'x': self.Vprofile.xArr[i],
                    #             'y': self.Vprofile.yArr[i], 'T': 0, 'Vph': 1, 'Vgr': 2, 'gamma': 3}
                dset = outgroup.create_dataset( name=name, shape=auxArr.shape, data=auxArr)
                dset.attrs.create(name = 'lon', data=lon, dtype='f')
                dset.attrs.create(name = 'lat', data=lat, dtype='f')
        FNULL.close()
        return
    
    def get_2D_map(self, outdir, pers=np.arange(10., 105., 5.), outtype='ph'):
        minlat  = self.attrs['minlat']
        maxlat  = self.attrs['maxlat']
        minlon  = self.attrs['minlon']
        maxlon  = self.attrs['maxlon']
        dlon    = self.attrs['dlon']; dlat = self.attrs['dlat']
        latarr  = minlat + np.arange( (maxlat-minlat)/dlat + 1)*dlat
        lonarr  = minlon + np.arange( (maxlon-minlon)/dlon + 1)*dlon
        ingroup = self['phase_vel']
        for per in pers:
            outfname = outdir+'/'+outtype+'V_%g.lst' %per
            outlat=np.array([])
            outlon=np.array([])
            outVarr=np.array([])
            for lat in latarr:
                for lon in lonarr:
                    name='%g_%g' %(lon, lat)
                    indisp  = ingroup[name][...]
                    outlat  = np.append(outlat, lat)
                    if lon < 0: lon+=360
                    outlon  = np.append(outlon, lon)
                    inpers  = indisp[0,:]
                    inphV   = indisp[1,:]
                    ingrV   = indisp[2,:]
                    if outtype=='ph': outV = inphV[inpers==per]
                    elif outtype=='gr': outV = ingrV[inpers==per]
                    outVarr=np.append(outVarr, outV)
            _write_txt(fname=outfname, outlat=outlat, outlon=outlon, outZ=outVarr)
        return
            
                    
        
        
        
