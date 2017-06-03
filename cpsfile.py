import numpy as np
import os, copy

def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
class DistFile(object):
    """
    An object for distance file of Computer Programs in Seismology.
    ========================================================================
    Parameters:
    distArr     - distance for origin point
    dtArr       - sampling interval in synthetic seismograms
    nptsArr     - npts array 
    T0Arr       - time of first sample is T0 + DIST/VRED
    VredArr     - see above
    ========================================================================
    """
    def __init__(self, distfname=None):
        try:
            self.read(distfname)
        except:
            self.distArr= np.array([])
            self.dtArr  = np.array([])
            self.nptsArr= np.array([])
            self.T0Arr  = np.array([])
            self.VredArr= np.array([])
        return
    
    def read(self, distfname):
        """
        Read Distance file 
        DIST DT NPTS T0 VRED
        """
        InArr       = np.loadtxt(distfname)
        self.distArr= InArr[:,0]
        self.dtArr  = InArr[:,1]
        self.nptsArr= InArr[:,2]
        self.T0Arr  = InArr[:,3]
        self.VredArr= InArr[:,4]
        return
    
    def write(self, distfname):
        outArr  = np.append(self.distArr, self.dtArr)
        outArr  = np.append(outArr, self.nptsArr)
        outArr  = np.append(outArr, self.T0Arr)
        outArr  = np.append(outArr, self.VredArr)
        outArr  = outArr.reshape(5, self.distArr.size)
        outArr  = outArr.T
        np.savetxt(distfname, outArr, fmt='%f %f %d %f %f')
        return
    
    def add(self, dist, dt=0.1, N2=14, T0=0.0, Vred=0.0):
        """Add a single distance point
        ============================================================
        Input Parameters:
        dist    - distance for origin point
        dt      - sampling interval in synthetic seismograms
        N2      - NPTS = 2**N2
        T0      - time of first sample is T0 + DIST/VRED
        Vred    - see above
        
        Output:
        self.distArr, dtArr, nptsArr, T0Arr, VredArr
        ============================================================
        """
        self.distArr= np.append(self.distArr, dist)
        self.dtArr  = np.append(self.dtArr, dt)
        self.nptsArr= np.append(self.nptsArr, 2**N2)
        self.T0Arr  = np.append(self.T0Arr, T0)
        self.VredArr= np.append(self.VredArr, Vred)
        return
    
    def addEqualDist(self, dist0, dD, Nd, dt=0.1, N2=14, T0=0.0, Vred=0.0):
        """Add equal distance list.
        ============================================================
        Input Parameters:
        dist0   - distance for origin point
        dD      - distance interval
        Nd      - number of distance point 
        dt      - sampling interval in synthetic seismograms
        N2      - NPTS = 2**N2
        T0      - time of first sample is T0 + DIST/VRED
        Vred    - see above
        
        Output:
        self.distArr, dtArr, nptsArr, T0Arr, VredArr
        ============================================================
        """
        self.distArr= np.append(self.distArr, np.arange(Nd)*dD+dist0 )
        self.dtArr  = np.append(self.dtArr, np.ones(Nd)*dt)
        self.nptsArr= np.append(self.nptsArr, np.ones(Nd)*2**N2)
        self.T0Arr  = np.append(self.T0Arr, np.ones(Nd)*T0)
        self.VredArr= np.append(self.VredArr, np.ones(Nd)*Vred)
        return
    
class DispCurve(object):
    """
    An object to handle single mode dispersion curve 
    ========================================================================
    Parameters:
    period  - period array
    Vph     - phase velocity array
    Vgr     - group velocity array
    energy  - energy integral array
    gamma   - anelastic attenuation coefficient array
    ellip   - Rayleigh wave ellipticity
    header  - header dictionary (type: RAYLEIGH or LOVE , mode: 0, 1, 2...)
    ========================================================================
    """
    def __init__(self, period=np.array([]), Vph=np.array([]), Vgr=np.array([]), energy=np.array([]), gamma=np.array([]),
                 ellip=np.array([]), header={'type': 'N/A', 'mode': -1}):
        self.period = period
        self.Vph    = Vph
        self.Vgr    = Vgr
        self.energy = energy
        self.gamma  = gamma
        self.ellip  = ellip
        self.header = header
        if period.size !=Vph.size and period.size !=Vgr.size:
            raise ValueError('Inconsistent dispersion curve!')
        return
    
    def gethdr(self, instr, verbose=True):
        """
        Get the header information
        """
        strLst=instr.split()
        self.header={'type': strLst[0], 'mode': int(strLst[4])}
        if verbose ==True:
            print 'Read dispersion curve for:', instr
        return
    
    def write(self, outfname, datatype='phase'):
        """
        Write dispersion curve to a txt file
        """
        if datatype=='phase':
            outArr=np.append(self.period, self.Vph)
            outArr=outArr.reshape((2, self.period.size))
            outArr=outArr.T
            np.savetxt(outfname, outArr, fmt='%g')
            print 'Write dispersion curve for:', self.header['type'], 'mode',self.header['mode']
        else:
            outArr=np.append(self.period, self.Vgr)
            outArr=outArr.reshape((2, self.period.size))
            outArr=outArr.T
            np.savetxt(outfname, outArr, fmt='%g')
            print 'Write dispersion curve for:', self.header['type'], 'mode',self.header['mode']
        return
    
    def InterpDisp(self, T0=5., dT=1., NT=155):
        """
        Interpolate dispersion curve, 5 ~ 159 second
        """
        Tinterp     = T0+np.arange(NT)*dT
        if self.Vph.size == self.period.size: self.Vph      = np.interp(Tinterp, self.period, self.Vph)
        if self.Vgr.size == self.period.size: self.Vgr      = np.interp(Tinterp, self.period, self.Vgr)
        if self.gamma.size == self.period.size: self.gamma  = np.interp(Tinterp, self.period, self.gamma)
        self.period = Tinterp
        return
        
class DispFile(object):
    """
    An object to handle a list of dispersion curves
    ========================================================================
    Parameters:
    DispLst - list of DispCurve object
    ========================================================================
    """
    def __init__(self, dispfname=None):
        self.DispLst    = []
        if os.path.isfile(dispfname):
            self.read(dispfname)
        return
    
    def get_tree(self, intree=None):
        """
        Get tree for ASDF database
        """
        if intree==None: outtree = {'ray':{}, 'love': {}}
        else: outtree = intree
        for disp in self.DispLst:
            if disp.period.size == 0 or disp.Vph.size == 0: continue
            s_tree = {'T': disp.period, 'Vph': disp.Vph}
            if disp.Vgr.size == disp.period.size:
                s_tree.update({'Vgr': disp.Vgr})
            if disp.energy.size == disp.period.size:
                s_tree.update({'energy': disp.energy})
            if disp.gamma.size == disp.period.size:
                s_tree.update({'gamma': disp.gamma})
            if disp.ellip.size == disp.period.size:
                s_tree.update({'ellip': disp.ellip})
            if disp.header['type'] == 'RAYLEIGH':
                outtree['ray'].update({disp.header['mode']: s_tree})
            elif disp.header['type'] == 'LOVE':
                outtree['love'].update({disp.header['mode']: s_tree})
        return outtree
        
    
    def read(self, dispfname=None):
        """
        read dispersion data in TXT format
        """
        with open(dispfname, 'r') as f:
            for line in f.readlines():
                cline   = line.split()
                if len(cline)==0: continue
                if len(cline)==5:
                    try:
                        self.DispLst.append(copy.deepcopy(dispcurve))
                        dispcurve   = DispCurve()
                        dispcurve.gethdr(line)
                    except:
                        dispcurve   = DispCurve()
                        dispcurve.gethdr(line)
                    continue
                if is_int(cline[0]):
                    try:
                        if len(cline)==3:
                            dispcurve.period    = np.append( dispcurve.period, float(cline[1]) )
                            dispcurve.Vph       = np.append( dispcurve.Vph, float(cline[2]) )
                        if len(cline)==8:
                            dispcurve.period    = np.append( dispcurve.period, float(cline[1]) )
                            dispcurve.Vph       = np.append( dispcurve.Vph, float(cline[3]) )
                            dispcurve.Vgr       = np.append( dispcurve.Vgr, float(cline[4]) )
                            dispcurve.energy    = np.append( dispcurve.energy, float(cline[5]) )
                            dispcurve.gamma     = np.append( dispcurve.gamma, float(cline[6]) )
                            dispcurve.ellip     = np.append( dispcurve.ellip, float(cline[7]) )
                        if len(cline)==7:
                            dispcurve.period    = np.append( dispcurve.period, float(cline[1]) )
                            dispcurve.Vph       = np.append( dispcurve.Vph, float(cline[3]) )
                            dispcurve.Vgr       = np.append( dispcurve.Vgr, float(cline[4]) )
                            dispcurve.energy    = np.append( dispcurve.energy, float(cline[5]) )
                            dispcurve.gamma     = np.append( dispcurve.gamma, float(cline[6]) )
                    except ValueError:
                        print 'WARNING: check period:' + cline[1]
                        continue
            self.DispLst.append(dispcurve)
        return
                
    def write(self, outfname, mode=0, T0=5., dT=1., NT=155, datatype='phase' ):
        if T0!=None and dT != None and NT !=None:
            self.DispLst[mode].InterpDisp(T0=T0, dT=dT, NT=NT )
        self.DispLst[mode].write(outfname=outfname, datatype=datatype)
        return
    
class EigenFunc(object):
    """
    An object to handle single period eigenfunction
    ========================================================================
    Parameters:
    -------------------------- eigenfunctions ------------------------------
    ur, tr, uz, tz  - Rayleigh wave eigenfunctions
    ut, tt          - Love wave eigenfunctions
    ------------------------ sensitivity kernels ---------------------------
    dcdh            - layer thickness
    dcda, dcdb      - P/S wave velocity
    dcdr            - density
    dcdav/dcdah     - PV/PH wave velocity
    dcdbv/dcdbh     - SV/SH wave velocity
    dcdn            - eta (eta = F/(A-2L); A, L, F are Love parameters)
    header          - header dictionary
                    (type, mode, T, C, U, energy, gamma, zref)
    ========================================================================
    """
    def __init__(self, ur=np.array([]), tr=np.array([]), uz=np.array([]), tz=np.array([]), ut=np.array([]), tt=np.array([]),
                dcdh=np.array([]), dcda=np.array([]), dcdb=np.array([]), dcdr=np.array([]),
                dcdav=np.array([]), dcdah=np.array([]), dcdbv=np.array([]), dcdbh=np.array([]), dcdn=np.array([]),
                header={'type': 'N/A', 'mode': -1, 'T': -1, 'C': -1, 'U': -1, 'energy': -1, 'gamma': 0.0, 'zref': -1}):
        self.ur     = ur
        self.tr     = tr
        self.uz     = uz
        self.tz     = tz
        self.ut     = ut
        self.tt     = tt
        self.dcdh   = dcdh
        self.dcda   = dcda
        self.dcdb   = dcdb
        self.dcdr   = dcdr
        self.dcdav  = dcdav
        self.dcdah  = dcdah
        self.dcdbv  = dcdbv
        self.dcdbh  = dcdbh
        self.dcdn   = dcdn
        self.header = header
        return
    
    def get_hdr(self, instr):
        """
        Get header information given a string
        """
        if instr[3]=='#':
            self.header['type'] = instr[0]
            self.header['mode'] = instr[4]
            return True
        if instr[0]=='T':
            try:
                self.header['T']    = float(instr[2])
                self.header['C']    = float(instr[5])
                self.header['U']    = float(instr[-1])
                return True
            except:
                raise ValueError('Check: '+' '.join(instr))
        if instr[0]=='AR=' or instr[0]=='AL=':
            try:
                self.header['energy']   = float(instr[1])
                self.header['gamma']    = float(instr[3])
                self.header['zref']     = float(instr[-1])
                return False
            except:
                raise ValueError('Check: '+' '.join(instr))
        raise ValueError('Check: '+' '.join(instr))
            
    
class derFile(object):
    """
    An object to handle a list of eigenfunctions
    ========================================================================
    Parameters:
    eigenLst    - list of DispCurve object
    -------------------------- model parameters ----------------------------
    HArr        - layer thickness
    rhoArr      - density
    VpArr/VsArr - P/S wave velocity
    QaArr/QbArr - P/S wave quality factor
    AArr, CArr, FArr
    LArr, NArr  - Love parameters
    ========================================================================
    """
    def __init__(self, derfname=None):
        self.eigenLst   = []
        self.HArr       = np.array([])
        self.rhoArr     = np.array([])
        self.VpArr      = np.array([])
        self.VsArr      = np.array([])
        self.QaArr      = np.array([])
        self.QbArr      = np.array([])
        self.AArr       = np.array([])
        self.CArr       = np.array([])
        self.FArr       = np.array([])
        self.LArr       = np.array([])
        self.NArr       = np.array([])
        self.isotropic  = None
        if os.path.isfile(derfname):
            self.read(derfname)
        return
    
    def get_tree(self, intree=None, updatemodel=False):
        """
        Get tree for ASDF database
        """
        if intree==None:
            outtree = {'model': {}, 'egn': {'ray':{}, 'love': {} } } 
        else:
            outtree = intree
        if updatemodel or len(outtree['model'].keys()) == 0:
            outtree.update({'model': {'isotropic': self.isotropic, 'H': self.HArr, 'rho': self.rhoArr, 'vs': self.VsArr, 'vp': self.VpArr, \
                'Qa': self.QaArr, 'Qb': self.QbArr, 'A': self.AArr, 'C': self.CArr, 'F': self.FArr, 'L': self.LArr, 'N': self.NArr}})
        for eigenf in self.eigenLst:
            ###
            # Update mode
            ###
            if eigenf.header['type'] == 'RAYLEIGH':
                try:
                    modetree = outtree['egn']['ray'][int(eigenf.header['mode'])]
                except:
                    outtree['egn']['ray'].update({int(eigenf.header['mode']): {}})
                    modetree = outtree['egn']['ray'][int(eigenf.header['mode'])]
            elif eigenf.header['type'] == 'LOVE':
                try:
                    modetree = outtree['egn']['love'][int(eigenf.header['mode'])]
                except:
                    outtree['egn']['love'].update({int(eigenf.header['mode']): {}})
                    modetree = outtree['egn']['love'][int(eigenf.header['mode'])]
            s_tree = {'ur': eigenf.ur, 'tr': eigenf.tr, 'uz': eigenf.uz, 'tz': eigenf.tz, 'ut': eigenf.ut, 'tt': eigenf.tt, \
                'dcdh': eigenf.dcdh, 'dcda': eigenf.dcda, 'dcdb': eigenf.dcdb, 'dcdr': eigenf.dcdr, 'dcdav': eigenf.dcdav, \
                'dcdah': eigenf.dcdah, 'dcdbv': eigenf.dcdbv, 'dcdbh': eigenf.dcdbh, 'dcdn': eigenf.dcdn, 'zref': eigenf.header['zref'],
                'C': eigenf.header['C'], 'U': eigenf.header['U'], 'energy': eigenf.header['energy'], 'gamma': eigenf.header['gamma']}
            modetree.update( {eigenf.header['T'] :  s_tree})
        return outtree
        
    
    def read(self, derfname=None):
        """
        read eigenfunction data in TXT format
        """
        readingmodel= False
        readinghdr  = False
        readingegn  = False
        with open(derfname, 'r') as f:
            for line in f.readlines():
                cline   = line.split()
                if len(cline)==0: continue
                ########
                # Read velocity model
                ########
                if len(cline)==1 and cline[0] == 'Model:':
                    readingmodel    = True
                    continue
                if cline[0] == 'LAYER':
                    if len(cline) == 10: isotropic = False; self.isotropic = False
                    elif len(cline) == 7: isotropic = True; self.isotropic = True
                    else: raise ValueError('Unrecognized model!')
                    continue
                if len(cline)!=10 and len(cline)!=7 and readingmodel: readingmodel = False
                if readingmodel:
                    self.HArr   = np.append(self.HArr, float(cline[1]))
                    if isotropic:
                        self.VpArr  = np.append(self.VpArr, float(cline[2]))
                        self.VsArr  = np.append(self.VsArr, float(cline[3]))
                        self.rhoArr = np.append(self.rhoArr, float(cline[4]))
                        self.QaArr  = np.append(self.QaArr, float(cline[5]))
                        self.QbArr  = np.append(self.QbArr, float(cline[6]))
                    else:
                        self.AArr   = np.append(self.AArr, float(cline[2]))
                        self.CArr   = np.append(self.CArr, float(cline[3]))
                        self.FArr   = np.append(self.FArr, float(cline[4]))
                        self.LArr   = np.append(self.LArr, float(cline[5]))
                        self.NArr   = np.append(self.NArr, float(cline[6]))
                        self.rhoArr = np.append(self.rhoArr, float(cline[7]))
                        self.QaArr  = np.append(self.QaArr, float(cline[8]))
                        self.QbArr  = np.append(self.QbArr, float(cline[9]))
                    continue
                ########
                # Read header information for eigenfunction
                ########
                if len(cline)==5 and cline[3]=='#':
                    readingegn  = False
                    try:
                        self.eigenLst.append(copy.deepcopy(eigenf))
                        eigenf      = EigenFunc()
                        readinghdr  = eigenf.get_hdr(instr=cline)
                    except:
                        eigenf      = EigenFunc()
                        readinghdr  = eigenf.get_hdr(instr=cline)
                    continue
                if readinghdr:
                    readinghdr = eigenf.get_hdr(instr=cline)
                    continue
                ########
                # Read eigenfunction data
                ########
                if cline[0]=='M':
                    readingegn = True
                    continue
                if readingegn:
                    # Isotropic Love wave
                    if len(cline) == 6:
                        if (not isotropic) or (not eigenf.header['type'] == 'LOVE'):
                            raise ValueError('Check: '+' '.join(cline))
                        eigenf.ut   = np.append( eigenf.ut, float(cline[1]) )
                        eigenf.tt   = np.append( eigenf.tt, float(cline[2]) )
                        eigenf.dcdh = np.append( eigenf.dcdh, float(cline[3]) )
                        eigenf.dcdb = np.append( eigenf.dcdb, float(cline[4]) )
                        eigenf.dcdr = np.append( eigenf.dcdr, float(cline[5]) )
                    # TI Love wave
                    if len(cline) == 7:
                        if isotropic or (not eigenf.header['type'] == 'LOVE'):
                            raise ValueError('Check: '+' '.join(cline))
                        eigenf.ut   = np.append( eigenf.ut, float(cline[1]) )
                        eigenf.tt   = np.append( eigenf.tt, float(cline[2]) )
                        eigenf.dcdh = np.append( eigenf.dcdh, float(cline[3]) )
                        eigenf.dcdbv= np.append( eigenf.dcdbv, float(cline[4]) )
                        eigenf.dcdbh= np.append( eigenf.dcdbh, float(cline[5]) )
                        eigenf.dcdr = np.append( eigenf.dcdr, float(cline[6]) )
                    # Isotropic Rayleigh wave
                    if len(cline) == 9:
                        if (not isotropic) or (not eigenf.header['type'] == 'RAYLEIGH'):
                            raise ValueError('Check: '+' '.join(cline))
                        eigenf.ur   = np.append( eigenf.ur, float(cline[1]) )
                        eigenf.tr   = np.append( eigenf.tr, float(cline[2]) )
                        eigenf.uz   = np.append( eigenf.uz, float(cline[3]) )
                        eigenf.tz   = np.append( eigenf.tz, float(cline[4]) )
                        eigenf.dcdh = np.append( eigenf.dcdh, float(cline[5]) )
                        eigenf.dcda = np.append( eigenf.dcda, float(cline[6]) )
                        eigenf.dcdb = np.append( eigenf.dcdb, float(cline[7]) )
                        eigenf.dcdr = np.append( eigenf.dcdr, float(cline[8]) )
                    # TI Rayleigh wave
                    if len(cline) == 12:
                        if isotropic or (not eigenf.header['type'] == 'RAYLEIGH'):
                            raise ValueError('Check: '+' '.join(cline))
                        eigenf.ur   = np.append( eigenf.ur, float(cline[1]) )
                        eigenf.tr   = np.append( eigenf.tr, float(cline[2]) )
                        eigenf.uz   = np.append( eigenf.uz, float(cline[3]) )
                        eigenf.tz   = np.append( eigenf.tz, float(cline[4]) )
                        eigenf.dcdh = np.append( eigenf.dcdh, float(cline[5]) )
                        eigenf.dcdav= np.append( eigenf.dcdav, float(cline[6]) )
                        eigenf.dcdah= np.append( eigenf.dcdah, float(cline[7]) )
                        eigenf.dcdn = np.append( eigenf.dcdn, float(cline[8]) )
                        eigenf.dcdbv= np.append( eigenf.dcdbv, float(cline[9]) )
                        eigenf.dcdbh= np.append( eigenf.dcdbh, float(cline[10]) )
                        eigenf.dcdr = np.append( eigenf.dcdr, float(cline[11]) )
            self.eigenLst.append(copy.deepcopy(eigenf))
        return
                
    # def write(self, outfname, mode=0, T0=5., dT=1., NT=155, datatype='phase' ):
    #     if T0!=None and dT != None and NT !=None:
    #         self.DispLst[mode].InterpDisp(T0=T0, dT=dT, NT=NT )
    #     self.DispLst[mode].write(outfname=outfname, datatype=datatype)
    #     return