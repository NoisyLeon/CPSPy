import numpy as np
import os
def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
class DistFile(object):
    def __init__(self, distfname=None):
        try:
            self.read(distfname);
        except:
            self.distArr=np.array([]);
            self.dtArr=np.array([]);
            self.nptsArr=np.array([]);
            self.T0Arr=np.array([]);
            self.VredArr=np.array([]);
        return;
    
    def read(self, distfname):
        """
        Read Distance file 
        DIST DT NPTS T0 VRED
        """
        InArr=np.loadtxt(distfname);
        self.distArr=InArr[:,0];
        self.dtArr=InArr[:,1];
        self.nptsArr=InArr[:,2];
        self.T0Arr=InArr[:,3];
        self.VredArr=InArr[:,4];
        return
    
    def write(self, distfname):
        outArr=np.append(self.distArr, self.dtArr);
        outArr=np.append(outArr, self.nptsArr);
        outArr=np.append(outArr, self.T0Arr);
        outArr=np.append(outArr, self.VredArr);
        outArr=outArr.reshape(5, self.distArr.size);
        outArr=outArr.T;
        np.savetxt(distfname, outArr, fmt='%f %f %d %f %f');
        return;
    
    def add(self, dist, dt=0.1, N2=14, T0=0.0, Vred=0.0):
        """Add a single distance point
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        dist     - distance for origin point
        dt         - sampling interval in synthetic seismograms
        N2        - NPTS = 2**N2
        T0         - time of first sample is T0 + DIST/VRED
        Vred     - see above
        
        Output:
        self.distArr, dtArr, nptsArr, T0Arr, VredArr
        -----------------------------------------------------------------------------------------------------
        """
        self.distArr=np.append(self.distArr, dist);
        self.dtArr=np.append(self.dtArr, dt);
        self.nptsArr=np.append(self.nptsArr, 2**N2);
        self.T0Arr=np.append(self.T0Arr, T0);
        self.VredArr=np.append(self.VredArr, Vred);
        return;
    
    def addEqualDist(self, dist0, dD, Nd, dt=0.1, N2=14, T0=0.0, Vred=0.0):
        """Add equal distance list.
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        dist0    - distance for origin point
        dD       - distance interval
        Nd       - number of distance point 
        dt         - sampling interval in synthetic seismograms
        N2        - NPTS = 2**N2
        T0         - time of first sample is T0 + DIST/VRED
        Vred     - see above
        
        Output:
        self.distArr, dtArr, nptsArr, T0Arr, VredArr
        -----------------------------------------------------------------------------------------------------
        """
        self.distArr=np.append(self.distArr, np.arange(Nd)*dD+dist0 );
        self.dtArr=np.append(self.dtArr, np.ones(Nd)*dt);
        self.nptsArr=np.append(self.nptsArr, np.ones(Nd)*2**N2);
        self.T0Arr=np.append(self.T0Arr, np.ones(Nd)*T0);
        self.VredArr=np.append(self.VredArr, np.ones(Nd)*Vred);
        return;
    
class DispCurve(object):
    def __init__(self, period=np.array([]), Vph=np.array([]), Vgr=np.array([]), header={'type': 'N/A', 'mode': -1}):
        self.period=period
        self.Vph=Vph
        self.Vgr=Vgr
        self.header=header
        if period.size !=Vph.size and period.size !=Vgr.size:
            raise ValueError('Inconsistent dispersion curve!')
        return
    
    def gethdr(self, instr, verbose=True):
        strLst=instr.split()
        self.header={'type': strLst[0], 'mode': int(strLst[4])}
        if verbose ==True:
            print 'Read dispersion curve for:', instr
        return
    
    def write(self, outfname, datatype='phase'):
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
        Tinterp=T0+np.arange(NT)*dT
        if self.Vph.size == self.period.size:
            self.Vph=np.interp(Tinterp, self.period, self.Vph)
        if self.Vgr.size == self.period.size:
            self.Vgr=np.interp(Tinterp, self.period, self.Vgr)
        self.period=Tinterp
        return
        
class DispFile(object):
    def __init__(self, dispfname=None):
        self.DispLst={}
        # self.ModeLst=np.array([])
        if os.path.isfile(dispfname):
            self.read(dispfname)
        return
    
    def read(self, dispfname=None):
        with open(dispfname, 'r') as f:
            for line in f.readlines():
                cline=line.split()
                if len(cline)==0:
                    continue
                if len(cline)==5:
                    try:
                        self.DispLst[dispcurve.header['mode']]=dispcurve
                        dispcurve=DispCurve()
                        dispcurve.gethdr(line)
                    except:
                        dispcurve=DispCurve()
                        dispcurve.gethdr(line)
                    continue
                if is_int(cline[0]):
                    dispcurve.period=np.append( dispcurve.period, float(cline[1]) )
                    dispcurve.Vph=np.append( dispcurve.Vph, float(cline[2]) )
            self.DispLst[dispcurve.header['mode']]=dispcurve
        return
                
    def write(self, outfname, mode=0, T0=5., dT=1., NT=155 ):
        if T0!=None and dT != None and NT !=None:
            self.DispLst[mode].InterpDisp(T0=T0, dT=dT, NT=NT )
        self.DispLst[mode].write(outfname)
        return
        
                    
    
    

    

    
    
    
class StaInfo(object):
    """
    An object contains a station information several methods for station related analysis.
    -----------------------------------------------------------------------------------------------------
    General Parameters:
    stacode     - station name
    network     - network
    chan        - channels for analysis
    x,y     - position for station
    dist    - distance
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self, stacode=None, network='CPS', x=None, y=None, dist=None, npts=None,dt=None, T0=None, Vred=None):

        self.stacode=stacode;
        self.network=network;
        self.x=x;
        self.y=y;
        self.dist=dist;
        self.npts=npts;
        self.dt=dt;
        self.T0=T0;
        self.Vred=Vred;
        return;
    
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of Traces in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self
    def ReadDistFile(self, dfname):
        """
        Read Distance file 
        DIST DT NPTS T0 VRED
        """
        InArr=np.loadtxt(dfname);
        distArr=InArr[:,0];
        dtArr=InArr[:,1];
        nptsArr=InArr[:,2];
        T0Arr=InArr[:,3];
        VredArr=InArr[:,4];
        for i in np.arange(distArr.size):
            self.append(StaInfo (stacode=stacode, dist=distArr[i], npts=nptsArr[i], dt=dtArr[i], T0=T0Arr[i], Vred=VredArr[i] ))
        
        return
    
    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'r')
        Sta=[]
        for lines in f.readlines():
            lines=lines.split()
            stacode=lines[0]
            network=lines[1]
            x=float(lines[2])
            z=float(lines[3])
            # if len(lines)==5:
            #     try:
            #         ccflag=int(lines[3])
            #         network=lines[4]
            #     except ValueError:
            #         ccflag=int(lines[4])
            #         network=lines[3]
            # if len(lines)==4:
            #     try:
            #         ccflag=int(lines[3])
            #     except ValueError:
            #         network=lines[3]
            netsta=network+'.'+stacode
            if Sta.__contains__(netsta):
                index=Sta.index(netsta)
                if abs(self[index].lon-lon) >0.01 and abs(self[index].lat-lat) >0.01:
                    raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                else:
                    print 'Warning: Repeated Station:' +netsta+' in Station List!'
                    continue
            Sta.append(netsta)
            self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
            f.close()
        return
    
    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'r')
        Sta=[]
        for lines in f.readlines():
            lines=lines.split()
            stacode=lines[0]
            network=lines[1]
            x=float(lines[2])
            z=float(lines[3])
            # if len(lines)==5:
            #     try:
            #         ccflag=int(lines[3])
            #         network=lines[4]
            #     except ValueError:
            #         ccflag=int(lines[4])
            #         network=lines[3]
            # if len(lines)==4:
            #     try:
            #         ccflag=int(lines[3])
            #     except ValueError:
            #         network=lines[3]
            netsta=network+'.'+stacode
            if Sta.__contains__(netsta):
                index=Sta.index(netsta)
                if abs(self[index].lon-lon) >0.01 and abs(self[index].lat-lat) >0.01:
                    raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                else:
                    print 'Warning: Repeated Station:' +netsta+' in Station List!'
                    continue
            Sta.append(netsta)
            self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ))
            f.close()
        return
    
    def WriteStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network x z
        """
        f = open(stafile, 'w')
        for sta in self.stations:
                tempstr='%s %s %1.5e %1.5e 0.0 0.0  \n' %( sta.stacode, sta.network, sta.x, sta.z )
                f.writelines(tempstr)
        f.close()
        return
    
    def AddSingle(self, x, z, stacode=None, network='MEM2D'):
        if stacode==None:
            stacode=str(int(x))+'S'+str(int(z));
        self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ));
        return;
    
    def HomoStaLst(self, xmin, xmax, dx, zmin, zmax, dz, network='MEM2D'):
        Nx = int((xmax-xmin)/dx)+1;
        Nz = int((zmax-zmin)/dz)+1;
        for i in np.arange(Nx):
            for j in np.arange(Nz):
                x=xmin+dx*i;
                z=zmin+dz*j;
                stacode=str(int(i))+'S'+str(int(j));
                self.append(StaInfo (stacode=stacode, network=network, x=x, z=z ));
        return;
    