import numpy as np

class DistFile(object):
    def __init__(self, dfname=None):
        try:
            self.read(dfname);
        except:
            self.distArr=np.array([]);
            self.dtArr=np.array([]);
            self.nptsArr=np.array([]);
            self.T0Arr=np.array([]);
            self.VredArr=np.array([]);
        return;
    
    def read(self, dfname):
        """
        Read Distance file 
        DIST DT NPTS T0 VRED
        """
        InArr=np.loadtxt(dfname);
        self.distArr=InArr[:,0];
        self.dtArr=InArr[:,1];
        self.nptsArr=InArr[:,2];
        self.T0Arr=InArr[:,3];
        self.VredArr=InArr[:,4];
        return
    
    def write(self, dfname):
        outArr=np.append(self.distArr, self.dtArr);
        outArr=np.append(outArr, self.nptsArr);
        outArr=np.append(outArr, self.T0Arr);
        outArr=np.append(outArr, self.VredArr);
        outArr=outArr.reshape(5, self.distArr.size);
        outArr=outArr.T;
        np.savetxt(dfname, outArr, fmt='%f %f %d %f %f');
        return;
    
    def add(self, dist, dt=0.1, N2=14, T0=-1.0, Vred=6.0):
        self.distArr=np.append(self.distArr, dist);
        self.dtArr=np.append(self.dtArr, dt);
        self.nptsArr=np.append(self.nptsArr, 2**N2);
        self.T0Arr=np.append(self.T0Arr, T0);
        self.VredArr=np.append(self.VredArr, Vred);
        return;
    
    def addEqualDist(self, dist0, dD, Nd, dt=0.1, N2=14, T0=-1.0, Vred=6.0):
        self.distArr=np.append(self.distArr, np.arange(Nd)*dD+dist0 );
        self.dtArr=np.append(self.dtArr, np.ones(Nd)*dt);
        self.nptsArr=np.append(self.nptsArr, np.ones(Nd)*2**N2);
        self.T0Arr=np.append(self.T0Arr, np.ones(Nd)*T0);
        self.VredArr=np.append(self.VredArr, np.ones(Nd)*Vred);
        return;
    
    
    
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
    