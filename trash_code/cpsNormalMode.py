import numpy as np

class DistFile(object):
    def __init__(self, dfname=None):
        try:
            self.ReadDistFile(self, dfname);
        except:
            self.distArr=np.array([]);
            self.dtArr=np.array([]);
            self.nptsArr=np.array([]);
            self.T0Arr=np.array([]);
            self.VredArr=np.array([]);
        return;
    
    def ReadDistFile(self, dfname):
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
