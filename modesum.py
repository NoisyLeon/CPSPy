
import numpy as np
import pyasdf
import vmodel

class modesumASDF(pyasdf.ASDFDataSet):
    
    def getmodel(self, inmodel=None):
        if isinstance(inmodel, vmodel.Model1d):
            self.model1d    = inmodel
        else:
            self.model1d    = vmodel.Model1d()
            self.model1d.ak135()
        return
    
    def getdfile(self, distfile=None):
        if isinstance(distfile, cpsfile.DistFile()):
            self.distfile   = distfile
        else os.path.isfile(distfile):
            self.distfile   = cpsfile.DistFile(distfname=distfile)
        return
    
    def run_prep96(self, workingdir):
        