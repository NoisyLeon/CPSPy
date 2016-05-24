
import os
import numpy as np


ak135Arr=np.loadtxt('ak135forvmodel.mod')

class Model1d(object):
    """
    An object to handle input 1d model for Computer Programs in Seismology.
    --------------------------------------------------------------------------------------------------------------
    Parameters:
    modelver           - model version
    modelname       - model name
    modelindex      - index indicating model type
                                1: 'ISOTROPIC', 2: 'TRANSVERSE ISOTROPIC', 3: 'ANISOTROPIC'
    modelunit         - KGS km, km/s, g/cm^3
    earthindex         - index indicating Earth type 1: 'FLAT EARTH', 2:'SPHERICAL EARTH'
    boundaryindex  - index indicating model boundaries 1: '1-D', 2: '2-D', 3: '3-D'
    Vindex              - index indicating nature of layer velocity
    HArr                  - layer thickness array
    VsArr, VpArr, rhoArr, QpArr, QsArr, etapArr, etasArr, frefpArr,  frefsArr
                             - model parameters
    DepthArr          - depth array
    --------------------------------------------------------------------------------------------------------------
    """
    def __init__(self, modelver='MODEL.01', modelname='TEST MODEL', modelindex=1, modelunit='KGS', earthindex=1,
            boundaryindex=1, Vindex=1, HArr=np.array([]), VsArr=np.array([]), VpArr=np.array([]), rhoArr=np.array([]),
            QpArr=np.array([]), QsArr=np.array([]), etapArr=np.array([]), etasArr=np.array([]), frefpArr=np.array([]),  frefsArr=np.array([])):
        self.modelver=modelver
        self.modelname=modelname
        modeldict={1: 'ISOTROPIC', 2: 'TRANSVERSE ISOTROPIC', 3: 'ANISOTROPIC'}
        self.modeltype=modeldict[modelindex]
        self.modelunit=modelunit
        earthdict={1: 'FLAT EARTH', 2:'SPHERICAL EARTH'}
        self.earthtype=earthdict[earthindex]
        boundarydict={1: '1-D', 2: '2-D', 3: '3-D'}
        self.boundarytype=boundarydict[boundaryindex]
        Vdict={1: 'CONSTANT VELOCITY', 2: 'VARIABLE VELOCITY'}
        self.Vtype=Vdict[Vindex]
        self.line08_11='LINE08\nLINE09\nLINE10\nLINE11\n'
        self.modelheader='H   VP  VS RHO  QP  QS ETAP ETAS FREFP FREFS'
        self.HArr=HArr
        self.VsArr=VsArr
        self.VpArr=VpArr
        self.rhoArr=rhoArr
        self.QpArr=QpArr
        self.QsArr=QsArr
        self.etapArr=etapArr
        self.etasArr=etasArr
        self.frefpArr=frefpArr
        self.frefsArr=frefsArr
        self.DepthArr=np.cumsum(self.HArr)
        return
    
    def ak135(self, modelname='AK135 CONTINENTAL MODEL'):
        """
        ak135 model
        """
        self.modelname = modelname
        self.HArr=ak135Arr[:,0]
        self.VpArr=ak135Arr[:,1]
        self.VsArr=ak135Arr[:,2]
        self.rhoArr=ak135Arr[:,3]
        self.QpArr=ak135Arr[:,4]
        self.QsArr=ak135Arr[:,5]
        self.etapArr=ak135Arr[:,6]
        self.etasArr=ak135Arr[:,7]
        self.frefpArr=ak135Arr[:,8]
        self.frefsArr=ak135Arr[:,9]
        self.DepthArr=np.cumsum(self.HArr)
        return

    def addlayer(self, H, vs, vp=None, rho=None, Qp=310., Qs=150., etap=0.0, etas=0.0, frefp=1.0, frefs=1.0,
                zmin=9999.):
        """ Add layer to the model
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        H                   - layer thickness
        vs                  - vs
        vp, rho          - vp, rho default is None by assuming Brocher Crust
        Qp, Qs          - quality factor
        etap, etas, frefp, frefs  - see the manual for computer programs in seismology
        zmin             -   top depth of the layer
                                1. default is 9999, which will simply append one layer to the model
                                2. if zmin < the bottom of preexisting top layer, it will append a new
                                    layer to the top (also replace some part of the preexisting top layer)
                                3. else, the code will replace some part of preexisting model( this part need further checking! )
        -----------------------------------------------------------------------------------------------------
        """
        if vp ==None:
            vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
        if rho==None:
            rho=1.6612*vp-0.4721*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**5
        if zmin >= self.DepthArr[-1]:
            self.HArr=np.append(self.HArr, H)
            self.VpArr=np.append(self.VpArr, vp)
            self.VsArr=np.append(self.VsArr, vs)
            self.rhoArr=np.append(self.rhoArr, rho)
            self.QpArr=np.append(self.QpArr, Qp)
            self.QsArr=np.append(self.QsArr, Qs)
            self.etapArr=np.append(self.etapArr, etap)
            self.etasArr=np.append(self.etasArr, etas)
            self.frefpArr=np.append(self.frefpArr, frefp)
            self.frefsArr=np.append(self.frefsArr, frefs)
            self.DepthArr=np.cumsum(self.HArr)
        elif zmin < self.DepthArr[0]:
            # discard layers with depth < H
            self.HArr=self.HArr[self.DepthArr >H]
            self.VpArr=self.VpArr[self.DepthArr >H]
            self.VsArr=self.VsArr[self.DepthArr >H]
            self.rhoArr=self.rhoArr[self.DepthArr >H]
            self.QpArr=self.QpArr[self.DepthArr >H]
            self.QsArr=self.QsArr[self.DepthArr >H]
            self.etapArr=self.etapArr[self.DepthArr >H]
            self.etasArr=self.etasArr[self.DepthArr >H]
            self.frefpArr=self.frefpArr[self.DepthArr >H]
            self.frefsArr=self.frefsArr[self.DepthArr >H]
            # change the thickness of the current first layer
            self.HArr[0]=(self.DepthArr[self.DepthArr >H])[0]-H
            # add H to the first layer
            self.HArr=np.append(H, self.HArr)
            self.VpArr=np.append(vp, self.VpArr)
            self.VsArr=np.append(vs, self.VsArr)
            self.rhoArr=np.append(rho, self.rhoArr)
            self.QpArr=np.append(Qp, self.QpArr)
            self.QsArr=np.append(Qs, self.QsArr)
            self.etapArr=np.append(etap, self.etapArr)
            self.etasArr=np.append(etas, self.etasArr)
            self.frefpArr=np.append(frefp, self.frefpArr)
            self.frefsArr=np.append(frefs, self.frefsArr)
            self.DepthArr=np.cumsum(self.HArr)
        else:
            zmax=zmin+H
            topArr=self.DepthArr-self.HArr
            # Top layers above zmin
            HArrT=self.HArr[self.DepthArr <=zmin]
            VpArrT=self.VpArr[self.DepthArr <=zmin]
            VsArrT=self.VsArr[self.DepthArr <=zmin]
            rhoArrT=self.rhoArr[self.DepthArr <=zmin]
            QpArrT=self.QpArr[self.DepthArr <=zmin]
            QsArrT=self.QsArr[self.DepthArr <=zmin]
            etapArrT=self.etapArr[self.DepthArr <=zmin]
            etasArrT=self.etasArr[self.DepthArr <=zmin]
            frefpArrT=self.frefpArr[self.DepthArr <=zmin]
            frefsArrT=self.frefsArr[self.DepthArr <=zmin]
            
            tH= zmin - (topArr[self.DepthArr >zmin]) [0]
            
            if tH !=0:
                tVp=(self.VpArr[self.DepthArr >zmin])[0]
                tVs=(self.VsArr[self.DepthArr >zmin])[0]
                trho=(self.rhoArr[self.DepthArr >zmin])[0]
                tQp=(self.QpArr[self.DepthArr >zmin])[0]
                tQs=(self.QsArr[self.DepthArr >zmin])[0]
                tetap=(self.etapArr[self.DepthArr >zmin])[0]
                tetas=(self.etasArr[self.DepthArr >zmin])[0]
                tfrefp=(self.frefpArr[self.DepthArr >zmin])[0]
                tfrefs=(self.frefsArr[self.DepthArr >zmin])[0]
                HArrT=np.append(HArrT, tH)
                VpArrT=np.append(VpArrT, tVp)
                VsArrT=np.append(VsArrT, tVs)
                rhoArrT=np.append(rhoArrT, trho)
                QpArrT=np.append(QpArrT, tQp)
                QsArrT=np.append(QsArrT, tQs)
                etapArrT=np.append(etapArrT, tetap)
                etasArrT=np.append(etasArrT, tetas)
                frefpArrT=np.append(frefpArrT, tfrefp)
                frefsArrT=np.append(frefsArrT, tfrefs)
            
            # Bottom layer bolow zmax
            HArrB=self.HArr[topArr >=zmax]
            VpArrB=self.VpArr[topArr >=zmax]
            VsArrB=self.VsArr[topArr >=zmax]
            rhoArrB=self.rhoArr[topArr >=zmax]
            QpArrB=self.QpArr[topArr >=zmax]
            QsArrB=self.QsArr[topArr >=zmax]
            etapArrB=self.etapArr[topArr >=zmax]
            etasArrB=self.etasArr[topArr >=zmax]
            frefpArrB=self.frefpArr[topArr >=zmax]
            frefsArrB=self.frefsArr[topArr >=zmax]
            bH=  (self.DepthArr[topArr <zmax]) [-1]- zmax
            print HArrT, tH, bH
            if bH !=0:
                bVp=(self.VpArr[topArr <zmax]) [-1]
                bVs=(self.VsArr[topArr <zmax]) [-1]
                brho=(self.rhoArr[topArr <zmax]) [-1]
                bQp=(self.QpArr[topArr <zmax]) [-1]
                bQs=(self.QsArr[topArr <zmax]) [-1]
                betap=(self.etapArr[topArr <zmax]) [-1]
                betas=(self.etasArr[topArr <zmax]) [-1]
                bfrefp=(self.frefpArr[topArr <zmax]) [-1]
                bfrefs=(self.frefsArr[topArr <zmax]) [-1]
                HArrB=np.append(bH, HArrB)
                VpArrB=np.append(bVp, VpArrB)
                VsArrB=np.append(bVs, VsArrB)
                rhoArrB=np.append(brho, rhoArrB)
                QpArrB=np.append(bQp, QpArrB)
                QsArrB=np.append(bQs, QsArrB)
                etapArrB=np.append(betap, etapArrB)
                etasArrB=np.append(betas, etasArrB)
                frefpArrB=np.append(bfrefp, frefpArrB)
                frefsArrB=np.append(bfrefs, frefsArrB)
            #####################################
            self.HArr=np.append(HArrT, H)
            self.VpArr=np.append(VpArrT, vp)
            self.VsArr=np.append(VsArrT, vs)
            self.rhoArr=np.append(rhoArrT, rho)
            self.QpArr=np.append(QpArrT, Qp)
            self.QsArr=np.append(QsArrT, Qs)
            self.etapArr=np.append(etapArrT, etap)
            self.etasArr=np.append(etasArrT, etas)
            self.frefpArr=np.append(frefpArrT, frefp)
            self.frefsArr=np.append(frefsArrT,frefs)
            #####################################
            self.HArr=np.append(self.HArr, HArrB)
            self.VpArr=np.append(self.VpArr, VpArrB)
            self.VsArr=np.append(self.VsArr, VsArrB)
            self.rhoArr=np.append(self.rhoArr, rhoArrB)
            self.QpArr=np.append(self.QpArr, QpArrB)
            self.QsArr=np.append(self.QsArr, QsArrB)
            self.etapArr=np.append(self.etapArr, etapArrB)
            self.etasArr=np.append(self.etasArr, etasArrB)
            self.frefpArr=np.append(self.frefpArr, frefpArrB)
            self.frefsArr=np.append(self.frefsArr, frefsArrB)
            self.DepthArr=np.cumsum(self.HArr)
        return


    def write(self, outfname, fmt='%g'):
        """
        Write profile to the Computer Programs in Seismology model format
        """
        with open(outfname, 'w') as f:
            f.write(self.modelver+'\n')
            f.write(self.modelname+'\n')
            f.write(self.modeltype+'\n')
            f.write(self.modelunit+'\n')
            f.write(self.earthtype+'\n')
            f.write(self.boundarytype+'\n')
            f.write(self.Vtype+'\n')
            f.write(self.line08_11)
            f.write(self.modelheader+'\n')
            for i in np.arange(self.HArr.size):
                tempstr='%f %f %f %f %f %f %f %f %f %f \n' %(self.HArr[i], self.VpArr[i], self.VsArr[i], self.rhoArr[i],
                        self.QpArr[i], self.QsArr[i], self.etapArr[i], self.etasArr[i], self.frefpArr[i], self.frefsArr[i])
                # tempstr=fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt+' '+fmt + ' \n' %(self.HArr[i], self.VpArr[i], self.VsArr[i], self.rhoArr[i],
                #         self.QpArr[i], self.QsArr[i], self.etapArr[i], self.etasArr[i], self.frefpArr[i], self.frefsArr[i])
                f.write(tempstr)
        return
            
    def read(self, infname, sep='\t'):
        """
        Read Computer Programs in Seismology model format
        """
        with open(infname, 'r') as f:
            self.modelver=(f.readline()).split('\n')[0]
            self.modelname=(f.readline()).split('\n')[0]
            self.modeltype=(f.readline()).split('\n')[0]
            self.modelunit=(f.readline()).split('\n')[0]
            self.earthtype=(f.readline()).split('\n')[0]
            self.boundarytype=(f.readline()).split('\n')[0]
            self.Vtype=(f.readline()).split('\n')[0]
            f.readline()
            f.readline()
            f.readline()
            f.readline()
            self.modelheader=(f.readline()).split('\n')[0]
            for line in f.readlines():
                cline=line.split()
                self.HArr=np.append(self.HArr, float(cline[0]))
                self.VpArr=np.append(self.VpArr, float(cline[1]))
                self.VsArr=np.append(self.VsArr, float(cline[2]))
                self.rhoArr=np.append(self.rhoArr, float(cline[3]))
                self.QpArr=np.append(self.QpArr, float(cline[4]))
                self.QsArr=np.append(self.QsArr, float(cline[5]))
                self.etapArr=np.append(self.etapArr, float(cline[6]))
                self.etasArr=np.append(self.etasArr, float(cline[7]))
                self.frefpArr=np.append(self.frefpArr, float(cline[8]))
                self.frefsArr=np.append(self.frefsArr, float(cline[9]))
        return
    
    def perturb(self, dm, zmin=-9999, zmax=9999, datatype='vs'):
        """ Add perturbation to the model given a depth interval
        -----------------------------------------------------------------------------------------------------
        Input Parameters:
        dm                - perturbed value (-1, 1)
        zmin, zmax   - 
        datatype       - data type for perturbation
        -----------------------------------------------------------------------------------------------------
        """
        index=(self.DepthArr<zmax) * (self.DepthArr>zmin)
        # print index
        if datatype=='vp':
            self.VpArr[index]=self.VpArr[index]*(1.+dm)
        if datatype=='vs':
            # print self.VsArr[index]*(1.+dm)
            self.VsArr[index]=self.VsArr[index]*(1.+dm)
        if datatype=='rho':
            self.rhoArr[index]=self.rhoArr[index]*(1.+dm)
        if datatype=='qp':
            self.QpArr[index]=self.QpArr[index]*(1.+dm)
        if datatype=='qs':
            self.QsArr[index]=self.QsArr[index]*(1.+dm)
        if datatype=='etap':
            self.etapArr[index]=self.etapArr[index]*(1.+dm)
        if datatype=='etas':
            self.etasArr[index]=self.etasArr[index]*(1.+dm)
        if datatype=='frefp':
            self.frefpArr[index]=self.frefpArr[index]*(1.+dm)
        if datatype=='frefs':
            self.frefsArr[index]=self.frefsArr[index]*(1.+dm)
        ztop=self.DepthArr[index][0]
        zbottom=self.DepthArr[index][-1]
        print 'Top:', ztop, 'km Bottom:', zbottom,'km'
        return
    
    
    
    

        
    
 
