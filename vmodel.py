
import os
import numpy as np
import matplotlib.pyplot as plt
import copy


class Model1d(object):
    """
    An object to handle input 1d model for Computer Programs in Seismology.
    ============================================================================================
    Parameters:
    modelver        - model version
    modelname       - model name
    modelindex      - index indicating model type
                        1: 'ISOTROPIC', 2: 'TRANSVERSE ISOTROPIC', 3: 'ANISOTROPIC'
    modelunit       - KGS km, km/s, g/cm^3
    earthindex      - index indicating Earth type 1: 'FLAT EARTH', 2:'SPHERICAL EARTH'
    boundaryindex   - index indicating model boundaries 1: '1-D', 2: '2-D', 3: '3-D'
    Vindex          - index indicating nature of layer velocity
    HArr            - layer thickness array
    VsArr, VpArr, rhoArr, QpArr, QsArr, etapArr, etasArr, frefpArr,  frefsArr
                    - model parameters
    DepthArr        - depth array
    ============================================================================================
    """
    def __init__(self, modelver='MODEL.01', modelname='TEST MODEL', modelindex=1, modelunit='KGS', earthindex=1,
            boundaryindex=1, Vindex=1, HArr=np.array([]), VsArr=np.array([]), VpArr=np.array([]), rhoArr=np.array([]),
            QpArr=np.array([]), QsArr=np.array([]), etapArr=np.array([]), etasArr=np.array([]), frefpArr=np.array([]),  frefsArr=np.array([])):
        self.modelver   = modelver
        self.modelname  = modelname
        modeldict       = {1: 'ISOTROPIC', 2: 'TRANSVERSE ISOTROPIC', 3: 'ANISOTROPIC'}
        self.modeltype  = modeldict[modelindex]
        self.modelunit  = modelunit
        earthdict       = {1: 'FLAT EARTH', 2:'SPHERICAL EARTH'}
        self.earthtype  = earthdict[earthindex]
        boundarydict    = {1: '1-D', 2: '2-D', 3: '3-D'}
        self.boundarytype=boundarydict[boundaryindex]
        Vdict           = {1: 'CONSTANT VELOCITY', 2: 'VARIABLE VELOCITY'}
        self.Vtype      = Vdict[Vindex]
        self.line08_11  = 'LINE08\nLINE09\nLINE10\nLINE11\n'
        self.modelheader= 'H   VP  VS RHO  QP  QS ETAP ETAS FREFP FREFS'
        self.HArr       = HArr
        self.VsArr      = VsArr
        self.VpArr      = VpArr
        self.rhoArr     = rhoArr
        self.QpArr      = QpArr
        self.QsArr      = QsArr
        self.etapArr    = etapArr
        self.etasArr    = etasArr
        self.frefpArr   = frefpArr
        self.frefsArr   = frefsArr
        self.DepthArr   = np.cumsum(self.HArr)
        return
    
    def copy(self):
        return copy.deepcopy(self)
    
    def ak135(self, modelname='AK135 CONTINENTAL MODEL'):
        """
        ak135 model
        """
        self.modelname  = modelname
        ak135Arr        = np.loadtxt('ak135forvmodel.mod')
        self.HArr       = ak135Arr[:,0]
        self.VpArr      = ak135Arr[:,1]
        self.VsArr      = ak135Arr[:,2]
        self.rhoArr     = ak135Arr[:,3]
        self.QpArr      = ak135Arr[:,4]
        self.QsArr      = ak135Arr[:,5]
        self.etapArr    = ak135Arr[:,6]
        self.etasArr    = ak135Arr[:,7]
        self.frefpArr   = ak135Arr[:,8]
        self.frefsArr   = ak135Arr[:,9]
        self.DepthArr   = np.cumsum(self.HArr)
        return
    
    def getmodel(self, modelname, HArr, VpArr, VsArr, rhoArr, QpArr, QsArr, etap=0., etas=0., frefp=1.0, fres=1.):
        """
        get model
        """
        self.modelname  = modelname
        self.HArr       = HArr
        self.VpArr      = VpArr
        self.VsArr      = VsArr
        self.rhoArr     = rhoArr
        self.QpArr      = QpArr
        self.QsArr      = QsArr
        self.etapArr    = etap*np.ones(HArr.size)
        self.etasArr    = etas*np.ones(HArr.size)
        self.frefpArr   = frefp*np.ones(HArr.size)
        self.frefsArr   = fres*np.ones(HArr.size)
        self.DepthArr   = np.cumsum(self.HArr)
        return
    
    def check_model(self):
        if (self.HArr[self.HArr<1.0]).size %2 !=0: raise ValueError('Invalid vertical profile! Check layer thicknesses!')
        tempH   = self.HArr[self.HArr<1.0]
        tempVs  = self.VsArr[self.HArr<1.0]
        tempVp  = self.VpArr[self.HArr<1.0]
        temprho = self.rhoArr[self.HArr<1.0]
        tempQs  = self.QsArr[self.HArr<1.0]
        tempQp  = self.QpArr[self.HArr<1.0]
        ind0    = np.arange(tempH.size/2, dtype=int)*2
        ind1    = ind0 + 1
        HArr    = np.append( (tempH[ind0] + tempH[ind1]), self.HArr[self.HArr>=1.0])
        VsArr   = np.append( (tempVs[ind0] + tempVs[ind1])/2., self.VsArr[self.HArr>=1.0])
        VpArr   = np.append( (tempVp[ind0] + tempVp[ind1])/2., self.VpArr[self.HArr>=1.0])
        rhoArr  = np.append( (temprho[ind0] + temprho[ind1])/2., self.rhoArr[self.HArr>=1.0])
        QsArr   = np.append( (tempQs[ind0] + tempQs[ind1])/2., self.QsArr[self.HArr>=1.0])
        QpArr   = np.append( (tempQp[ind0] + tempQp[ind1])/2., self.QpArr[self.HArr>=1.0])
        if HArr.size <= 200:
            self.getmodel(modelname=self.modelname, HArr=HArr, VpArr=VpArr, VsArr=VsArr,
                        rhoArr=rhoArr, QpArr=QpArr, QsArr=QsArr)
        else:
            self.getmodel(modelname=self.modelname, HArr=HArr[:200], VpArr=VpArr[:200], VsArr=VsArr[:200],
                        rhoArr=rhoArr[:200], QpArr=QpArr[:200], QsArr=QsArr[:200])
        
        return
        
        
        

    def addlayer(self, H, vs, vp=None, rho=None, Qp=310., Qs=150., etap=0.0, etas=0.0, frefp=1.0, frefs=1.0,
                zmin=9999.):
        """ Add layer to the model
        =========================================================================================================================
        Input Parameters:
        H               - layer thickness
        vs              - vs
        vp, rho         - vp, rho default is None by assuming Brocher Crust
        Qp, Qs          - quality factor
        etap, etas, frefp, frefs
                        - see the manual for computer programs in seismology
        zmin            -   top depth of the layer
                            1. default is 9999, which will simply append one layer to the model
                            2. if zmin < the bottom of preexisting top layer, it will append a new
                                layer to the top (also replace some part of the preexisting top layer)
                            3. else, the code will replace some part of preexisting model( this part need further checking! )
        =========================================================================================================================
        """
        if vp  is None:
            vp=0.9409+2.0947*vs-0.8206*vs**2+0.2683*vs**3-0.0251*vs**4
        if rho is None:
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
        elif zmin <= 0.:
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


    def write(self, outfname):
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
    
    def read_layer_txt(self, infname):
        vpind=None; rhoind=None
        with open(infname) as f:
            inline = f.readline()
            if inline.split()[0] =='#':
                c=inline.split()[1:]
            c=np.array(c)
            z0ind = np.where(c=='z0')[0][0]
            Hind  = np.where(c=='H')[0][0]
            vsind = np.where(c=='vs')[0][0]
            try: vpind  = np.where(c=='vp')[0][0]
            except: pass
            try: rhoind = np.where(c=='rho')[0][0]
            except: pass
        inArr = np.loadtxt(infname)
        z0Arr = inArr[:, z0ind]
        HArr  = inArr[:, Hind]
        vsArr = inArr[:, vsind]
        # if vsind != None: vsArr = inArr[:, vsind]
        if vpind != None: vpArr = inArr[:, vpind]
        if rhoind != None: rhoArr = inArr[:, rhoind]
        for i in xrange(z0Arr.size):
            z0 = z0Arr[i]; H = HArr[i]; vs=vsArr[i]
            if vpind==None and rhoind ==None: self.addlayer( H=H, vs=vs, vp=None, rho=None, zmin=z0)
            elif vpind==None and rhoind !=None:
                rho=rhoArr[i]; self.addlayer( H=H, vs=vs, vp=None, rho=rho, zmin=z0)
            elif vpind!=None and rhoind ==None:
                vp=vpArr[i]; self.addlayer( H=H, vs=vs, vp=vp, rho=None, zmin=z0)
            elif vpind!=None and rhoind !=None:
                vp=vpArr[i]; rho=rhoArr[i]; self.addlayer( H=H, vs=vs, vp=vp, rho=rho, zmin=z0)
        return
    
    def perturb(self, dm, zmin=0, zmax=9999, datatype='vs'):
        """ Add perturbation to the model given a depth range
        =======================================================
        Input Parameters:
        dm          - perturbed value (-1, 1)
        zmin, zmax  - depth range for the perturbation
        datatype    - data type for perturbation
        =======================================================
        """
        topArr = self.DepthArr - self.HArr
        bottomArr = self.DepthArr
        indexArr= np.arange(self.HArr.size)
        indexT =  topArr >=zmin
        indexB = bottomArr <=zmax
        if np.any(topArr == zmin):
            Ht =0.
        else:
            print 'Will add top layer !'
            Ht = ( topArr [topArr >zmin ]) [0] - zmin
            indext=(indexArr [indexT])[0] - 1
            z1 = zmin
        if np.any(bottomArr == zmax):
           Hb = 0.
        else:
            print 'Will add bottom layer !'
            if zmax < bottomArr [0]:
                Hb = zmax
                indexb = 0
                z2=0.
            else:
                Hb = zmax - ( bottomArr[bottomArr < zmax ]) [-1]
                indexb=(indexArr [indexB])[-1] + 1
                z2 = ( bottomArr [bottomArr < zmax ]) [-1]
        index = indexT*indexB
        if datatype=='vp':
            self.VpArr[index]=self.VpArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb]*(1.+dm), rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext]*(1.+dm), rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='vs':
            self.VsArr[index]=self.VsArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb]*(1.+dm), vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext]*(1.+dm), vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='rho':
            self.rhoArr[index]=self.rhoArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb]*(1.+dm),
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext]*(1.+dm),
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='qp':
            self.QpArr[index]=self.QpArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb]*(1.+dm), Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext]*(1.+dm), Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='qs':
            self.QsArr[index]=self.QsArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb]*(1.+dm), etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext]*(1.+dm), etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='etap':
            self.etapArr[index]=self.etapArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb]*(1.+dm), etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext]*(1.+dm), etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='etas':
            self.etasArr[index]=self.etasArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb]*(1.+dm),
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext]*(1.+dm),
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='frefp':
            self.frefpArr[index]=self.frefpArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb]*(1.+dm), frefs=self.frefsArr[indexb], zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext]*(1.+dm), frefs=self.frefsArr[indext], zmin=z1)
        if datatype=='frefs':
            self.frefsArr[index]=self.frefsArr[index]*(1.+dm)
            if Hb !=0:
                self.addlayer(H=Hb, vs=self.VsArr[indexb], vp=self.VpArr[indexb], rho=self.rhoArr[indexb],
                        Qp=self.QpArr[indexb], Qs=self.QsArr[indexb], etap=self.etapArr[indexb], etas=self.etasArr[indexb],
                            frefp=self.frefpArr[indexb], frefs=self.frefsArr[indexb]*(1.+dm), zmin=z2)
            if Ht !=0.:
                self.addlayer(H=Ht, vs=self.VsArr[indext], vp=self.VpArr[indext], rho=self.rhoArr[indext],
                        Qp=self.QpArr[indext], Qs=self.QsArr[indext], etap=self.etapArr[indext], etas=self.etasArr[indext],
                            frefp=self.frefpArr[indext], frefs=self.frefsArr[indext]*(1.+dm), zmin=z1)
        return
    
    def getArr4plot(self, zmin=-9999, zmax=9999, datatype='vs'):
        """
        Get data array for plotting purpose
        """
        index=(self.DepthArr<=zmax) * (self.DepthArr>zmin)
        depthArr=self.DepthArr[index]
        if datatype=='vp':
            dataArr=self.VpArr[index]
        if datatype=='vs':
            dataArr=self.VsArr[index]
        if datatype=='rho':
            dataArr=self.rhoArr[index]
        if datatype=='qp':
            dataArr=self.QpArr[index]
        if datatype=='qs':
            dataArr=self.QsArr[index]
        if datatype=='etap':
            dataArr=self.etapArr[index]
        if datatype=='etas':
            dataArr=self.etasArr[index]
        if datatype=='frefp':
            dataArr=self.frefpArr[index]
        if datatype=='frefs':
            dataArr=self.frefsArr[index]
        dataArr=np.repeat(dataArr, 2)
        depthArr=np.repeat(depthArr, 2)
        dataArr=np.append(dataArr, dataArr[-1])
        zmin=max(0., zmin)
        depthArr=np.append(0., depthArr)
        return dataArr, depthArr
    
    def plotvsak135(self, zmin=0, zmax=200., datatype='vs'):
        ak135model = Model1d()
        ak135model.ak135()
        dataArr, depthArr=self.getArr4plot(zmax=zmax, datatype=datatype)
        dataArr0, depthArr0=ak135model.getArr4plot(zmax=zmax, datatype=datatype)
        fig, ax= plt.subplots(figsize=(8,12))
        plt.plot(dataArr0, depthArr0, 'b-', lw=5, label='ak135')
        plt.plot(dataArr, depthArr, 'k--', lw=5, label='user model' )
        plt.ylabel('Depth (km)', fontsize=30)
        plt.xlabel(datatype+' (km/s)', fontsize=30)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        plt.yticks(np.append(0., self.DepthArr))
        plt.legend(loc='lower left', fontsize=25)
        plt.ylim((0, zmax))
        plt.gca().invert_yaxis()
        plt.show()
        return
    

    
class vprofile(object):
    def __init__(self, vsArr=np.array([]), vpArr=np.array([]), rhoArr=np.array([]), RmaxArr=np.array([]), RminArr=np.array([]),
                 z0Arr=np.array([]), HArr=np.array([]), xArr=np.array([]), yArr=np.array([]), dtypeArr=np.array([]), infname=None ):
        """
        An object to handle vertical profile, will be used by CPSPy
        ========================================================================================
        output txt format:
        vs/dvs    vp/dvp    rho/drho    Rmax    Rmin    z0    H    x    y    dtype
        
        dtype:
        0. absolute value for vs/vp/rho
        1. percentage value for dvs/dvp/drho 
        ========================================================================================
        """
        self.vsArr=vsArr
        self.vpArr=vpArr
        self.rhoArr=rhoArr
        self.RmaxArr=RmaxArr
        self.RminArr=RminArr
        self.z0Arr=z0Arr
        self.HArr=HArr
        self.xArr=xArr
        self.yArr=yArr
        self.dtypeArr=dtypeArr
        if infname!=None:
            self.read(infname=infname)
        return
    
    def read(self, infname):
        inArr=np.loadtxt(infname)
        try:
            self.vsArr=inArr[:,0]
            self.vpArr=inArr[:,1]
            self.rhoArr=inArr[:,2]
            self.RmaxArr=inArr[:,3]
            self.RminArr=inArr[:,4]
            self.z0Arr=inArr[:,5]
            self.HArr=inArr[:,6]
            self.xArr=inArr[:,7]
            self.yArr=inArr[:,8]
            self.dtypeArr=inArr[:,9]
        except:
            self.vsArr=np.array([inArr[0]])
            self.vpArr=np.array([inArr[1]])
            self.rhoArr=np.array([inArr[2]])
            self.RmaxArr=np.array([inArr[3]])
            self.RminArr=np.array([inArr[4]])
            self.z0Arr=np.array([inArr[5]])
            self.HArr=np.array([inArr[6]])
            self.xArr=np.array([inArr[7]])
            self.yArr=np.array([inArr[8]])
            self.dtypeArr=np.array([inArr[9]])
        return
    
    def write(self, outfname):
        N=self.vsArr.size
        outArr=np.append(self.vsArr, self.vpArr)
        outArr=np.append(outArr, self.rhoArr)
        outArr=np.append(outArr, self.RmaxArr)
        outArr=np.append(outArr, self.RminArr)
        outArr=np.append(outArr, self.z0Arr)
        outArr=np.append(outArr, self.HArr)
        outArr=np.append(outArr, self.xArr)
        outArr=np.append(outArr, self.yArr)
        outArr=np.append(outArr, self.dtypeArr)
        outArr=outArr.reshape(10, N)
        outArr=outArr.T
        self.check()
        np.savetxt(outfname, outArr, fmt='%g', header='vs/dvs vp/dvp rho/drho Rmax Rmin z0 H x y dtype(0: m, 1: dm)')
        return
    
    def check(self):
        N = self.vsArr.size
        for i in xrange(N):
            x=self.xArr[i]
            y=self.yArr[i]
            Rmax=self.RmaxArr[i]
            H=self.HArr[i]
            index = (self.xArr==x)* (self.yArr==y) * (self.RmaxArr==Rmax) * (self.HArr==H)
            if np.any(index):
                warnings.warn('Profile of same geometry and location exists.', UserWarning, stacklevel=1)
                return
        return
    
    def add(self, Rmax, Rmin, z0, H, x, y, dtype, vs=None, vp=None, rho=None):
        """
        Two kinds of input:
        1. Rmax or x is an array, the function will assume N vertical profile input(N=Rmax.size or x.size)
        2. If not, the function will first check if there is any previous profile with same x, y, Rmax. If yes, change
            vs, vp or rho in the existing profile; if not, add a new profile
        """
        if isinstance(x, np.ndarray) or isinstance(Rmax, np.ndarray):
            if vs ==None or vp == None or rho == None:
                raise ValueError('For array input, all model parameters must be assigned !')
            try:
                N=x.size
            except:
                N=Rmax.size
            ##################################################
            # Making all the input variables to be numpy array of size N
            ##################################################
            if not isinstance(Rmax, np.ndarray):
                Rmax=np.ones(N)*Rmax
            if not isinstance(Rmin, np.ndarray):
                Rmin=np.ones(N)*Rmin
            if not isinstance(z0, np.ndarray):
                z0=np.ones(N)*z0
            if not isinstance(H, np.ndarray):
                H=np.ones(N)*H
            if not isinstance(x, np.ndarray):
                x=np.ones(N)*x
            if not isinstance(y, np.ndarray):
                y=np.ones(N)*y
            if not isinstance(dtype, np.ndarray):
                dtype=np.ones(N)*dtype
            if not isinstance(vs, np.ndarray):
                vs=np.ones(N)*vs
            if not isinstance(vp, np.ndarray):
                vp=np.ones(N)*vp
            if not isinstance(rho, np.ndarray):
                rho=np.ones(N)*rho
        else:
            index = (self.xArr==x)* (self.yArr==y) * (self.RmaxArr==Rmax) * (self.HArr==H)
            # change model parameters if profile with same geometry/location exist
            if np.any(index):
                if vs !=None:
                    self.vsArr[index]=vs
                if vp !=None:
                    self.vpArr[index]=vp
                if rho !=None:
                    self.rhoArr[index]=rho
                return
            if dtype==0 and (vs ==None or vp == None or rho == None):
                raise ValueError('For absolute value profile, all model parameters must be assigned !')
            if vs==None:
                vs=0.
            if vp == None:
                vp=0.
            if rho == None:
                rho=0.
        self.vsArr=np.append( self.vsArr, vs)
        self.vpArr=np.append( self.vpArr, vp)
        self.rhoArr=np.append( self.rhoArr, rho)
        self.xArr=np.append( self.xArr, x)
        self.yArr=np.append( self.yArr, y)
        self.RmaxArr=np.append( self.RmaxArr, Rmax)
        self.RminArr=np.append( self.RminArr, Rmin)
        self.HArr=np.append( self.HArr, H)
        self.z0Arr=np.append( self.z0Arr, z0)
        self.dtypeArr=np.append( self.dtypeArr, dtype)
        return
    

        
    