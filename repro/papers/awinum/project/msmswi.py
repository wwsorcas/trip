import os
import sys

import linalg
import data
import vcl
import vcalg
import vpm
import segyvc
import rsfvc
import asg
import awi
import tmp

class smoother(vcl.LinearOperator):

    def __init__(self, dom):
        RSFROOT = os.getenv('RSFROOT')
        smoother = os.path.join(RSFROOT,'bin/sfsmooth')
        self.cmd = smoother + ' rect1=10 rect2=10 repeat=2 '
        self.dom = dom
        if not isinstance(self.dom,rsfvc.Space):
            raise Exception('grid smoother: provided domain not rsfvc.Space')

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.dom

    def applyFwd(self,x,y):
        try:
            ret = os.system(self.cmd + ' < ' + x.data + ' > ' + y.data)
            if ret != 0:
                raise Exception('command failed: ' + self.cmd)
        except Exception as ex:
            print(ex)
            raise Exception('called from smoother.applyFwd')

    def applyAdj(self,x,y):
        try:
            ret = os.system(self.cmd + ' < ' + x.data + ' > ' + y.data)
            if ret != 0:
                raise Exception('command failed: ' + self.cmd)
        except Exception as ex:
            print(ex)
            raise Exception('called from smoother.applyAdj')

    def myNameIs():
        print('2D Grid Smoother using sfsmooth')
        print('rect1=10 rect2=10 repeat=2')
                  
try:

    # arg1 = job choice: 'vals', 'grad0', 'grad', or gradtest (default)
    # arg2 = frequency index
    # arg3 = lens index
    # arg4 = alpha
    # arg5 = sigma
    # arg6 = rho
    
    myalpha = 0.0001
    mysigma = 0.00001
    myrho = 0.0025

    wsc=1.0e+3
    mykmax=1000

    NTU = 251
    NTD = 626
    NX  = 401
    NZ  = 201
    DELRTW =-1000
    DT  = 8
    DX  = 20
    F1  = 1.0
    F2  = 2.5
    F3  = 7.5
    F4  = 12.5
    LENSRADD = 0.4
    LENSFAC = 0.5
    NTR = 201
    NSHOT = 1
    DRX = 20
    DSX = 200
    RX=2000
    SX=3000

    print('build spaces, models')
    
    data.rechdr(file='baru.su',nt=NTU,dt=DT,\
                    ntr=NTR,rx=RX,rz=1000,sx=SX,sz=3000,drx=DRX,delrt=-1000,
                    nshot=NSHOT, dsx=DSX)

    usp=segyvc.Space('baru.su')
    u0 = vcl.Vector(usp,'baru.su')

    # bulk modulus with less focussing lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=NX, nz=NZ,
                   dx=DX, dz=DX, lensfac=LENSFAC, lensradd=LENSRADD)

    # homogeneous bulk modulus 
    data.model(bulkfile='m0.rsf', bulk=4.0, nx=NX, nz=NZ,
                   dx=DX, dz=DX, lensfac=1.0)

    # bandpass filter source at sx=4200, sz=3000 (single trace)
    data.bpfiltgather(file='wstar.su',nt=NTU,dt=DT,s=1.0,
                          f1=F1,f2=F2,f3=F3,f4=F4,
                          ntr=NSHOT, sxstart=SX, szstart=3000,
                          dsx=DSX, dsz=0)
    linalg.scale('wstar.su',wsc)
    
    # create zero data file with same source position, rz=500, rx=[2000,6000]
    data.rechdr(file='d.su',nt=NTD,dt=DT,ntr=NTR,
                    rx=RX, rz=1000.0, sx=SX, sz=3000, drx=DRX,
                    nshot=NSHOT, dsx=DSX)
    
    bulksp = rsfvc.Space('m.rsf')
    datasp = segyvc.Space('d.su')

    # target model
    m  = vcl.Vector(bulksp,'m.rsf')
    # homog model
    m0 = vcl.Vector(bulksp,'m0.rsf')

    Winv = smoother(dom=bulksp)

    F = asg.fsbop(dom=bulksp, rng=datasp, \
                buoyancy='bym.rsf', source_p='wstar.su', \
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

    d = vcl.Vector(datasp,'d.su')
    d.copy(F(m))
    
    #    sfcp = os.path.join(RSFROOT,'bin/sfcp')

    # max number of descent steps
    nsteps=12
    # max number of line search steps
    nstepint=10
    # lower G-A limit
    minfac = 0.1
    # upper G-A limit
    maxfac = 0.9
    
    mswiargs = dict(dom=usp, sim=F, data=d, alpha=myalpha, sigma=mysigma, kmax=mykmax, rho=myrho, verbose=2)
    sdargs = dict(Winv=Winv)
    lsargs = dict(lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9, lsverbose=0)

    Jx = tmp.lsopt(m0, awi.mswi, DD=tmp.DDwgrad, descmax=12, desceps=0.01, descverbose=1, lsargs=lsargs, jetargs=mswiargs, ddargs=sdargs)  
    
except Exception as ex:
    print(ex)

        

