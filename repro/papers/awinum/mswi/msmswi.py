import os
import sys
import scenv

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

    def __init__(self, dom, rect1=10, rect2=10, repeat=2):
        RSFROOT = os.getenv('RSFROOT')
        smoother = os.path.join(RSFROOT,'bin/sfsmooth')
        self.dom = dom
        if not isinstance(self.dom,rsfvc.Space):
            raise Exception('grid smoother: provided domain not rsfvc.Space')
        self.cmd = smoother + ' rect1=' + str(rect1) + ' rect2=' + str(rect2) + ' repeat=' +str(repeat)
        
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

# arguments:
# output bulk mod, output filter,
# input bulk mod, input data, buoyancy, source wavelet
# fdpars
# smpars
# mswipars
# descpars
# lspars

def mswiopt(bulkmodout, filterout,
                bulkmodin, datain, filterin, buoyancy, source, 
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,
                boundstest=True, noisy=False,
                rect1=10, rect2=10, repeat=2,
                alpha=0.0001, sigma=0.00001, rho=0.0025, kmax=1000, verbose=0,
                descmax=12, desceps=0.01, descverbose=1, descout=None,
                lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9, lsverbose=1
                ):
    try:
        
        bulksp = rsfvc.Space(bulkmodin)
        datasp = segyvc.Space(datain)
        usp = segyvc.Space(filterin)
        
        m  = vcl.Vector(bulksp,bulkmodout)
        f = vcl.Vector(usp,filterout)
        m0 = vcl.Vector(bulksp,bulkmodin)
        d = vcl.Vector(datasp,datain)

        Winv = smoother(bulksp,rect1,rect2,repeat)
        
        F = asg.fsbop(bulksp, datasp, buoyancy, source, \
                    order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                    nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy \
                    )
        
        mswiargs = dict(dom=usp, sim=F, data=d, alpha=alpha, sigma=sigma, kmax=kmax, rho=rho, verbose=verbose)
        archargs = dict(archivepath=None, point='rsf', grad='rsf', filter='su', dataerr='su')
        sdargs = dict(Winv=Winv)
        lsargs = dict(lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9)
        
        Jx = vcalg.lsopt(m0, awi.mswi, SD=vcalg.SDwgrad, LS=vcalg.btls, descmax=descmax, desceps=desceps,\
                             descverbose=descverbose, descout=descout,\
                             lsverbose=lsverbose, lsargs=lsargs, jetargs=mswiargs, ddargs=sdargs, archargs=archargs)
        m.copy(Jx.point())
        f.copy(Jx.filter())
        
    except Exception as ex:
        print(ex)

if __name__ == '__main__':
    args = ", ".join(sys.argv[1:])
    cmd = 'mswiopt(' + args + ')'
    print(cmd)
    exec(cmd)



        

