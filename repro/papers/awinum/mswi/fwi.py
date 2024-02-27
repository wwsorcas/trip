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

def fwiopt(bulkmodout, bulkmodin, datain, buoyancy, source, 
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,
                boundstest=True, noisy=False,
                rect1=10, rect2=10, repeat=2,
                descmax=12, desceps=0.01, descverbose=1, descout=None,
                lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9, lsverbose=1
                ):
    try:
        
        bulksp = rsfvc.Space(bulkmodin)
        datasp = segyvc.Space(datain)
        
        m  = vcl.Vector(bulksp,bulkmodout)
        m0 = vcl.Vector(bulksp,bulkmodin)
        d = vcl.Vector(datasp,datain)

        Winv = smoother(bulksp,rect1,rect2,repeat)
        
        F = asg.fsbop(bulksp, datasp, buoyancy, source, \
                    order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                    nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy \
                    )

        fwiargs = dict(F=F, b=d)
        archargs = dict(archivepath=None, point='rsf', grad='rsf', dataerr='su')
        sdargs = dict(Winv=Winv)
        lsargs = dict(lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9)
        
        Jx = vcalg.lsopt(m0, vcl.LeastSquaresJet, SD=vcalg.SDwgrad, LS=vcalg.btls, descmax=descmax, desceps=desceps, \
                             descverbose=descverbose, descout=descout, \
                             lsverbose=lsverbose, lsargs=lsargs, jetargs=fwiargs, ddargs=sdargs, archargs=archargs)
        m.copy(Jx.point())
        
    except Exception as ex:
        print(ex)

if __name__ == '__main__':
    args = ", ".join(sys.argv[1:])
    cmd = 'fwiopt(' + args + ')'
    exec(cmd)



        

