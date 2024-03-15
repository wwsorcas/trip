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

def fwiopt(bulkmodout, bulkmodin, datain, buoyancy, source, 
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,
                boundstest=True, noisy=False, partask=0,
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

        Winv = rsfvc.smoother(bulksp,rect1,rect2,repeat)
        
        F = asg.fsbop(bulksp, datasp, buoyancy, source, \
                    order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                    nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy, partask \
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



        

