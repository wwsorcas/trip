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


# arguments:
# output bulk mod, output filter,
# input bulk mod, input data, buoyancy, source wavelet
# fdpars
# smpars
# mswipars
# descpars
# lspars

def modmswiopt(bulkmodout, filterout,
                bulkmodin, cmodin, datain, filterin, buoyancy, source, 
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,\
                boundstest=True, noisy=False, partask=0,\
                rect1=10, rect2=10, repeat=2,\
                alpha=0.0001, sigma=0.00001, rho=0.0025, kmax=1000, verbose=0,\
                etar=None, ratminus=None, ratplus=None,\
                descmax=12, desceps=0.01, descverbose=1, descout=None,\
                lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9,\
                lsverbose=1
                ):
    try:
        
        bulksp = rsfvc.Space(bulkmodin)
        cxsp    = rsfvc.Space(cmodin)
        datasp = segyvc.Space(datain)
        usp = segyvc.Space(filterin)
        
        m  = vcl.Vector(bulksp,bulkmodout)
        f = vcl.Vector(usp,filterout)
        m0 = vcl.Vector(bulksp,bulkmodin)
        d = vcl.Vector(datasp,datain)

        Winv = rsfvc.smoother(cxsp,rect1,rect2,repeat)

        F = asg.fsbop(bulksp, datasp, buoyancy, source, \
                    order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                    nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy, partask \
                    )

        # inverse map expanded c from bulkmod
        inv = asg.invrntobulkfb(bulksp, cxsp, buoyancy,
            cmin=cmin, cmax=cmax, dmin=dmin, dmax=dmax)
        # initial expanded velo model, corresponding to m
        cx0 = inv(m0)
        
        # forward map expanded c to bulkmod
        fwd = asg.rntobulkfb(cxsp, bulksp, buoyancy,
            cmin=cmin, cmax=cmax, dmin=dmin, dmax=dmax)

        # composition: modeling on expanded velo space
        Fmod = vcl.comp(F, fwd)

        mswiargs = dict(dom=usp, sim=Fmod, data=d, alpha=alpha, sigma=sigma, \
                            kmax=kmax, rho=rho, verbose=verbose, \
                            etar=etar, ratminus=ratminus, ratplus=ratplus)
                            
        archargs = dict(archivepath=None, point='rsf', grad='rsf', filter='su', dataerr='su')
        sdargs = dict(Winv=Winv)
        lsargs = dict(lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9)
        
        Jx = vcalg.lsopt(cx0, awi.mswi, SD=vcalg.SDwgrad, LS=vcalg.btls,\
                             descmax=descmax, desceps=desceps,\
                             descverbose=descverbose, descout=descout,\
                             lsverbose=lsverbose, lsargs=lsargs, jetargs=mswiargs,\
                             ddargs=sdargs, archargs=archargs)
        m.copy(fwd(Jx.point()))
        f.copy(Jx.filter())
        
    except Exception as ex:
        print(ex)

if __name__ == '__main__':
    args = ", ".join(sys.argv[1:])
    cmd = 'modmswiopt(' + args + ')'
    print(cmd)
    exec(cmd)



        

