import os
import vcl
import segyvc
import rsfvc
import asg

# order=2, sampord=1, nsnaps=20,\
#                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
#                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,
#                boundstest=True, noisy=False

def asgsim(bulkmod, data_p, buoyancy, source_p, order=2, sampord=1, nsnaps=20,\
               cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
               nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,
               boundstest=True, noisy=False, partask=0):
    
    try:
    
        bulksp = rsfvc.Space(bulkmod)
        datasp = segyvc.Space(data_p)

        m = vcl.Vector(bulksp, bulkmod)
        d = vcl.Vector(datasp, data_p)

        F = asg.fsbop(bulksp, datasp, buoyancy, source_p, \
                          order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                          nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy, partask \
                          )

        d.copy(F(m))
        
    except Exception as ex:
        print(ex)
        raise Exception('called from sim')


#########################################################################

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
                lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9,
                lsverbose=1):
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
        
        Jx = vcalg.lsopt(m0, vcl.LeastSquaresJet, SD=vcalg.SDwgrad, \
                             LS=vcalg.btls, descmax=descmax, desceps=desceps, \
                             descverbose=descverbose, descout=descout, \
                             lsverbose=lsverbose, lsargs=lsargs, \
                             jetargs=fwiargs, ddargs=sdargs, archargs=archargs)
        m.copy(Jx.point())
        
    except Exception as ex:
        print(ex)

#########################################################################

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

def mswiopt(bulkmodout, filterout,
                bulkmodin, datain, filterin, buoyancy, source, 
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0,\
                boundstest=True, noisy=False, partask=0,\
                rect1=10, rect2=10, repeat=2,\
                alpha=0.0001, sigma=0.00001, rho=0.0025, kmax=1000, verbose=0,\
                etar=None, ratminus=None, ratplus=None, upmax=1,\
                descmax=12, desceps=0.01, descverbose=1, descout=None,\
                lsmax=10, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9,\
                lsverbose=1
                ):
    try:
        
        bulksp = rsfvc.Space(bulkmodin)
        datasp = segyvc.Space(datain)
        usp = segyvc.Space(filterin)
        
        m  = vcl.Vector(bulksp,bulkmodout)
        f = vcl.Vector(usp,filterout)
        m0 = vcl.Vector(bulksp,bulkmodin)
        d = vcl.Vector(datasp,datain)

        Winv = rsfvc.smoother(bulksp,rect1,rect2,repeat)
        
        F = asg.fsbop(bulksp, datasp, buoyancy, source, \
                    order, sampord, nsnaps, cfl, cmin, cmax, dmin, dmax,\
                    nl1, nr1, nl2, nr2, pmlampl, boundstest, noisy, partask \
                    )
        
        mswiargs = dict(dom=usp, sim=F, data=d, alpha=alpha, sigma=sigma, \
                            kmax=kmax, rho=rho, verbose=verbose, \
                            etar=etar, ratminus=ratminus, ratplus=ratplus, upmax=upmax)
        archargs = dict(archivepath=None, point='rsf', grad='rsf', \
                            filter='su', dataerr='su')
        sdargs = dict(Winv=Winv)
        lsargs = dict(lsmax=lsmax, mured=mured, muinc=muinc, \
                          gammared=gammared, gammainc=gammainc)
        
        Jx = vcalg.lsopt(m0, awi.mswi, SD=vcalg.SDwgrad, LS=vcalg.btls,\
                             descmax=descmax, desceps=desceps,\
                             descverbose=descverbose, descout=descout,\
                             lsverbose=lsverbose, lsargs=lsargs, \
                             jetargs=mswiargs,\
                             ddargs=sdargs, archargs=archargs)
        m.copy(Jx.point())
        f.copy(Jx.filter())
        
    except Exception as ex:
        print(ex)

#########################################################################

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
                etar=None, ratminus=None, ratplus=None, upmax=1,\
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
                            etar=etar, ratminus=ratminus, ratplus=ratplus, upmax=upmax)
                            
        archargs = dict(archivepath=None, point='rsf', grad='rsf', filter='su', dataerr='su')
        sdargs = dict(Winv=Winv)
        lsargs = dict(lsmax=lsmax, mured=mured, muinc=muinc, \
                          gammared=gammared, gammainc=gammainc)
        
        Jx = vcalg.lsopt(cx0, awi.mswi, SD=vcalg.SDwgrad, LS=vcalg.btls,\
                             descmax=descmax, desceps=desceps,\
                             descverbose=descverbose, descout=descout,\
                             lsverbose=lsverbose, lsargs=lsargs, jetargs=mswiargs,\
                             ddargs=sdargs, archargs=archargs)
        m.copy(fwd(Jx.point()))
        f.copy(Jx.filter())
        
    except Exception as ex:
        print(ex)

#########################################################################

import numpy as np
import math
import data
import linalg
import sys

def cam(rad=800, inner=2.0, outer=4.0, n1=201, n2=401, d1=20, d2=20):
    try:

        # create rsf files
        data.rsffile('cambulk.rsf', 'Bulk_modulus', 'GPa', n1, n2, d1, d2, val=1.0)
        data.rsffile('cambuoy.rsf', 'Buoyancy', 'cc/g', n1, n2, d1, d2, val=1.0)

        # create array
        z = np.ones(n1*n2)
        cam = z.reshape(n2,n1)

        # loop through array scaling as appropriate
        for i in range(n1):
            z = i*20.0 - 2000.0
            for j in range(n2):
                x = j*20 - 4000.0
                d = math.sqrt(x*x + z*z)
                if d<rad:
                    cam[j,i] *= inner
                else:
                    cam[j,i] *= outer

        # write bulkmod to rsf file
        linalg.ndarraytorsfdata(cam,'cambulk.rsf')

    except Exception as ex:
        print(ex)
        raise Exception('called from cam')

#########################################################################


        






