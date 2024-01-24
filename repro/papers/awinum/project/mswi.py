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

ORDER=2
SAMPORD=1
NSNAPS=20
CFL=0.5
CMIN=1.0
CMAX=3.0
DMIN=0.8
DMAX=3.0
NL1=250
NL2=250
NR1=250
NR2=250
PMLAMPL=1.0

KMAX=1000
VERBOSE=2

def mswigrad(bulk, buoy, wavelet, data, filt, grad. alpha, sigma, rho):
    try: 

        bulksp = rsfvc.Space(bulk)
        datasp = segyvc.Space(data)
        filtsp = segyvc.Space(filt)

        m = vcl.Vector(bulksp, bulk)
        d = vcl.Vector(datasp, data)
        g = vcl.Vector(bulksp, grad)
        u = vcl.Vector(filtsp, filt)
        
        F = asg.fsbop(dom=bulksp, rng=datasp, \
                    buoyancy=buoy, source_p=wavelet, \
                    order=ORDER, sampord=SAMPORD, nsnaps=NSNAPS,\
                    cfl=CFL, cmin=CMIN, cmax=CMAX, dmin=DMIN, dmax=DMAX,\
                    nl1=NL1, nr1=NR1, nl2=NL2, nr2=NR2, pmlampl=PMLAMPL)

        j = awi.mswi(dom=filtsp, sim=F, mod=m, data=d, alpha=alpha,
                         sigma=sigma, kmax=KMAX, rho=rho, verbose=VERBOSE)

        u.copy(j.filter())
        g.copy(j.gradient())

    except Exception as ex:
        print(ex)
        raise Exception('called from mswigrad')

try:
        
    # argv1 = rsf file defining model space and bulk mod
    # argv2 = rsf file defining buoyancy
    # argv3 = su file defining wavelet
    # argv4 = su file defining data space and data
    # argv5 = su file defining adaptive filter space
    # argv6 = adaptive filter (out)
    # argv7 = rsf file defining gradient (out)
    # argv8 = alpha
    # argv9 = sigma
    # argv10 = rho

    os.system('/bin/cp ' + sys.argv[5] + ' ' + sys.argv[6])
    os.system('sfcp ' + sys.argv[1] + ' ' + sys.argv[7])

    mswigrad(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
        sys.argv[6],sys.argv[7], float(sys.argv[8]), float(sys.argv[9]),
        float(sys.argv[10]))
                 
