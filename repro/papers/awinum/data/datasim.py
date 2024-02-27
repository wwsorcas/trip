import scenv
import os
import vcl
import segyvc
import rsfvc
import asg

import sys

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

if __name__ == '__main__':
    args = ", ".join(sys.argv[1:])
    cmd = 'asgsim(' + args + ')'
    exec(cmd)


