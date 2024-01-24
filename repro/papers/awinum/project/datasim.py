import scenv
import os
import vcl
import segyvc
import rsfvc
import asg

import sys
ORDER   = int(sys.argv[1])
SAMPORD = int(sys.argv[2])
NSNAPS  = int(sys.argv[3])
CFL     = float(sys.argv[4])
CMIN    = float(sys.argv[5])
CMAX    = float(sys.argv[6])
DMIN    = float(sys.argv[7])
DMAX    = float(sys.argv[8])
NL1     = int(sys.argv[9])
NL2     = int(sys.argv[10])
NR1     = int(sys.argv[11])
NR2     = int(sys.argv[12])
PMLAMPL = float(sys.argv[13])

dfile = sys.argv[14]
wfile = sys.argv[15]
mfile = sys.argv[16]
bfile = sys.argv[17]

bulksp = rsfvc.Space(mfile)
datasp = segyvc.Space(dfile)

m = vcl.Vector(bulksp, mfile)
d = vcl.Vector(datasp, dfile)

F = asg.fsbop(dom=bulksp, rng=datasp, \
        buoyancy=bfile, source_p=wfile, \
        order=ORDER, sampord=SAMPORD, nsnaps=NSNAPS,\
        cfl=CFL, cmin=CMIN, cmax=CMAX,dmin=DMIN, dmax=DMAX,\
        nl1=NL1, nr1=NR1, nl2=NL2, nr2=NR2, pmlampl=PMLAMPL)

d.copy(F(m))





