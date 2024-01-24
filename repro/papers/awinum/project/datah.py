import scenv
import os
import vcl
import segyvc
import rsfvc
import asg

bulksph = rsfvc.Space('m0h.rsf')
datasph = segyvc.Space('hh.su')

m0h = vcl.Vector(bulksp,'mhl.rsf')
os.system('/bin/cp hh.su d0h.su')
d0h = vcl.Vector(datasp, 'd0h.su')

# homog model
mhhom = vcl.Vector(bulksp,'mhhom.rsf')
os.system('/bin/cp hh.su dhhom.su')
dhhom  = vcl.Vector(datasp, 'dhhom.su')

F = asg.fsbop(dom=bulksph, rng=datasph, \
        buoyancy='bym0h.rsf', source_p='wh.su', \
        order=2, sampord=1, nsnaps=20,\
        cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
        nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

d0h.copy(F(m0h))
dhhom.copy(F(mhhom))

