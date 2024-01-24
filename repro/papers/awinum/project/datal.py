import scenv
import os
import vcl
import segyvc
import rsfvc
import asg

bulkspl = rsfvc.Space('m0l.rsf')
dataspl = segyvc.Space('hl.su')

m0l = vcl.Vector(bulkspl,'m0l.rsf')
m1l = vcl.Vector(bulkspl,'m1l.rsf')
os.system('/bin/cp hl.su d0l.su')
os.system('/bin/cp hl.su d1l.su')
d0l = vcl.Vector(dataspl, 'd0l.su')
d1l = vcl.Vector(dataspl, 'd1l.su')

# homog model
mlhom = vcl.Vector(bulkspl,'mlhom.rsf')
os.system('/bin/cp hl.su dlhom.su')
dlhom  = vcl.Vector(dataspl, 'dlhom.su')

F = asg.fsbop(dom=bulkspl, rng=dataspl, \
        buoyancy='bym0l.rsf', source_p='wl.su', \
        order=2, sampord=1, nsnaps=20,\
        cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
        nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

d0l.copy(F(m0l))
d1l.copy(F(m1l))
dlhom.copy(F(mlhom))





