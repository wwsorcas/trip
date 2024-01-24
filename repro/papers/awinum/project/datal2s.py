import scenv
import os
import vcl
import segyvc
import rsfvc
import asg

bulkspl = rsfvc.Space('m0l.rsf')
dataspl2s = segyvc.Space('hl2s.su')

m0l = vcl.Vector(bulkspl,'m0l.rsf')
m1l = vcl.Vector(bulkspl,'m1l.rsf')
os.system('/bin/cp hl2s.su d0l2s.su')
os.system('/bin/cp hl2s.su d1l2s.su')
d0l2s = vcl.Vector(dataspl2s, 'd0l2s.su')
d1l2s = vcl.Vector(dataspl2s, 'd1l2s.su')

# homog model
mlhom = vcl.Vector(bulkspl,'mlhom.rsf')
os.system('/bin/cp hl2s.su dlhom2s.su')
dlhom2s  = vcl.Vector(dataspl2s, 'dlhom2s.su')

F = asg.fsbop(dom=bulkspl, rng=dataspl2s, \
        buoyancy='bym0l.rsf', source_p='wl2s.su', \
        order=2, sampord=1, nsnaps=20,\
        cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
        nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

d0l2s.copy(F(m0l))
d1l2s.copy(F(m1l))
dlhom2s.copy(F(mlhom))





