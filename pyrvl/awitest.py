import linalg
import data
import vcl
import vcalg
import segyvc
import rsfvc
import asg
import os
import awi

    ################## domain, range spaces ###################

data.rechdr(file='baru.su',nt=251,dt=8.0,\
            ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=-1000)

usp=segyvc.Space('baru.su')
baru=vcl.Vector(usp,'baru.su')

# bulk modulus with less focussing lens
data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)

# homogeneous bulk modulus 
data.model(bulkfile='m0.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=1.0)

# bandpass filter source at sx=4200, sz=3000 (single trace)
data.bpfilt(file='wstar.su',nt=251,dt=8.0,s=1.0,f1=1.0,f2=2.5,f3=7.5,f4=12,sx=4200,sz=3000)
linalg.scale('wstar.su',1.0e+3)

# create zero data file with same source position, rz=500, rx=[2000,6000]
data.rechdr(file='g.su',nt=626,dt=8.0,ntr=201,
                rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0)

bulksp = rsfvc.Space('m.rsf')
datasp = segyvc.Space('g.su')

F = asg.fsbop(dom=bulksp, rng=datasp, \
            buoyancy='bym.rsf', source_p='wstar.su', \
            order=2, sampord=1, nsnaps=20,\
            cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
            nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

print('compute lens data')
m = vcl.Vector(bulksp,'m.rsf')
Fm = F(m)
print('compute hom data')
m0 = vcl.Vector(bulksp,'m0.rsf')
Fm0 = F(m0)

print('construct AWI penalty operator')
op = awi.awiop(dom=usp, d=Fm, p=Fm, precond=1, alpha=1.0e-4, kmax=50, eps=0.01, rho=0.01, verbose=2)
#op.myNameIs()

rhs = vcl.Vector(op.getRange())
rhs[0].copy(Fm)
e = vcl.Vector(op.getRange())

vcalg.conjgrad(x=baru, b=rhs, A=op, kmax=50, eps=0.01, rho=0.01, verbose=2)



#x1 = vcl.transp(op)*rhs
#y1 = op*x1
#x2 = vcl.transp(op)*y1
#y2 = op*x2
#x3 = vcl.transp(op)*y2

#print('x1 dot op^Ty2 = ' + str(x1.dot(x3)))
#print('op x1 dot y2 = ' + str(y1.dot(y2)))
