import data
import vcl
import rsfvc
import segyvc
import asg
import linalg
import os
# bulk modulus with lens
data.model(bulkfile='mstar.rsf', bulk=4.0, nx=401, nz=201,\
               dx=20, dz=20, lensfac=0.7)
# bulk modulus without lens
data.model(bulkfile='m0.rsf', bulk=4.0, nx=401, nz=201,\
               dx=20, dz=20, lensfac=1.0)
# non-transient storage for gradient
data.model(bulkfile='g.rsf', bulk=4.0, nx=401, nz=201,\
               dx=20, dz=20, lensfac=1.0)
# bandpass filter source at sx=4200, sz=3000 (single trace)
data.bpfiltgather(file='wstar.su',nt=251,dt=8.0,s=1.0,\
                    f1=1.0,f2=2.5,f3=7.5,f4=12,\
                    ntr=1,sxstart=4200,szstart=3000,dsx=1000,dsz=0)
# create zero data file with same source position, rz=500, rx=[2000,6000]
data.rechdr(file='data.su',nt=626,dt=8.0,rxmin=2000,rxmax=6000,\
                ntr=201,rz=1000,sx=4200,sz=3000,delrt=0,nshot=1,dsx=1000)
# domain and range spaces
bulksp = rsfvc.Space('mstar.rsf')
datasp = segyvc.Space('data.su')
# wrap bulk modulus in Vector
mstar = vcl.Vector(bulksp,'mstar.rsf')
# instantiate modeling operator
F = asg.fsbop(dom=bulksp, rng=datasp, \
                  buoyancy='bymstar.rsf', source_p='wstar.su', \
                  order=2, sampord=1, nsnaps=20,\
                  cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                  nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)
# evaluate F[mstar], create noise-free data
dstar = F(mstar)
data = vcl.Vector(datasp,'data.su')
data.copy(dstar)
#dstar.myNameIs()
# least-squares function
JFWI = vcl.LeastSquares(F,dstar)
#JFWI.myNameIs()
# wrap bulk modulus in Vector
m0 = vcl.Vector(bulksp,'m0.rsf')
# convex combination of mstar, m0
t=0.7
m = vcl.Vector(bulksp)
m.copy(m0)
m.linComb(t,mstar,1-t)
# wrap gradient file in Vector
g = vcl.Vector(bulksp,'g.rsf')
# compute gradient, copy to non-transient storage
gtmp=JFWI.gradient(m)
g.copy(gtmp)
linalg.simplot(gtmp.data,addcb=True,clip=1e-12)
