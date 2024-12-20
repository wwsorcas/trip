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

# file for sim output
data.model(bulkfile='dm.rsf', bulk=4.0, nx=401, nz=201,\
               dx=20, dz=20, lensfac=0.7)
               
# bandpass filter source at sx=4200, sz=3000 (single trace)
data.bpfilt(file='wstar.su',nt=251,dt=8.0,s=1.0,\
                f1=1.0,f2=2.5,f3=7.5,f4=12,sx=4200,sz=3000)
                
# create zero data file with same source position, rz=500, rx=[2000,6000]
data.rechdr(file='data.su',nt=626,dt=8.0,rxmin=2000,rxmax=6000,\
                ntr=201,rz=1000,sx=4200,sz=3000)

# instantiate modeling operator
args = 'bulkmod=mstar.rsf buoyancy=bymstar.rsf ' + \
  'bulkmod_b1=dm.rsf buoyancy_b1=bydm.rsf ' + \
  'data_p=data.su source_p=wstar.su deriv=1 adjoint=1 ' +\
  'order=2 sampord=1 nsnaps=20 cfl=0.5 cmin=1.0 cmax=3.0 dmin=0.8 dmax=3.0 ' +\
  'nl1=250 nr1=250 nl2=250 nr2=250 pmlampl=1.0'

TRIP = os.getenv('TRIP')
cmd = os.path.join(TRIP,'iwave/asg/main/sim.x') + ' ' + args
print(cmd)

ret=os.system(cmd + ' >& jnk')
print('ret=' + str(ret))

