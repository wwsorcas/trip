import scenv
import data
import linalg

import sys

NTU = int(sys.argv[1])
NTD = int(sys.argv[2])
DELRTW = float(sys.argv[3])
DT  = float(sys.argv[4])
F1  = float(sys.argv[5])
F2  = float(sys.argv[6])
F3  = float(sys.argv[7])
F4  = float(sys.argv[8])
NS  = int(sys.argv[9])
DSX = float(sys.argv[10])

u = sys.argv[11]
w = sys.argv[12]
h = sys.argv[13]

### 2 shot headers
# adaptive filter prototype
data.rechdr(file=u,nt=NTU,dt=DT,
            ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=DELRTW,
            nshot=NS, dsx=DSX)

wsc=1.0e+3
data.bpfiltgather(file=w,nt=NTU,dt=DT,s=1.0,
            f1=F1,f2=F2,f3=F3,f4=F4,
            ntr=NS,sxstart=4200,szstart=3000,dsx=DSX,dsz=0)
linalg.scale(w,wsc)

# create zero data files with same source position, rz=500, rx=[2000,6000]
data.rechdr(file=h,nt=NTD,dt=DT, ntr=201,
            rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0,delrt=0,
            nshot=NS, dsx=DSX)



