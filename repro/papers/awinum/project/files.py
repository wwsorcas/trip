import scenv
import data
import linalg

NTU = [251,501]
NTD = [626,1251]
NX  = [401, 801]
NZ  = [201, 401]
DELRTW = [-1000, -500]
DT  = [8, 4]
DX  = [20, 10]
F1  = [1.0, 2.0]
F2  = [2.5, 5.0]
F3  = [7.5, 15.0]
F4  = [12.5, 25.0]

# adaptive filter prototype
data.rechdr(file='ul.su',nt=NTU[0],dt=DT[0],
            ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=DELRTW[0])
data.rechdr(file='uh.su',nt=NTU[1],dt=DT[1],
            ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=DELRTW[1])

wsc=1.0e+3
data.bpfilt(file='wl.su',nt=NTU[0],dt=DT[0],s=1.0,
            f1=F1[0],f2=F2[0],f3=F3[0],f4=F4[0],
            sx=4200,sz=3000)
data.bpfilt(file='wh.su',nt=NTU[1],dt=DT[1],s=1.0,
            f1=F1[1],f2=F2[1],f3=F3[1],f4=F4[1],
            sx=4200,sz=3000)
linalg.scale('wl.su',wsc)
linalg.scale('wh.su',wsc)

# create zero data files with same source position, rz=500, rx=[2000,6000]
data.rechdr(file='hl.su',nt=NTD[0],dt=DT[0],ntr=201,
            rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0,delrt=0)
data.rechdr(file='hh.su',nt=NTD[1],dt=DT[1],ntr=201,
            rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0,delrt=0)

# bulk modulus with smaller focussing lens
data.model(bulkfile='m0l.rsf', bulk=4.0, nx=NX[0], nz=NZ[0],
   dx=DX[0], dz=DX[0], lensfac=0.7)
data.model(bulkfile='m0h.rsf', bulk=4.0, nx=NX[1], nz=NZ[1],
   dx=DX[1], dz=DX[1], lensfac=0.7)

# bulk modulus with larger focussing lens
data.model(bulkfile='m1l.rsf', bulk=4.0, nx=NX[0], nz=NZ[0],
   dx=DX[0], dz=DX[0], lensfac=0.5, lensradd=0.4)
data.model(bulkfile='m1h.rsf', bulk=4.0, nx=NX[1], nz=NZ[1],
   dx=DX[1], dz=DX[1], lensfac=0.5, lensradd=0.4)

# homogeneous bulk modulus
data.model(bulkfile='mlhom.rsf', bulk=4.0, nx=NX[0], nz=NZ[0],
   dx=DX[0], dz=DX[0], lensfac=1.0)
data.model(bulkfile='mhhom.rsf', bulk=4.0, nx=NX[1], nz=NZ[1],
   dx=DX[1], dz=DX[1], lensfac=1.0)




