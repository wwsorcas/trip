import vcl
import npvc
import vcalg

sp = npvc.Space(4)
x = vcl.Vector(sp)
b = vcl.Vector(sp)

x.data[0]=-1.2
x.data[1]=1.0
x.data[2]=-1.2
x.data[3]=1.0
b.data[0]=0.0
b.data[1]=-1.0
b.data[2]=0.0
b.data[3]=-1.0

F = npvc.DoubleRosie(sp)

vcalg.trgn(x, b, F, imax=40, eps=0.001, kmax=10, rho=1.e-6, Delta=10.0, \
            mured=0.5, muinc=1.8, \
            gammared=0.1, gammainc=0.95, gnverbose=0, cgverbose=0)

print('\nsolution estimate:')
x.myNameIs()
