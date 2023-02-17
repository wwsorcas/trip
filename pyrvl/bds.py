import vcl
import npvc
import vcalg
import numpy as np

sp = npvc.Space(4)
x = vcl.Vector(sp)
u = vcl.Vector(sp)
l = vcl.Vector(sp)

u.data=2.0*np.ones((4,1))
l.data=-2.0*np.ones((4,1))

x.data[0]=-1.2
x.data[1]=1.0
x.data[2]=-1.2
x.data[3]=1.0

print(npvc.testbounds(u,l,x))

x.data[2]=3.0
                   
print(npvc.testbounds(u,l,x))

x.data[2]=-1.2

ful = npvc.ulbounds(sp,u,l)
iul = npvc.invulbounds(sp,u,l)

xx = iul(x)

yy = ful(xx)

print(xx.data)

print(yy.data)

FB = npvc.DoubleRosieWithBounds(sp,u,l)

print(FB(x).data)

b = vcl.Vector(sp)

b.data[0]=0.0
b.data[1]=-1.0
b.data[2]=0.0
b.data[3]=-1.0

FBC = vcl.comp(FB,ful)

vcalg.trgn(xx, b, FBC, imax=40, eps=0.001, kmax=10, rho=1.e-6, \
           Delta=10.0, mured=0.5, muinc=1.8, \
           gammared=0.1, gammainc=0.95, \
           gnverbose=1, cgverbose=0)

x = ful(xx)

print('\nsolution estimate:')
x.myNameIs()  
