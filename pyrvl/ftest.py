import vcl
import npvc
from vcl import transp
from vcl import Vector
from npvc import Space as npSpace

dom = npSpace(2)
rng = npSpace(3)
f = npvc.OpExpl1(dom,rng)
x = Vector(dom)
y = Vector(rng)
x.data[0]=1
x.data[1]=-2
print('input vector:')
x.myNameIs()
f.apply(x,y)
print('output of apply method:')
y.myNameIs()
z = f(x)
print('output of overloaded eval')
z.myNameIs()
dfx = f.deriv(x)
print('output of deriv method:')
dfx.myNameIs()
dx = Vector(dom)
dy = Vector(rng)
dx.data[0]=1
dx.data[1]=2
print('input of deriv overloaded eval')
dx.myNameIs()
print('output of deriv overloaded eval = input of adjoint eval')
#dy=dfx(dx)
dy=dfx*dx
dy.myNameIs()
print('output of deriv adjoint overloaded eval')
#dz=transp(dfx)(dy)
dz=transp(dfx)*dy
dz.myNameIs()




