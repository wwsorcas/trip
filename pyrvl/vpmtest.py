import vcl
import npvc
import numpy as np
import vpm

######################

n=4
sp0 = npvc.Space(n)
sp1 = npvc.Space(n)
sp2 = npvc.Space(2*n)

sp = vcl.ProductSpace([sp0, sp1])

f = vpm.expl1(sp0)

x = vcl.Vector(f.getDomain())

print('value at origin = ' + str(f(x)))

x.data[0] = 1.0
x.data[1] = 0.2
x.data[2] =-0.5
x.data[3] =-0.4

print('value at [1.0,0.2,-0.5,-0.4] = ' + str(f(x)))
print('gradient at [1.0,0.2,-0.5,-0.4] = ')
print(f.gradient(x).data)

x.data[0] = 1.0
x.data[1] = 0.0
x.data[2] = 0.0
x.data[3] = 0.0

fx = vcl.StandardJet(x, fcn=f)

print('value at [1.0,0.0,0.0,0.0] = ' + str(fx.value()))
print('gradient at [1.0,0.0,0.0,0.0] = ')
print(fx.gradient().data)

for i in range(11):
    x.data[0]=i*0.1
    fx = vcl.StandardJet(x, fcn=f)
    print('value at x[0]=' + str(x.data[0]) + ' = ' + str(fx.value()))

#############

print('\nvpmcgjet evaluation:')

g = vpm.sepexpl1(sp,sp2)

kmax = 10
eps  = 1.e-3
rho  = 1.e-3
b = vcl.Vector(g.getRange())
# nonzeros in b only in top half
for i in range(n):
    b.data[i]=(-i-1)**(i+1)

gx = vpm.vpmcgjet(x, g, b, kmax, eps, rho, verbose=2)

print('\nvpmcgjet value = ' + str(gx.value()))
print('\ngradient = ')
print(gx.gradient().data)
