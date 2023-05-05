import vcl
import numpy as np
import npvc
import vcalg

mat = np.zeros((20,20))
for i in range(20):
    mat[i,i]=i+1

sp = npvc.Space(20)

A = npvc.MatrixOperator(sp,sp,mat)

xarr = np.ones((20,1))
x0 = vcl.Vector(sp,xarr)

b = A*x0

x = vcl.Vector(sp)
r = vcl.Vector(sp)
kmax = 30
rho = 1.e-3
verbose = 2

vcalg.bcg(x, b, A, kmax, rho, verbose, r)

print('\nDelta = 4.4')

Delta = 4.45
vcalg.bcg(x, b, A, kmax, rho, verbose, r, Delta)

print('x = ')
print(x.data)
print('r = ')
print(r.data)
print('Ax = ')
print((A*x).data)








