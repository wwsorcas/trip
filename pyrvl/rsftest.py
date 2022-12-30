from vcl import Vector
from rsfvc import Space as rsfSpace
import data

# bulk modulus with lens
data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)

# space
sp = rsfSpace('m.rsf')
print('\n**** rsf space')
sp.myNameIs()

# vector
print('\n**** rsf vector')
m = Vector(sp)
m.myNameIs()

# link
print('\n*** link to m.rsf')
m.link('m.rsf')
m.myNameIs()

# dup
print('\n**** duplicate m.rsf')
m1 = m.dup()
m1.myNameIs()
