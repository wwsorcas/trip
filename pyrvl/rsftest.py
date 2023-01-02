from vcl import Vector
from rsfvc import Space as rsfSpace
import data
import op
import os

# bulk modulus with lens
data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)

# bulk modulus with different grid
data.model(bulkfile='m1.rsf', bulk=4.0, nx=301, nz=151, dx=20, dz=20, lensfac=0.7)

def runtests():

    # space
    sp = rsfSpace('m.rsf')
    print('\n**** rsf space')
    #sp.myNameIs()

    # vector
    print('\n**** rsf vector')
    m = Vector(sp,'m.rsf')
    m.myNameIs()

    # link
    print('\n*** link to m.rsf')
    m.link('m.rsf')
    #m.myNameIs()

    # dup
    print('\n**** duplicate m.rsf')
    mdup = m.dup()
    #mdup.myNameIs()

    print("\n**** correct copy")
    try:
        mdup.copy(m)
    except Exception as ex:
        print(ex)
        
    # different space - same geom
    sp0 = rsfSpace('m.rsf')
    
    # vector in different space
    print('\n*** vector in different space')
    m0 = Vector(sp0)
    #m0.myNameIs()
    
    # incorrect copy
    print('*** incorrect copy')
    try:
        m0.copy(m)
    except Exception as ex:
        print(ex)
        
    # incorrect dot
    print('*** incorrect dot')
    try:
        d=m0.dot(m)
    except Exception as ex:
        print(ex)
        
    #del m0

    # yet another space
    sp1 = rsfSpace('m1.rsf')
    m1 = Vector(sp1)
    
    # incorrect lincomb
    print('*** incorrect lincomb')
    try:
        m1.linComb(1.0,m)
    except Exception as ex:
        print(ex)
        
    #del m1

    print('*** correct lincomb')
    try:
        mdup.linComb(1.0,m)
    except Exception as ex:
        print(ex)
        
    #del mdup
    
#print('*** clean all rsf and su files from current directory')
#op.cleanup()

runtests()






