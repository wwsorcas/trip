### development test file for vcl2.0, a "more python" version of vcl

import vcl2
import numpy as np
import npvc2
#import segyvc
import op
import linalg
from functools import partial

def resid(x,y,mat=None,b=None):
    try:
        if mat is None:
            raise Exception('mat = (matrix) required')
        if b is None:
            raise Exception('b = (vector) required')
        if y is None:
            print('matmult')
            print('residual y=mat*x - b')
            print('matrix = ')
            print(mat)
            print('b = ')
            print(b)
            return
        y=mat@x-b
    except Exception as ex:
        print(ex)
        raise Exception('called from resid')
    else:
        return y

def squarem(x,y):
    if y is None:
        print('squarem')
        print('square input vector element-wise')
        return
    y=x*x
    return y

def dsquarem(x,dx,dy):
    if dy is None:
        print('dsquarem')
        print('derivative of square input vector element-wise')
        return
    dy=2*x*dx
    return dy

def conv(input,output,kernel=None):
    try:
        if kernel is None:
            raise Exception('kernel argument defined as None')            
        if output is None:
            print('conv op rep')
            print('kernel = ' + kernel)
            return
        if not op.convop(kernel,input,output,adj=0):
            raise Exception('Error: op.convop call failed')
    except Exception as ex:
        print(ex)
        raise Exception('called from conv')
    else:
        return output       
    
try:
    sp = npvc2.Space(4)
    f = vcl2.StdFunction(sp,sp,squarem)

    x = vcl2.Vector(sp)
    for i in range(4):
        x.data[i]=i
    print('squarem:')
    y=f(x)
    print(y.data)
except Exception as ex:
    print(ex)

try:

    ### np expl
    mat = np.zeros((3,2))
    mat[0,0] = 1.0
    mat[0,1] = 0.0
    mat[1,0] = 3.0
    mat[1,1] = -2.0
    mat[2,0] = 0.0
    mat[2,1] = 2.0

    insp = npvc2.Space(2)
    outsp = npvc2.Space(3)
    mataux = {'mat' : mat}

    f = vcl2.StdFunction(insp,outsp,npvc2.matmult,aux=mataux)

    x = vcl2.Vector(insp)
    x.data[0]=1
    x.data[1]=-1

    y = f(x)

    print('\nmatmult:')
    print(y.data)

    f.myNameIs()

    print('\ntry to multiply 3x2 by 3x1')
    z = f(y)
except Exception as ex:
    print(ex)

try:

    A = vcl2.StdLinearOperator(insp,outsp,npvc2.matmult,npvc2.adjmult,aux=mataux)

    y=A*x

    print('\nlinear matmult')
    print(y.data)

    AT = vcl2.transp(A)
    z = AT*y
    
    print('\ntransp matmult')
    print(z.data)

    print('\ncompare with numpy')
    print(mat.T@y.data)

    print('\ntry to multiply adjoint of 3x2 by 2x1')
    w = AT*x
    
except Exception as ex:
    print(ex)

try:
    resaux = {
        'mat' : mat,
        'b' : y.data
    }
    f = vcl2.StdFunction(insp,outsp,resid,aux=resaux)
#    f = vcl2.StdFunction(insp,outsp,resid,mat=mat,b=y.data)

    z = f(x)

    print('\nresid:')
    print(z.data)
    print('A*x')
    print(mat@x.data)
    print('b')
    print(y.data)
    f.myNameIs()
except Exception as ex:
    print(ex)

print('\n#####################\n')

def printem(**aux):
    for arg in aux:
        print(arg + ' = ' + aux[arg])

printem(a = '1', b = '2')

stuff = {'a':'1', 'b':'2'}

printem(**stuff)

print('\n#####################\n')

def passem(**jnk):
    printem(**jnk)

passem(a = '1', b = '2')
passem(**stuff)

print('\n########## ex ###########\n')

def exprintem(a, b):
    print('a = ' + a)
    print('b = ' + b)

ina = '1'
inb = '2'

exprintem(a=ina, b=inb)

print('stuff')

exprintem(**stuff)

print('wrong order')
try:
    exprintem(b=inb, a=ina)
except Exception as ex:
    print(ex)

print('wrong arg')
try:
    exprintem(c=ina)
except Exception as ex:
    print(ex)

print('missing arg')
try:
    exprintem(b=ina)
except Exception as ex:
    print(ex)

ina = '0'
inb = '0'

print('no args')
try:
    exprintem()
except Exception as ex:
    print(ex)


class foo():

    def __init__(self,**jnk):
        self.jnk = jnk

    def doit(self):
        print('\n########## doit ###########\n')        
        printem(**self.jnk)

print('same inputs')
jnkobj = foo(**stuff)
jnkobj.doit()

print('no inputs')
emp = foo()
emp.doit()

