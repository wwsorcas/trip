### development test file for vcl2.0, a "more python" version of vcl

import vcl
import numpy as np
import npvc
import segyvc
import op
import linalg
from functools import partial

def matmult(x,y,mat=None):
    try:
        if mat is None:
            raise Exception('mat argument defined as None')
        if y is None:
            print('matmult')
            print('matrix multiplication function: matrix = ')
            print(mat)
            return
        y=mat@x
    except Exception as ex:
        print(ex)
        raise Exception('called from matmult')
    else:
        return y

def adjmult(x,y,mat=None):
    try:
        if mat is None:
            raise Exception('mat argument defined as None')
        if y is None:
            print('matmult')
            print('matrix multiplication function: matrix = ')
            print(mat)
            return
        y=mat.T@x
    except Exception as ex:
        print(ex)
        raise Exception('called from matmult')
    else:
        return y


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
    
# simple nonlinear op class
class StdFunction:
    __doc__ = '''
    domain and range are vcl.Space objects

    apply is a function with signature apply(x,y,aux), 
    where
        x is a data object for input vector in the domain
        y is a data object for output vector in the range
        aux is a dictionary of additional arguments (if any)

    applyFwd is a function with signature 
        applyFwd(x,dx,dy,aux)
    where
        x is a data object for input reference vector in the domain
        dx is a data object for input perturbation vector in the domain
        dy is a data object for output perturbation vector in the range
        aux is a dictionary listing auxiliary arguments for apply, passed via
            kwargs

    applyAdj is a function with signature 
        applyAdj(x,dy,dx,aux)
    where
        x is a data object for input reference vector in the domain
        dy is a data object for input perturbation vector in the range
        dx is a data object for output perturbation vector in the domain
        aux is a dictionary listing auxiliary arguments for apply, passed via
            kwargs

    '''

    def __init__(self,dom,rng,apply,applyFwd=None,applyAdj=None,aux=None):
        self.dom = dom
        self.rng = rng
        self.apply = apply
        self.applyFwd = applyFwd
        self.applyAdj = applyAdj
        self.aux = aux
        
    def __call__(self,x):
        try:
            print('call')
            if x.space != self.dom:
                raise Exception('Error: input vec not in domain')
            print('y')
            y = vcl.Vector(self.rng)
            if self.aux is None:
                print('data')
                y.data=self.apply(x.data,y.data)
            else:
                y.data=self.apply(x.data,y.data,**self.aux)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.StdFunction.operator()')
        else:
            return y

    def deriv(self,x):
        try:
            if self.applyFwd is None:
                raise Exception('deriv construction requires non-None applyFwd')
            if self.applyAdj is None:
                raise Exception('deriv construction requires non-None applyAdj')
            return LinearFunction(self.dom,self.rng,
                                      partial(self.applyFwd,x.data),
                                      partial(self.applyAdj,x.data),
                                      aux = self.aux)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.LinearFunction.deriv')

    def myNameIs(self):
        print('\nvcl.Function instance')
        print('*** domain')
        self.dom.myNameIs()
        print('*** range')
        self.rng.myNameIs()
        print('*** apply function')
        self.apply(None,None,**self.aux)
        if self.applyFwd is not None:
            print('*** deriv: applyFwd supplied')
        if self.applyAdj is not None:
            print('*** deriv: applyAdj supplied')

class LinearFunction(StdFunction):

    '''
    special case of Function in which action is presumed to be linear
    and adjoint is supplied

    applyFwd and applyAdj have same signature as apply in Function
    '''

    def __init__(self,dom,rng,applyFwd,applyAdj,aux=None):

        self.dom = dom
        self.rng = rng
        self.applyFwd = applyFwd
        self.applyAdj = applyAdj
        self.aux = aux    

    def __call__(self,x):
        return __mul__(x)

    def __mul__(self,x):
        try:
            if x.space != self.dom:
                raise Exception('Error: input vec not in domain')
            y = vcl.Vector(self.rng)
            if self.aux is None:
                y.data=self.applyFwd(x.data,y.data)
            else:
                y.data=self.applyFwd(x.data,y.data,**self.aux)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector operator()')
        else:
            return y

    def myNameIs(self):
        print('\nvcl.LinearFunction instance')
        print('*** domain')
        self.dom.myNameIs()
        print('*** range')
        self.rng.myNameIs()
        print('*** apply function')
        self.applyFwd(None,None,**self.aux)

class transp(LinearFunction):

    def __init__(self, lf):
        self.dom = lf.rng
        self.rng = lf.dom
        self.applyFwd = lf.applyAdj
        self.applyAdj = lf.applyFwd
        self.aux = lf.aux

    def myNameIs(self):
        print('\nTranspose of Linear Function:')
        self.lf.myNameIs()

try:
    sp = npvc.Space(4)
    f = StdFunction(sp,sp,squarem)

    x = vcl.Vector(sp)
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

    insp = npvc.Space(2)
    outsp = npvc.Space(3)
    mataux = {'mat' : mat}

    f = StdFunction(insp,outsp,matmult,aux=mataux)

    x = vcl.Vector(insp)
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

    A = LinearFunction(insp,outsp,matmult,adjmult,aux=mataux)

    y=A*x

    print('\nlinear matmult')
    print(y.data)

    AT = transp(A)
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
    f = StdFunction(insp,outsp,resid,aux=resaux)

    z = f(x)

    print('\nresid:')
    print(z.data)
    f.myNameIs()
except Exception as ex:
    print(ex)


try:
    ### segy expl
    domsp = segyvc.Space('wstar.su')
    input = vcl.Vector(domsp,'wstar.su')

    rngsp = segyvc.Space('barw.su')
    output = vcl.Vector(rngsp)

    convaux = {'kernel' : 'baru.su'}

    g = StdFunction(domsp,rngsp,conv,aux=convaux)

    output = g(input)

    linalg.copy(output.data,'savedoutput.su')

    print('\nconvolutional operator')
    g.myNameIs()
    
except Exception as ex:
    print(ex)
    


