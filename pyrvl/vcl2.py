# abstract space class
import math
import os
from functools import partial
from abc import ABC, abstractmethod

class Space(ABC):

    # returns a data container handle ("data") - pointer, file name, whatever -
    # consistent with the data structure of vectors in this space
    @abstractmethod
    def	getData(self):
    	pass

    # returns true if x is a data object associated with this space
    def isData(self,x):
    	pass

    # a, b are scalars, x, y are data objects
    @abstractmethod
    def linComb(self, a, x, y, b=1.0):
        pass

    # x, y are data objects
    @abstractmethod    
    def dot(self,x,y):
        pass

    # convenience functions

    # x is data object
    @abstractmethod    
    def scale(self,x,c):
        pass

    # copies x onto y
    @abstractmethod        
    def copy(self,x,y):
        pass
    
    # for use in Vector destructor
    @abstractmethod        
    def cleanup(self,x):
        pass

    # for use in Vector.myNameIs
    @abstractmethod                
    def printData(self,x):
        pass

    @abstractmethod
    def myNameIs(self):
        pass

class ProductSpace(Space):

    def __init__(self,SpaceList):
        if not isinstance(SpaceList,list):
            raise Exception('constructor arg not list')
        self.spl=[]
        for i in range(0,len(SpaceList)):
            sp = SpaceList[i]
            if not isinstance(sp,Space):
                raise Exception('arglist[' + str(i) + ' not Space')
            self.spl.append(sp)

    def getData(self):
        dl=[]
        for i in range(0,len(self.spl)):
            dl.append(self.spl[i].getData())
        return dl

    def isData(self,x):
        try:
            if not isinstance(x,list):
                raise Exception('isData arg not list')
            if len(x) != len(self.spl):
                raise Exception('isData arg is list but wrong length')
            ret = True
            for i in range(0,len(self.spl)):
                ret = ret and self.spl[i].isData(x[i])
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.ProductSpace.isData')
        else:
            return ret


    def __getitem__(self,i):
        try:
            if not isinstance(i,int):
                raise Exception('Error: index not int')
            if i<0 or i>len(self.spl)-1:
                raise Exception('Error: index out of range [0,' + \
                                    str(len(self.spl)-1) + ']')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.ProductSpace operator[]')
        else:
            return self.spl[i]            

    def linComb(self,a,x,y,b=1.0):
        for i in range(0,len(self.spl)):
            y[i]=spl[i].linComb(a,x[i],y[i],b)
        return y

    def dot(self,x,y):
        ret = 0.0
        for i in range(0,len(self.spl)):
            ret += spl[i].dot(x[i],y[i])
        return ret

    def scale(self,x,c):
        for i in range(0,len(self.spl)):
            x[i]-self.spl[i].scale(x[i],c)
        return x
    
    def copy(self,x,y):
        for i in range(0,len(self.spl)):
            y[i]=self.spl[i].copy(x[i],y[i])
        return y
    
    def cleanup(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].cleanup(x[i])
    
    def printData(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].printData(x[i])
            
    def myNameIs(self):
        print('vcl.ProductSpace')
        for i in range(0,len(self.spl)):
            print('*** component ' + str(i))
            self.spl[i].myNameIs()
        
# implemented vector class

class Vector:

    def __init__(self, sp, data=None):
        try:
            self.space = sp
            if data is None:
                self.data = sp.getData()
                self.own = True
            else:
                if self.space.isData(data):
                    self.data = data
                    self.own = False
                else:
                    raise Exception('Error: data is not compatible with Space')
            self.comp = []
            if isinstance(self.space,ProductSpace):
                for i in range(0,len(self.space.spl)):
                    self.comp.append(Vector(self.space.spl[i],self.data[i]))
            else:
                self.comp.append(self)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector constructor')

    def __del__(self):
        try:
            if self.own:
                # print('call self.space.cleanup')
                # print('os=')
                # print(os)
                self.space.cleanup(self.data)
        except Exception as ex:
            print(str(ex))
            raise Exception('called from vcl.Vector destructor')
        
    def link(self,x):
        try:
            if not self.space.isData(x):
                raise Exception('Error: attempt to wrap non-data object')
#            if not self.own:
#                raise Exception('Error: already managing external data')
        except Exception as ex:
            print(ex)
            raise Exception("called from vcl.Vector.link")
        else:
            if self.own:
                self.space.cleanup(self.data)            
            self.data = x
            self.own = False
        
    def linComb(self,a,x,b=1.0):
        try:
            if x.space != self.space:
                raise Exception('Error: vector summand not in same space')
            self.data=self.space.linComb(a,x.data,self.data,b)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.linComb')
        
    def dot(self,x):
        try:
            if x.space != self.space:
                raise Exception('Error: vector factor not in same space')
            dotp = self.space.dot(self.data,x.data)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.dot')
        else:
            return dotp
        
    def scale(self,c):
        try:
            self.data = self.space.scale(self.data,c)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.scale')

    # copy data of argument onto self data
    def copy(self,x):
        try:
            if x.space != self.space:
                raise Exception('Error: vector source not in same space')
            self.data = self.space.copy(x.data,self.data)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.copy')

    # instantiate new vector, copy self data, return
    def dup(self):
        try:
            x = Vector(self.space)
            x.data = self.space.copy(self.data,x.data)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.dup')
        else:
            return x

    def norm(self):
        return math.sqrt(self.dot(self))
    
    # if not in ProdSp, 0th component is self
    # else data is list, so wrap ith component via link
    def component(self,i):
        try:
            if not isinstance(i,int):
                raise Exception('Error: index not int')
            if not isinstance(self.space,ProductSpace) and i>0:
                raise Exception('Error: not product and i>0')
            if i<0 or i>len(self.space.spl)-1:
                raise Exception('Error: index out of range [0,' + \
                                    str(len(self.space.spl)-1) + ']')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.component')
        else:
            return self.comp[i]

    # magic version
    def __getitem__(self,i):
        try:
            return self.component(i)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector operator[]')
            
    def myNameIs(self):
        print('Vector in space:')
        self.space.myNameIs()
        print('Data object:')
        self.space.printData(self.data)
        # print('Owned by vector = ' + str(self.own))


# simple nonlinear op class
class Function(ABC):

    def __init__(self,dom,rng):
        self.dom = dom
        self.rng = rng

    @abstractmethod
    def __call__(self,x):
    	pass

    # should return linear op
    @abstractmethod
    def deriv(self,x):
        pass
    
    @abstractmethod
    def myNameIs(self):
        pass

# most functions will be built this way
class StdFunction(Function):
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
            if x.space != self.dom:
                raise Exception('Error: input vec not in domain')
            y = vcl.Vector(self.rng)
            if self.aux is None:
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

# function composition f after g
class FunctionComp(Function):

    def __init__(self,f,g):
        try:
            if f.dom != g.rng:
                raise Exception('Error: domain of output fcn ne range of input fcn')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.comp constructor')
        else:
            self.f = f
            self.g = g
	    self.dom = g.dom
	    self.rng = f.rng

    def __call__(self,x):
        try:
            if x.space != self.dom:
                raise Exception('Error: input vec not in domain')
            y = self.g(x)
            z = self.f(y)
	except Exception as ex:
	    print(ex)
	    raise Exception('called from FunctionComp operator()')
	else:
	    return z
    
    def deriv(self,x):
        try:
            gx = self.g(x)
            dfgx = self.f.deriv(gx)
            dg=self.g.deriv(x)
            return lopcomp(dfgx, dg) 
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.comp.raw_deriv')

    def myNameIs(self):
        print('vcl.comp object: f composed with (after) g')
        print('composition of input function (g)')
        self.g.myNameIs()
        print('and output function (f)')
        self.f.myNameIs()

class LinearFunction(Function):

    def __call__(self,x):
        return __mul__(x)

    # application - *
    @abstractmethod
    def __mul__(self,x):
    	pass

    # application of adjoint 
    @abstractmethod
    def tmul(self,x):
    	pass

    def deriv(self,x):
        return self

class StdLinearFunction(LinearFunction):

    '''
    applyFwd and applyAdj have same signature as apply in Function
    '''

    def __init__(self,dom,rng,applyFwd,applyAdj,aux=None):

        self.dom = dom
        self.rng = rng
        self.applyFwd = applyFwd
        self.applyAdj = applyAdj
        self.aux = aux    

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

    def tmul(self,x):
        try:
            if x.space != self.rng:
                raise Exception('Error: input vec not in range')
            y = vcl.Vector(self.dom)
            if self.aux is None:
                y.data=self.applyAdj(x.data,y.data)
            else:
                y.data=self.applyAdj(x.data,y.data,**self.aux)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.StdLinearFunction.tmul')
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
        self.lf = lf

    def __mul__(self,x):
        try:
            if x.space != self.dom:
                raise Exception('Error: input vec not in domain')
            y = self.lf.tmul(x)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.transp operator()')
        else:
            return y

    def tmul(self,x):
        try:
            if x.space != self.rng:
                raise Exception('Error: input vec not in range')
            y = self.lf*x
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.transp.tmul')
        else:
            return y
	    

    def myNameIs(self):
        print('\nTranspose of Linear Function:')
        self.lf.myNameIs()

class lopcomp(LinearOperator):

    def __init__(self,a,b):
        try:
            if a.dom != b.rng:
                raise Exception('Error: domain of output fcn ne range of input fcn')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp constructor')
        else:
            self.a = a
            self.b = b
	    self.dom = b.dom
	    self.rng = a.rng
    
    def __mul__(self,x,y):
        try:
            z = self.b*x
            y = self.a*z
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp.applyFwd')
	else:
	    return y

    def tmul(self,x:
        try:
            z=transp(self.a)*x
            y=transp(self.b)*z
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp.applyAdj')
        else:
	    return y

    def myNameIs(self):
        print('vcl.lopcomp object: a composed with (after) b')
        print('composition of input lin op (b)')
        self.b.myNameIs()
        print('and output lin op (a)')
        self.a.myNameIs()


