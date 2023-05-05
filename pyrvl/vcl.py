# abstract space class
import math
import os
from abc import ABC, abstractmethod

class Space(ABC):

    # returns a data container handle ("data") - pointer, file name, whatever -
    # consistent with the data structure of vectors in this space
    @abstractmethod
    def getData(self):
        pass

    # returns True if data object is compatible, else False
    @abstractmethod
    def isData(self,x):
        pass
    
    # a, b are scalars, x, y are data objects
    @abstractmethod
    def raw_linComb(self, a, x, y, b=1.0):
        pass

    def linComb(self,a,x,y,b=1.0):
        try:
            if not self.isData(x):
                raise Exception('Error: 1st arg not space data object')
            if not self.isData(y):
                raise Exception('Error: 2nd arg not sapce data object')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Space.linComb')
        else:
            return self.raw_linComb(a,x,y ,b)

    # x, y are data objects
    @abstractmethod    
    def raw_dot(self,x,y):
        pass

    def dot(self,x,y):
        try:
            if not self.isData(x):
                raise Exception('Error: 1st arg not space data object')
            if not self.isData(y):
                raise Exception('Error: 2nd arg not space data object')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Space.dot')
        else:
            return self.raw_dot(x,y)
        
    # x is data object
    @abstractmethod    
    def raw_scale(self,x,c):
        pass

    def scale(self,x,c):
        try:
            if not self.isData(x):
                raise Exception('Error: arg not space data object')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Space.zero')
        else:            
            return self.raw_scale(x,c)

    # raw_copy should implement the un-Python like behaviour of
    # altering its second argument, by copying the first onto it.
    # so like numpy.copyto, not numpy.copy. 
    @abstractmethod        
    def raw_copy(self,x,y):
        pass

    # Note that this function returns None.
    def copy(self,x,y):
        try:
            if not self.isData(x):
                raise Exception('Error: 1st arg not space data object')
            if not self.isData(y):
                raise Exception('Error: 2nd arg not sapce data object')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Space.copy')
        else:        
            self.raw_copy(x,y)

    # for use in vector destructor - x is data
    # there is no "raw_cleanup" because this should be
    # called only by Vector destructor therefore only
    # on space data objects - so no error checking necessary.
    @abstractmethod        
    def cleanup(self,x):
        pass
    
    @abstractmethod                
    def raw_printData(self,x):
        pass

    def printData(self,x):
        try:
            if not self.isData(x):
                raise Exception('Error: arg not space data object')
        except Exception as ex:
            print(ex)
            raise Exception('from vcl.Space.printData')
        else:            
            self.raw_printData(x)
            
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

    def raw_linComb(self,a,x,y,b=1.0):
        for i in range(0,len(self.spl)):
            y[i]=self.spl[i].linComb(a,x[i],y[i],b)
        return y

    def raw_dot(self,x,y):
        ret = 0.0
        for i in range(0,len(self.spl)):
            ret += self.spl[i].dot(x[i],y[i])
        return ret

    def raw_scale(self,x,c):
        for i in range(0,len(self.spl)):
            x[i]-self.spl[i].scale(x[i],c)
        return x
    
    def raw_copy(self,x,y):
        for i in range(0,len(self.spl)):
            self.spl[i].copy(x[i],y[i])
    
    def cleanup(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].cleanup(x[i])
    
    def raw_printData(self,x):
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
        except Exception as ex:
            print(ex)
            raise Exception("called from vcl.Vector.link")
        else:
            print('in link: x = ' + x)
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
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.copy')
        else:
            self.space.copy(x.data,self.data)

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

    @abstractmethod
    def getDomain(self):
        pass
    
    @abstractmethod
    def getRange(self):
        pass

    @abstractmethod
    def apply(self,x,y):
        pass

    def __call__(self,x):
        try:
            #os.system('ls /var/tmp')
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            y = Vector(self.getRange())
            #os.system('ls /var/tmp')
            self.apply(x,y)
            #y.myNameIs()
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector operator()')
        else:
            return y

    # should return linear op
    @abstractmethod
    def raw_deriv(self,x):
        pass
    
    def deriv(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
        except Exception as ex:
            print(ex)
            # self.MyNameIs()
            raise Exception('called from vcl.Function.deriv')
        else:        
            return self.raw_deriv(x)

    @abstractmethod
    def myNameIs(self):
        pass

# function composition f after g
class comp(Function):

    def __init__(self,f,g):
        try:
            if f.getDomain() != g.getRange():
                raise Exception('Error: domain of output fcn ne range of input fcn')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.comp constructor')
        else:
            self.f = f
            self.g = g
    
    def getDomain(self):
        return self.g.getDomain()

    def getRange(self):
        return self.f.getRange()

    def apply(self,x,y):
        try:
            z = self.g(x)
            self.f.apply(z,y)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.comp.apply')
        else:
            return y

    def raw_deriv(self,x):
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
    
class LinearOperator(Function):

    def apply(self,x,y):
        self.applyFwd(x,y)

    @abstractmethod
    def applyFwd(self,x,y):
        pass

    # alternate evaluation syntax
    def __mul__(self,x):
        try:
#            print('in mul before:')
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            y = Vector(self.getRange())
#            y.myNameIs()
            self.applyFwd(x,y)
#            print('in mul after:')
#            y.myNameIs()
#            print('exit mul')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.LinearOperator operator*')
        else:
            return y

    @abstractmethod
    def applyAdj(self,x,y):
        pass
    
    def raw_deriv(self,x):
        return self

class transp(LinearOperator):

    def __init__(self,op):
        self.op = op

    def getDomain(self):
        return self.op.getRange()

    def getRange(self):
        return self.op.getDomain()

    def applyFwd(self,x,y):
        try:
#            print('transp.applyFwd')
#            print('|x|=' + str(x.norm()))
            self.op.applyAdj(x,y)
#            print('|y|=' + str(y.norm()))
#            print('exit transp.applyFwd')
        except Exception as ex:
            print(ex)
            raise Exception('called from transp:applyFwd')
        else:
            return y

    def applyAdj(self,x,y):
        try:
            self.op.applyFwd(x,y)
        except Exception as ex:
            print(ex)
            raise Exception('called from transp:raw_applyFwd')
        else:
            return y

    def myNameIs(self):
        print('adjoint operator of:')
        self.op.myNameIs()
    
class lopcomp(LinearOperator):

    def __init__(self,a,b):
        try:
#            print('lopcomp constructor')
#            print('first arg')
#            a.myNameIs()
#            print('second arg')
#            b.myNameIs()
            if a.getDomain() != b.getRange():
                raise Exception('Error: domain of output fcn ne range of input fcn')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp constructor')
        else:
            self.a = a
            self.b = b
    
    def getDomain(self):
        return self.b.getDomain()

    def getRange(self):
        return self.a.getRange()

    def applyFwd(self,x,y):
        try:
            z = self.b*x
            self.a.applyFwd(z,y)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp.applyFwd')

    def applyAdj(self,x,y):
        try:
            z=transp(self.a)*x
            self.b.applyAdj(z,y);
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.lopcomp.applyAdj')

    def myNameIs(self):
        print('vcl.lopcomp object: a composed with (after) b')
        print('composition of input lin op (b)')
        self.b.myNameIs()
        print('and output lin op (a)')
        self.a.myNameIs()
           
# makes a list of linear ops with same range act like a linear op on product
# space
# row of linear ops, acts on a column of vectors (vector in product space)
# to-do: make 1x1 case transparent
class RowLinearOperator(LinearOperator):

    def __init__(self,dom,rng,oplist):
        try:
            self.prod = True
            if not isinstance(dom,ProductSpace):
                self.prod = False
                if len(oplist) > 1:
                    raise Exception('Error: dom, oplist not both len=1')
            if len(dom.spl) != len(oplist):
                raise Exception('Error: dom, oplist diff lens')
            for i in range(0,len(oplist)):
                if dom.spl[i] != oplist[i].getDomain():
                    raise Exception('Error: op[' + str(i) + '] has wrong domain')
                if rng != oplist[i].getRange():
                    raise Exception('Error: op[' + str(i) + '] has wrong range')
        except Exception as ex:
            print(ex)
            raise Exception('called from RowLinearOp constructor') 
        else:
            self.dom = dom
            self.rng = rng
            if self.prod:
                self.oplist = oplist
            else:
                self.oplist = []
                self.oplist.append(oplist[0])

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def applyFwd(self,x,y):
        y.scale(0.0)
        yy = Vector(self.rng)
        for i in range(0,len(self.oplist)):
            if self.prod:
                xx = x.component(i)
                self.oplist[i].applyFwd(xx,yy)
                y.linComb(1.0,yy)
            else:
                self.oplist[0].applyFwd(x,y)
    
    def applyAdj(self,x,y):
        for i in range(0,len(self.oplist)):
            if self.prod:
                yy = y.component(i)
                self.oplist[i].applyAdj(x,yy)
            else:
                self.oplist[0].applyAdj(x,y)

    # magic component access
    def __getitem__(self,i):
        try:
            if i<0 or i>len(self.oplist)-1:
                raise Exception('Error: index out of range [0, ' + \
                                    str(len(self.oplist)-1) + ']')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.RowLinearOperator operator[]')
        else:
            return self.oplist[i]
        
    def myNameIs(self):
        print('RowLinearOperator length = ' + str(len(self.oplist)))
        for i in range(0,len(self.oplist)):
            print('*** Component ' + str(i) + ':')
            self.oplist[i].myNameIs()


class Jet:

    def __init__(self,f,x):
        try:
            if f.getDomain() != x.space:
                raise Exception('Jet: vec not in fcn dom')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Jet constr')
        else:
            self.point=vcl.Vector(f.getDomain())
            self.point.copy(x)
            self.value=vcl.Vector(f.getRange())
            f.apply(point,value)
            self.deriv=f.deriv(point)

    def getPoint(self):
        return self.point

    def getValue(self):
        return self.value

    def getDeriv(self):
        return self.deriv

    def myNameIs(self):
        print('vcl.Jet instance')
        print('point of evaluation:')
        self.point.myNameIs()
        print('value')
        self.value.myNameIs()
        print('derivative')
        self.deriv.myNameIs()

# simple nonlinear function class
class ScalarFunction(ABC):

    @abstractmethod
    def getDomain(self):
        pass
    
    @abstractmethod
    def value(self,x):
        pass

    def __call__(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.ScalarFunction operator()')
        else:        
            return self.value(x)

    # returns vector
    @abstractmethod
    def raw_gradient(self,x):
        pass
    
    def gradient(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            return self.raw_gradient(x)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.ScalarFunction gradient')
        
    # returns linear operator. Note required for some first order methods
    # (steepest descent, BFGS) so has default implementation returning None
    def raw_Hessian(self,x):
        return None
    
    def Hessian(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            H = self.raw_Hessian(x)
            if H is None:
                raise Exception('Error: Hessian not implemented for this class')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.ScalarFunction Hessian')
        else:
            return H
            
    def myNameIs(self):
        pass

# function x -> 0.5*|f(x)-b|^2
class LeastSquares(ScalarFunction):
    '''
    Simple nonlinear least squares function. Gauss-Newton Hessian supplied 
    by default, can be overridden. No attempt at retaining intermediate data
    common between methods - there is no assured way to do so. If such
    retention is advantageous, write a ScalarJet subclass instead - that
    allows for retention of intermediate data without risk of incoherence, 
    absent malicious interference.
    '''

    def __init__(self,f,b):
        self.f = f
        self.b = b

    def getDomain(self):
        return self.f.getDomain()

    def value(self,x):
        res=Vector(self.f.getRange())
        res=self.f(x)
        res.linComb(-1.0,self.b)
        return 0.5*res.dot(res)

    def raw_gradient(self,x):
        res = self.f(x)
        res.linComb(-1.0,self.b)
        df=self.f.deriv(x)
        return transp(df)*res

    def raw_Hessian(self,x):
        return NormalOp(self.f.deriv(x))

    def myNameIs(self):
        print('Least Squares Function')
        print('*** operator:')
        self.f.myNameIs()
        print('*** rhs vector')
        self.b.myNameIs()

class NormalOp(LinearOperator):

    def __init__(self,op):
        self.op = op

    def getDomain(self):
        return self.op.getDomain()

    def getRange(self):
        return self.op.getDomain()

    def applyFwd(self,x,y):
        #z = Vector(self.op.getRange())
        #self.op.applyFwd(x,z)
        #self.op.applyAdj(z,y)
        y = transp(op)*(op*x)

    def applyAdj(self,x,y):
        self.applyFwd(x,y)

    def myNameIs(self):
        print('normal operator of')
        self.op.myNameIs()

# separable function, aka linear-function-valued function. Not
# realized as a function itself, but as a repository of the
# components necessary to create the two restrictions that
# figure in VPM
class SepFunction(ABC):

    @abstractmethod    
    def getDomain(self):
        pass
    
    @abstractmethod    
    def getRange(self):
        pass

    # linear in x1
    @abstractmethod
    def applyFwd(self,x0,x1,y):
        pass

    # adjoint of x1 -> applyFwd(x0,x1,y)
    @abstractmethod
    def applyAdj(self,x0,y,x1):
        pass

    # returns partial derivative of applyFwd(x0,x1,y) in x0
    @abstractmethod
    def applyFwdDeriv(self, x0, dx0, x1, y):
        pass

    # adjoint partial derivative of applyFwd(x0,x1,y) in x0
    @abstractmethod    
    def applyAdjDeriv(self, x0, y, x1, dx0):
        pass

    @abstractmethod
    def myNameIs(self):
        pass

# linear restriction of separable function
class LinearRestriction(LinearOperator):

    def __init__(self, s, x0):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x0, Vector):
                raise Exception('second input not Vector')
            if x0.space != s.getDomain()[0]:
                raise Exception('second input not Vector in comp 0 of SepFunction domain')
        except Exception as ex:
            print(ex)
            raise Exception('called from LinearRestriction constructor')
        else:
            self.s = s
            self.x0 = x0

    def getDomain(self):
        return self.s.getDomain()[1]

    def getRange(self):
        return self.s.getRange()

    def applyFwd(self,x,y):
        self.s.applyFwd(x0,x,y)

    def applyAdj(self,x,y):
        self.s.applyAdj(x0,x,y)

    def myNameIs(self):
        print('Linear Restriction of Separable Function:')
        self.s.myNameIs()
        print('at nonlinear component vector:')
        self.x0.myNameIs()
        
# derivative of nonlinear restriction of separable function
class DerivNonlinearRestriction(LinearOperator):

    def __init__(self, s, x0, x1):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x0, Vector):
                raise Exception('second input not Vector')
            if x0.space != s.getDomain()[0]:
                raise Exception('second input not Vector in comp 0 of SepFunction domain')
            if not isinstance(x1, Vector):
                raise Exception('second input not Vector')
            if x1.space != s.getDomain()[1]:
                raise Exception('second input not Vector in comp 1 of SepFunction domain')
        except Exception as ex:
            print(ex)
            raise Exception('called from DerivNonlinearRestriction constructor')
        else:
            self.s = s
            self.x0 = x0
            self.x1 = x1

    def getDomain(self):
        return self.s.getDomain()[0]

    def getRange(self):
        return self.s.getRange()

    def applyFwd(self,dx0,y):
        self.s.applyFwdDeriv(self.x0, dx0, self.x1, y)

    def applyAdj(self,y,dx0):
        self.s.applyAdjDeriv(self.x0, y, self.x1, dx0)

    def myNameIs(self):
        print('Derivative of Nonlinear Restriction of Separable Function:')
        self.s.myNameIs()
        print('at nonlinear component vector:')
        self.x0.myNameIs()
        print('and linear component vector:')
        self.x1.myNameIs()

# Nonlinear restriction of separable function
class NonlinearRestriction(Function):

    def __init__(self, s, x1):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x1, Vector):
                raise Exception('second input not Vector')
            if x1.space != s.getDomain()[1]:
                raise Exception('second input not Vector in comp 1 of SepFunction domain')
        except Exception as ex:
            print(ex)
            raise Exception('called from LinearRestriction constructor')
        else:
            self.s = s
            self.x1 = x1

    def getDomain(self):
        return self.s.getDomain()[0]

    def getRange(self):
        return self.s.getRange()

    def apply(self,x0,y):
        self.s.applyFwd(x0,self.x1,y)

    def deriv(self,x0):
        return DerivNonlinearRestriction(self.s, x0, self.x1)

    def myNameIs(self):
        print('Nonlinear Restriction of Separable Function:')
        self.s.myNameIs()
        print('at nonlinear component vector:')
        self.x0.myNameIs()        
        
class ScalarJet(ABC):
    '''
    Abstract 2-jet class for scalar functions. Stores  
    copy of evaluation point to reduce dependence on calling
    unit. This choice makes accidental incoherence of attributes
    unlikely.
    '''

    @abstractmethod
    def point(self):
        pass

    @abstractmethod
    def value(self):
        pass

    @abstractmethod
    def gradient(self):
        pass

    @abstractmethod
    def Hessian(self):
        pass

    @abstractmethod
    def myNameIs(self):
        pass

    

