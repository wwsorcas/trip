# abstract space class

import math
from abc import ABC, abstractmethod

class VCLExcept(Exception):
    pass

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

    # convenience functions
    @abstractmethod        
    def raw_copy(self,x,y):
        pass
    
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
            return self.raw_copy(x,y)

    # for use in vector destructor - x is data
    @abstractmethod        
    def raw_cleanup(self,x):
        pass
    
    def cleanup(self,x):
        try:
            if not self.isData(x):
                raise Exception('Error: arg not space data object')
        except Exception as ex:
            print(ex)
            raise Exception('from vcl.Space.printData')
        else:
            self.raw_cleanup(x)        

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
#            print(str(i))
            dl.append(self.spl[i].getData())
        return dl

    def isData(self,x):
        try:
            if not isinstance(x,list):
                raise Exception('isData arg not list')
            if len(x) < len(self.spl):
                raise Exception('isData arg is list but too short')
            ret = True
            for i in range(0,len(self.spl)):
                ret = ret and self.spl[i].isData(x[i])
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.ProductSpace.isData')
        else:
            return ret

    def raw_linComb(self,a,x,y,b=1.0):
        for i in range(0,len(self.spl)):
            y[i]=spl[i].linComb(a,x[i],y[i],b)
        return y

    def raw_dot(self,x,y):
        ret = 0.0
        for i in range(0,len(self.spl)):
            ret += spl[i].dot(x[i],y[i])
        return ret

    def raw_scale(self,x,c):
        for i in range(0,len(self.spl)):
            x[i]-self.spl[i].scale(x[i],c)
        return x
    
    def raw_copy(self,x,y):
        for i in range(0,len(self.spl)):
            y[i]=self.spl[i].copy(x[i],y[i])
        return y
    
    def raw_cleanup(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].cleanup(x[i])
    
    def raw_printData(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].printData(x[i])
            
    def myNameIs(self):
        print('vcl.ProductSpace')
        for i in range(0,len(self.spl)):
            print('*** compontent ' + str(i))
            self.spl[i].myNameIs()
        
# implemented vector class

class Vector:

    def __init__(self, sp):
        self.space = sp
        self.data = sp.getData()
        self.own = True

    def __del__(self):
        if self.own:
            self.space.cleanup(self.data)

    def link(self,x):
        try:
            if not self.space.isData(x):
                raise Exception('Error: attempt to wrap non-data object')
            if not self.own:
                raise Exception('Error: already managing external data')
        except Exception as ex:
            print(ex)
            raise Exception("called from vcl.Vector.link")
        else:
            self.space.cleanup(self.data)            
            self.data=x
            self.own = False
        
    def linComb(self,a,x,b=1.0):
        try:
            self.data=self.space.linComb(a,x.data,self.data,b)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.linComb')
        
    def dot(self,x):
        try:
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

    # this is assignment, but assignment "=" has fixed meaning
    # in Python, so this requires a function.
    def copy(self,x):
        try:
            self.data = self.space.copy(x.data,self.data)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.copy')            

    def norm(self):
        return math.sqrt(self.dot(self))
    
    # if not in ProdSp, 0th component is self
    # else data is list, so wrap ith component via link
    def component(self,i):
        try:
            if not isinstance(i,int):
                raise Exception('Error: input not int')
            if not isinstance(self.space,ProductSpace) and i>0:
                raise Exception('Error: not product and i>0')
            if i<0 or i>len(self.space.spl)-1:
                raise Exception('Error: index out of range')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.component')
        else:
            if not isinstance(self.space,ProductSpace):
                return self
            else:
                c = Vector(self.space.spl[i])
                c.link(self.data[i])
                return c
            
    def myNameIs(self):
        print('Vector in space:')
        self.space.myNameIs()
        print('Data object:')
        self.space.printData(self.data)

# simple nonlinear op class
class Function(ABC):

    @abstractmethod
    def getDomain(self):
        pass
    
    @abstractmethod
    def getRange(self):
        pass

    @abstractmethod
    def raw_apply(self,x,y):
        pass

    def apply(self,x,y):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            if y.space != self.getRange():
                raise Exception('Error: output vec not in range')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.Function.apply')
        else:        
            self.raw_apply(x,y)

    #### first attempt at overloaded evaluation
    def __call__(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            y = Vector(self.getRange())
            self.raw_apply(x,y)
        except Exception as ex:
            print(ex)
            raise Exception('called from Function()')
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
            self.MyNameIs()
            raise Exception('called from vcl.Function.deriv')
        else:        
            return self.raw_deriv(x)

#    def raw_partialDeriv(self,x,i):
#        pass

    def partialDeriv(self,x,i):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')            
            if not isinstance(self.getDomain(),ProductSpace):
                raise Exception('Error: domain not product space')
            if not isinstance(i,int):
                raise Exception('Error: index arg not integer')
            if i<0 or i>len(self.getDomain().spl)-1:
                raise Exception('Error: index ' + str(i) + ' out of range')
            dfx = self.deriv(x)
            if not isinstance(dfx,RowLinearOperator):
                raise Exception('Error: derivative not implemented as RowLinOp')
        except Exception as ex:
            print(ex)
            self.myNameIs()
            raise Exception('called from vcl.Function.partialDeriv')
        else:
            return dfx.oplist[i]
        
    @abstractmethod
    def myNameIs(self):
        pass
    
class LinearOperator(Function):

    def raw_apply(self,x,y):
        self.raw_applyFwd(x,y)

    @abstractmethod
    def raw_applyFwd(self,x,y):
        pass

    def applyFwd(self,x, y):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            if y.space != self.getRange():
                raise Exception('Error: output vec not in range')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.LinearOperator.applyFwd')
        else:        
            self.raw_applyFwd(x,y)

    # alternate evaluation syntax
    def __mul__(self,x):
        return self(x)

    @abstractmethod
    def raw_applyAdj(self,x,y):
        pass
    
    def applyAdj(self,x, y):
        try:
            if y.space != self.getDomain():
                raise Exception('Error: output vec not in domain')
            if x.space != self.getRange():
                raise Exception('Error: input vec not in range')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.LinearOperator.applyAdj')
        else:        
            self.raw_applyAdj(x,y)     

    def raw_deriv(self,x):
        return self

class transp(LinearOperator):

    def __init__(self,op):
        self.op = op

    def getDomain(self):
        return self.op.getRange()

    def getRange(self):
        return self.op.getDomain()

    def raw_applyFwd(self,x,y):
        self.op.applyAdj(x,y)
        return y

    def raw_applyAdj(self,x,y):
        self.op.applyFwd(x,y)
        return y

    def myNameIs(self):
        print('adjoint operator of:')
        self.op.myNameIs()
    

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

    def raw_applyFwd(self,x,y):
        y.scale(0.0)
        yy = Vector(self.rng)
        for i in range(0,len(self.oplist)):
            if self.prod:
                xx = x.component(i)
                self.oplist[i].applyFwd(xx,yy)
                y.linComb(1.0,yy)
            else:
                self.oplist[0].applyFwd(x,y)
    
    def raw_applyAdj(self,x,y):
        for i in range(0,len(self.oplist)):
            if self.prod:
                yy = y.component(i)
                self.oplist[i].applyAdj(x,yy)
            else:
                self.oplist[0].applyAdj(x,y)

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
    def raw_value(self,x):
        pass

    def value(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.ScalarFunction.apply')
        else:        
            return self.raw_value(x)

    # should return vector
    @abstractmethod
    def raw_gradient(self,x,g):
        pass
    
    def gradient(self,x,g):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
            if g.space != self.getDomain():
                raise Exception('Error: output vec not in domain')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.Function.deriv')
        else:        
            self.raw_gradient(x,g)
        
    @abstractmethod
    def myNameIs(self):
        pass

# function x -> 0.5*|f(x)-b|^2
class LeastSquares(ScalarFunction):

    def __init__(self,f,b):
        self.f = f
        self.b = b

    def getDomain(self):
        return self.f.getDomain()

    def raw_value(self,x):
        res=Vector(self.f.getRange())
        self.f.apply(x,res)
        res.linComb(-1.0,self.b)
        return 0.5*res.dot(res)

    def raw_gradient(self,x,g):
        res=Vector(self.f.getRange())
        self.f.apply(x,res)
        res.linComb(-1.0,self.b)
        df=self.f.deriv(x)
        df.applyAdj(res,g)

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

    def raw_applyFwd(self,x,y):
        z = Vector(self.op.getRange())
        self.op.applyFwd(x,z)
        self.op.applyAdj(z,y)

    def raw_applyAdj(self,x,y):
        self.applyFwd(x,y)

    def myNameIs(self):
        print('normal operator of')
        self.op.myNameIs()


