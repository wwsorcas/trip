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
            self.raw_linComb(a,x,y,b)

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
    def raw_zero(self,x):
        pass

    def zero(self,x):
        try:
            if not self.isData(x):
                raise Exception('Error: arg not space data object')
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Space.zero')
        else:            
            self.raw_zero(x)

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
            self.raw_copy(x,y)

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
            print(str(i))
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
            spl[i].linComb(a,x[i],y[i],b)

    def raw_dot(self,x,y):
        ret = 0.0
        for i in range(0,len(self.spl)):
            ret += spl[i].dot(x[i],y[i])
        return ret

    def raw_zero(self,x):
        for i in range(0,len(self.spl)):
            self.spl[i].zero(x[i])

    def raw_copy(self,x,y):
        for i in range(0,len(self.spl)):
            self.spl[i].copy(x[i],y[i])
        
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

    # for this to work the data object type must have
    # assignment semantics
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
            self.space.linComb(a,x.data,self.data,b)
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
        
    def zero(self):
        self.space.zero(self.data)

    def copy(self,x):
        try:
            self.space.copy(x.data,self.data)
        except Exception as ex:
            print(ex)
            raise Exception('called from vcl.Vector.copy')            

    def norm(self):
        return math.sqrt(self.dot(self))

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

    def raw_partialDeriv(self,x,i):
        pass

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
        except Exception as ex:
            print(ex)
            self.myNameIs()
            raise Exception('called from vcl.Function.partialDeriv')
        else:
            return self.raw_partialDeriv(self,x,i)
        
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
    def raw_apply(self,x):
        pass

    def apply(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.ScalarFunction.apply')
        else:        
            return self.raw_apply(x)

    # should return vector
    @abstractmethod
    def raw_grad(self,x):
        pass
    
    def grad(self,x):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')
        except Exception as ex:
            print(ex)
            self.MyNameIs()
            raise Exception('called from vcl.Function.deriv')
        else:        
            return self.raw_grad(x)

    def raw_partialGrad(self,x,i):
        pass

    def partialGrad(self,x,i):
        try:
            if x.space != self.getDomain():
                raise Exception('Error: input vec not in domain')            
            if not isinstance(self.getDomain(),ProductSpace):
                raise Exception('Error: domain not product space')
            if not isinstance(i,int):
                raise Exception('Error: index arg not integer')
            if i<0 or i>len(self.getDomain().spl)-1:
                raise Exception('Error: index ' + str(i) + ' out of range')
        except Exception as ex:
            print(ex)
            self.myNameIs()
            raise Exception('called from vcl.ScalarFunction.partialGrad')
        else:
            return self.raw_partialGrad(self,x,i)
        
    @abstractmethod
    def myNameIs(self):
        pass

    


