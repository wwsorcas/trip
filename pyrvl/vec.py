# abstract space class

import math
from abc import ABC, abstractmethod

class Space(ABC):

    # returns a data container handle ("data") - pointer, file name, whatever -
    # consistent with the data structure of vectors in this space
    @abstractmethod
    def getData(self):
        pass

    # returns True if data object is compatible, else False
    @abstractmethod
    def compat(self,x):
        pass
    
    # a, b are scalars, x, y are data objects
    @abstractmethod
    def linComb(self, a, x, y, b=1.0):
        pass

    # x, y are data objects
    @abstractmethod    
    def dot(self,x,y):
        pass

    # x is data object
    @abstractmethod    
    def zero(self,x):
        pass

    # removes data objects compatible with this space
    @abstractmethod
    def cleanup(self,x):
        pass

    # prints name of a compatible data object
    @abstractmethod
    def printData(self,x):
        pass
    
    @abstractmethod
    def myNameIs(self):
        pass

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
        if self.space.compat(x):
            if self.own:
                self.space.cleanup(self.data)            
                self.data=x
            else:
                print("Vector.link: already managing non-owned data")
                return False
        else:
            print("Vector.link: data object not compatible with space")
            return False
        self.own = False
        return True
        
    def linComb(self,a,x,b=1.0):
        self.space.linComb(a,x.data,self.data,b)

    def dot(self,x):
        return self.space.dot(self.data,x.data)

    def zero(self):
        self.space.zero(self.data)

    def copy(self,x):
        self.space.copy(x.data,self.data)

    def norm(self):
        return math.sqrt(self.dot(self))

    def myNameIs(self):
        print('Vector in space:')
        self.space.myNameIs()
        print('Data object:')
        self.space.printData(self.data)

# simple linear op class

class LinearOperator:

    @abstractmethod
    def getDomain():
        pass
    
    @abstractmethod
    def getRange():
        pass
    
    @abstractmethod
    def applyFwd(x, y):
        pass

    @abstractmethod
    def applyAdj(x, y):
        pass
        
    @abstractmethod
    def myNameIs():
        pass

