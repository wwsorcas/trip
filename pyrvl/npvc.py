import vcl
import numpy as np
import sys

# numpy vector space
class Space(vcl.Space):

    def __init__(self,n):
        self.dim = n

    def getData(self):
        return np.zeros(self.dim).reshape(self.dim,1)
    
    def isData(self,x):
        return (isinstance(x,np.ndarray) and x.shape == (self.dim,1))

    def raw_linComb(self,a,x,y,b=1.0):
        y = a*x + b*y

    def raw_dot(self,x,y):
        return np.dot(x.T,y)[0][0]

    def raw_zero(self,x):
        zip=0.0
        x=zip*x

    # convenience functions
    def raw_copy(self,x,y):
        y=x.copy()

    # for use in vector destructor - x is data 
    def raw_cleanup(self,x):
        pass

    def raw_printData(self,x):
        print(x)
       
    def myNameIs(self):
        print('npvc.Space of dimension ' + str(self.dim))
    
class MatrixOperator(vcl.LinearOperator):    

    def __init__(self,dom,rng,mat):

        self.dom = dom
        self.rng = rng
        matuple  = mat.shape
        if matuple[0] != rng.dim:
            print('MatrixOperator constructor: num rows ne rng dim')
            return False
        if matuple[1] != dom.dim:
            print('MatrixOperator constructor: num cols ne dom dim')
            return False
        self.mat = np.copy(mat)
    
    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.rng
    
    def raw_applyFwd(self,x, y):
        y.data = self.mat@x.data
            
    def raw_applyAdj(self,x, y):
        y.data = self.mat.T@x.data

    def myNameIs(self):
        print('NUMPY Matrix Operator with matrix:')
        print(self.mat)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()

class OpExpl1(vcl.Function):

    def __init__(self,dom,rng):
        try:
            self.dom = dom
            self.rng = rng
            if dom.dim != 2 or rng.dim !=3:
                raise Exception('Error: input dims wrong')
        except Exception as ex:
            print(ex)
            raise Exception('called from npvc:OpExpl1')
        
    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def raw_apply(self,x,y):
        y.data[0] = x.data[0]*x.data[1]
        y.data[1] = -x.data[1]+x.data[0]*x.data[0]
        y.data[2] = x.data[1]*x.data[1]

    def raw_deriv(self,x):
        mat = np.zeros((3,2))
        mat[0,0] = x.data[1]
        mat[0,1] = x.data[0]
        mat[1,0] = 2*x.data[0]
        mat[1,1] = -1.0
        mat[2,0] = 0.0
        mat[2,1] = 2*x.data[1]
        return npvc.MatrixOperator(self.dom,self.rng,mat)

    def myNameIs(self):
        print('OpExpl1: npvc example of vcl.Function class')
        
    

        

