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
        y = (a*x + b*y)
        return y
        
    def raw_dot(self,x,y):
        return np.dot(x.T,y)[0][0]

    def raw_scale(self,x,c):
        x=c*x
        return x

    # convenience functions
    def raw_copy(self,x,y):
        y=np.copy(x)
        return y

    # nothing to see here
    def cleanup(self,x):
        pass

    def raw_printData(self,x):
        print(x)
       
    def myNameIs(self):
        print('npvc.Space of dimension ' + str(self.dim))
    
class MatrixOperator(vcl.LinearOperator):    

    def __init__(self,dom,rng,mat):

        try: 
            self.dom = dom
            self.rng = rng
            matuple  = mat.shape
            if matuple[0] != rng.dim:
                raise Exception('num rows ne rng dim')
            if matuple[1] != dom.dim:
                raise Exception(' num cols ne dom dim')
            self.mat = np.copy(mat)
        except Exception as ex:
            print(ex)
            raise Exception('called from MatrixOperator constructor')
    
    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.rng
    
    def applyFwd(self,x, y):
        y.data = self.mat@x.data
        return y

    def applyAdj(self,x, y):
        y.data = self.mat.T@x.data
        return y

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

    def apply(self,x,y):
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
        return MatrixOperator(self.dom,self.rng,mat)

    def myNameIs(self):
        print('OpExpl1: npvc example of vcl.Function class')
        print('implements (x0,x1) -> (x0*x1, -x1+x0^2, x1^2)')
        print('domain = R^2, range = R^3')
        
class OpExpl2(vcl.Function):

    def __init__(self,dom,rng):
        try:
            # check that dom = R^1 oplus R^1
            if isinstance(dom,vcl.ProductSpace):
                if len(dom.spl) == 2:
                    if not isinstance(dom.spl[0],Space):
                        raise Exception('Error: dom.spl[0] not npvc.Space')
                    if not isinstance(dom.spl[1],Space):
                        raise Exception('Error: dom.spl[1] not npvc.Space')
                    if dom.spl[0].dim != 1 or dom.spl[1].dim !=1:
                        raise Exception('Error: dom factor dims wrong')
                else:
                    raise Exception('Error: dom has wrong num of factors')
            else:
                raise Exception('Error: dom not prod sp')
            # check that rng = R^3 
            if rng.dim !=3:
                raise Exception('Error: rng dim wrong')
        except Exception as ex:
            print(ex)
            raise Exception('called from npvc:OpExpl1')
        else:
            self.dom = dom
            self.rng = rng
        
    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def apply(self,x,y):
        y.data[0] = x.data[0][0]*x.data[1][0]
        y.data[1] = -x.data[1][0]+x.data[0][0]*x.data[0][0]
        y.data[2] = x.data[1][0]*x.data[1][0]
        return y

    def raw_deriv(self,x):
        oplist = []
        mat = np.zeros((3,1))
        mat[0] = x.data[1][0]
        mat[1] = 2*x.data[0][0]
        mat[2] = 0.0
        oplist.append(MatrixOperator(self.dom.spl[0],self.rng,mat))
        mat[0] = x.data[0][0]
        mat[1] = -1.0
        mat[2] = 2*x.data[1][0]
        oplist.append(MatrixOperator(self.dom.spl[1],self.rng,mat))
        return vcl.RowLinearOperator(self.dom,self.rng,oplist)

    def myNameIs(self):
        print('OpExpl2: npvc example of vcl.Function class')
        print('implements ((x0),(x1)) -> (x0*x1, -x1+x0^2, x1^2)')
        print('domain = R^1 oplus R^1, range = R^3')

# two copies of Rosenbrack function
class DoubleRosie(vcl.Function):

    def __init__(self,dom):
        try:
            self.dom = dom
            if dom.dim != 4:
                raise Exception('Error: input dims wrong')
        except Exception as ex:
            print(ex)
            raise Exception('called from npvc:DoubleRosie')
        
    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.dom

    def apply(self,x,y):
        y.data[0]=10*(x.data[1]-x.data[0]*x.data[0]);
        y.data[1]=-x.data[0];
        y.data[2]=2*(x.data[3]-x.data[2]*x.data[2]);
        y.data[3]=-x.data[2];
        return y
        
    def raw_deriv(self,x):
        mat = np.zeros((4,4))
        mat[0,0] = -20*x.data[0]
        mat[0,1] = 10
        mat[1,0] = -1.0
        mat[1,1] = 0.0
        mat[2,2] = -4*x.data[2]
        mat[2,3] = 2
        mat[3,2] = -1.0
        mat[3,3] = 0.0
        return MatrixOperator(self.dom,self.dom,mat)

    def myNameIs(self):
        print('DoubleRosie: npvc example of vcl.Function class')
        print('Rosenbrock 2x2, twice, scale factors 10 and 2')
        print('domain = R^4, range = R^4')
    

        
