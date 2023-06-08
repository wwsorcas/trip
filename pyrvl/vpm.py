from abc import ABC, abstractmethod
import vcl
import vcalg
import numpy as np
import typing

class SepFunction(ABC):
    '''
    separable function, that is, a function of two args that is linear 
    in one of them:
    
    (x0,x1) \rightarrow A(x0)x1, 

    in which A is a linear-op-value function
    of x0. Since linear operators do not have a natural Hilbert space 
    structure, it is not possible to realize this concept as a vcl.Function,
    which by convention has a Hilbert space (really, vcl.Space) range. 
    Instead, this class has twp methods, each returning a linear operator:
    the function x0 -> A(x0), and the x0-partial derivative of the function 
    (x0,x1) -> A(x0)x1, that is, (x0,x1) -> D_x0(A(x0)x1).
    '''

    # return produce space in which (x0,x1) lies
    @abstractmethod
    def getDomain(self):
        pass

    # range of linear operator A(x0)
    @abstractmethod
    def getRange(self):
        pass
    
    @abstractmethod
    def opfcn(self, x0) -> vcl.LinearOperator:
        pass

    @abstractmethod
    def derfcn(self, x0, x1):
        pass
    
    @abstractmethod
    def myNameIs(self):
        pass        

class vpmjet(vcl.ScalarJet):

    '''
    Implementation of variable projection reduction for nonlinear
    separable least squares problems, using a vcl.Function S to supply
    approximate solution of the linear least squares subproblem.

    Parameters:
        x0 (vcl.Vector): evaluation point in nonlinear subspace 
            of domain
        F (SepFunction): separable function object
        b (vcl.Vector): rhs vector in LS objective function
        S (vcl.LSsolver): LS solution function
        kmax (int): max number of CG iterations
        eps (float): residual tolerance
        rho (float): normal residual tolerance
        verbose (int):: verbosity level for CGNE
        
    '''

    def __init__(self, x0, F, b, S):
        try:
            if not isinstance(x0, vcl.Vector):
                raise Exception('first arg not vcl.Vector')
            if not isinstance(F, SepFunction):
                raise Exception('second arg not SepFunction')
            if not isinstance(b, vcl.Vector):
                raise Exception('third arg not vcl.Vector')
            if not isinstance(S, vcl.LSSolver):
                raise Exception('fourth arg not vcl.LSSolver')
            if x0.space != F.getDomain()[0]:
                raise Exception('first arg not in domain of second arg')
            if b.space != F.getRange():
                raise Exception('third arg not vector in range of second arg')
        except Exception as ex:
            print(ex)
            raise Exception('called from vpm.vpmcg constructor')
        else:
            self.x = x0
            self.F = F
            self.b = b
            self.S = S
            # storage for residual 
            self.e = None
            self.w = None
            self.g = None
            self.v = None

    def point(self):
        return self.x

    def value(self):
        try:
            if self.e is None:
                [self.w, self.e] = self.S.solve(self.F.opfcn(self.x),self.b)
            if self.v is None:
                self.v = 0.5*(self.e.norm()**2)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmjet.value')
        else:
            return self.v

    def gradient(self):
        try:
            if self.e is None:
                [self.w, self.e] = self.S.solve(self.F.opfcn(self.x),self.b)
            if self.g is None:
                #G = Der(self.F,self.x,self.w)
                # r = b-Ax, grad = DA^T(Ax-b)
                self.g = vcl.transp(self.F.derfcn(self.x,self.w))*self.e
#                self.g.scale(-1.0)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmjet.gradient')
        else:
            return self.g

    # place holder
    def Hessian(self):
        return None

    def myNameIs(self):
        print('vpmjet blah blah')

############## EXAMPLES ##############

import npvc

class expl1(vcl.ScalarFunction):

    def __init__(self,sp):
        try:
            if not isinstance(sp,npvc.Space):
                raise Exception('input not npvc.Space')
            self.n = sp.dim
            self.dom = sp
            self.mat0 = np.zeros((2*self.n,self.n))
            self.mat1 = np.zeros((2*self.n,self.n))            
            for i in range(self.n):
                self.mat0[i][i] = 1.0
                self.mat1[self.n+i][i] = i+1
            self.rhs = np.zeros((2*self.n,1))
            for i in range(self.n):
                self.rhs[i]=(-i-1)**(i+1)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmtest.expl1 constructor')
        
    def getDomain(self):
        return self.dom

    def value(self,x):
        xs = x.norm()**2/(1 + x.norm()**2)
        mat = self.mat0 + xs*self.mat1
        [w, res, rk, s] = np.linalg.lstsq(mat,self.rhs, rcond=None)
        r = mat @ w - self.rhs
        return 0.5*np.dot(r.T,r)[0][0]

    def raw_gradient(self,x):
        xs = x.norm()**2/(1+x.norm()**2)
        mat = self.mat0 + xs*self.mat1
        [w, res, rk, s] = np.linalg.lstsq(mat,self.rhs, rcond=None)
        r = mat @ w - self.rhs
        g = vcl.Vector(self.getDomain())
#        p = 0.0
#        for i in range(self.n):
#            p = p + (i+1)*w[i]*r[n+i]
        p = np.dot(r.T,self.mat1 @ w)[0][0]
        np.copyto(g.data, 2.0*x.data*p/((1+x.norm()**2)**2))
        return g

    def raw_Hessian(self,x):
        return None

    def myNameIs(self):
        print('vpm expl 1, n = ' + str(self.n))

class jetexpl1(vcl.ScalarJet):

    def __init__(self,x):
        try:
            if not isinstance(x,vcl.Vector):
                raise Exception('input not vcl.Vector')
            if not isinstance(x.space,npvc.Space):
                raise Exception('input Vector not member of npvc.Space')
            self.x = x
            self.n = x.space.dim
            self.mat0 = np.zeros((2*self.n,self.n))
            self.mat1 = np.zeros((2*self.n,self.n))            
            for i in range(self.n):
                self.mat0[i][i] = 1.0
                self.mat1[self.n+i][i] = i+1
            self.rhs = np.zeros((2*self.n,1))
            for i in range(self.n):
                self.rhs[i]=(-i-1)**(i+1)
            ## trigger for internal computation
            self.w = None
            self.r = None
            self.v = None
            self.g = None
            self.H = None
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmtest.expl1 constructor')
        
    def point(self):
        return self.x

    # internal computations
    def nobodysbiz(self):
        if self.r is None:
            print('--> nobodysbiz')
            xs = self.x.norm()**2/(1 + self.x.norm()**2)
            mat = self.mat0 + xs*self.mat1
            [self.w, res, rk, s] = np.linalg.lstsq(mat,self.rhs, rcond=None)
            self.r = mat @ self.w - self.rhs        

    def value(self):
        self.nobodysbiz()
        if self.v is None:
            self.v = 0.5*np.dot(self.r.T,self.r)[0][0]
        return self.v
    
    def gradient(self):
        self.nobodysbiz()
        if self.g is None:
            self.g = vcl.Vector(self.x.space)
            p = np.dot(self.r.T,self.mat1 @ self.w)[0][0]
            np.copyto(self.g.data,
                          2.0*self.x.data*p/((1+self.x.norm()**2)**2))
        return self.g

    # temporary
    def Hessian(self,x):
        return None

    def myNameIs(self):
        print('vpm jet expl 1, n = ' + str(self.n))
        
###################

class sepexpl1(SepFunction):

    def __init__(self,dom,rng):
        try:
            if not isinstance(dom,vcl.ProductSpace):
                raise Exception('input dom not ProductSpace')
            if not isinstance(dom.spl[0],npvc.Space):
                raise Exception('input dom[0] not npvc.Space')
            if not isinstance(dom.spl[1],npvc.Space):
                raise Exception('input dom[1] not npvc.Space')
            if not isinstance(rng,npvc.Space):
                raise Exception('input rng not npvc.Space')
            if dom.spl[0].dim != dom.spl[1].dim:
                raise Exception('dom.spl[0].dim != dom.spl[1].dim')
            if rng.dim != 2*dom.spl[1].dim:
                raise Exception('rng.dim != 2*dom[1].dim')
            self.n = dom.spl[1].dim
            self.dom = dom
            self.rng = rng
            self.mat0 = np.zeros((2*self.n,self.n))
            self.mat1 = np.zeros((2*self.n,self.n))
            for i in range(self.n):
                self.mat0[i][i] = 1.0
                self.mat1[self.n+i][i] = i+1
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmtest.expl1 constructor')

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng
        
    def opfcn(self, x0):
        xs = x0.norm()**2/(1 + x0.norm()**2)
        mat = self.mat0 + xs*self.mat1
        return npvc.MatrixOperator(self.dom.spl[1],self.rng,mat)
        
    # returns partial derivative of applyFwd(x0,x1,y) in x0
    def derfcn(self, x0, x1):
#        dxs = 2*x0.dot(dx0)/((1+x0.norm()**2)**2)
# m x 1
        f1 = self.mat1 @ x1.data
# 1 x n
        f2 =  2*x0.data.T/((1+x0.norm()**2)**2)
        mat = f1 @ f2
#        mat = dxs*self.mat1
#        np.copyto(y.data,mat @ x1.data)
        return npvc.MatrixOperator(self.dom.spl[0],self.rng,mat)

    def myNameIs(self):
        print('sepexpl1')



        



            

    
