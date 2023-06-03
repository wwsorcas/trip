from abc import ABC, abstractmethod
import vcl
import vcalg
import numpy as np

class SepFunction(ABC):
    '''
    separable function, that is, a function of two args that is linear 
    in one of them:
    
    (x0,x1) \rightarrow A(x0)x1, 

    in which A is a linear-op-value function
    of x0. Not realized as a function itself, but as a repository of the
    components necessary to create the two restrictions that
    figure in VPM. 
    '''
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
class LinearRestriction(vcl.LinearOperator):

    def __init__(self, s, x0):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), vcl.ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x0, vcl.Vector):
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
        self.s.applyFwd(self.x0,x,y)

    def applyAdj(self,x,y):
        self.s.applyAdj(self.x0,x,y)

    def myNameIs(self):
        print('Linear Restriction of Separable Function:')
        self.s.myNameIs()
        print('at nonlinear component vector:')
        self.x0.myNameIs()
        
# derivative of nonlinear restriction of separable function
class DerivNonlinearRestriction(vcl.LinearOperator):

    def __init__(self, s, x0, x1):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), vcl.ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x0, vcl.Vector):
                raise Exception('second input not Vector')
            if x0.space != s.getDomain()[0]:
                raise Exception('second input not Vector in comp 0 of SepFunction domain')
            if not isinstance(x1, vcl.Vector):
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
class NonlinearRestriction(vcl.Function):

    def __init__(self, s, x1):
        try:
            if not isinstance(s, SepFunction):
                raise Exception('first input not SepFunction')
            if not isinstance(s.getDomain(), vcl.ProductSpace):
                raise Exception('domain of first input not ProductSpace')
            if len(s.getDomain().spl) != 2:
                raise Exception('number of components of domain != 2')
            if not isinstance(x1, vcl.Vector):
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
        print('at Linear component vector:')
        self.x1.myNameIs()        

class vpmcgjet(vcl.ScalarJet):

    '''
    Implementation of variable projection reduction for nonlinear
    separable least squares problems, using conjugate gradient
    approximate solution of the linear least squares subproblem

    Parameters:
        x0 (vcl.Vector): evaluation point in nonlinear subspace 
            of domain
        F (SepFunction): separable function object
        b (vcl.Vector): rhs vector in LS objective function
        kmax (int): max number of CG iterations
        eps (float): residual tolerance
        rho (float): normal residual tolerance
        verbose (int):: verbosity level for CGNE
        
    '''

    def __init__(self, x0, F, b, kmax, eps, rho, verbose=0):
        try:
            if not isinstance(x0, vcl.Vector):
                raise Exception('first arg not vcl.Vector')
            if x0.space != F.getDomain()[0]:
                raise Exception('first arg not in domain of second arg')
            if not isinstance(F, SepFunction):
                raise Exception('second arg not vcl.SepFunction')
            if b.space != F.getRange():
                raise Exception('third arg not vector in range of second arg')
        except Exception as ex:
            print(ex)
            raise Exception('called from vpm.vpmcg constructor')
        else:
            self.x = x0
            self.F = F
            self.b = b
            self.kmax = kmax
            self.eps = eps
            self.rho = rho
            self.verbose = verbose
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
                A = LinearRestriction(self.F,self.x)
                self.e = vcl.Vector(A.getRange())
                self.w = vcl.Vector(A.getDomain())
                vcalg.conjgrad(self.w, self.b, A, self.kmax,
                                   self.eps, self.rho, self.verbose,
                                   e=self.e)
            if self.v is None:
                self.v = 0.5*(self.e.norm()**2)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmcg.value')
        else:
            return self.v

    def gradient(self):
        try:
            if self.e is None:
                A = LinearRestriction(self.F,self.x)
                self.e = vcl.Vector(A.getRange())
                self.w = vcl.Vector(A.getDomain())
                vcalg.conjgrad(self.w, self.b, A, self.kmax,
                                   self.eps, self.rho, self.verbose,
                                   e=self.e)
            if self.g is None:
                G = DerivNonlinearRestriction(self.F,self.x,self.w)
                # r = b-Ax, grad = DA^T(Ax-b)
                self.g = vcl.transp(G)*self.e
                self.g.scale(-1.0)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmcg.gradient')
        else:
            return self.g

    # place holder
    def Hessian(self):
        return None

    def myNameIs(self):
        print('vpmcgjet blah blah')

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

    def applyFwd(self,x0,x1,y):
        xs = x0.norm()**2/(1 + x0.norm()**2)
        mat = self.mat0 + xs*self.mat1
        np.copyto(y.data,mat @ x1.data)
        
    def applyAdj(self,x0,y,x1):
        xs = x0.norm()**2/(1 + x0.norm()**2)
        mat = self.mat0 + xs*self.mat1
        np.copyto(x1.data, mat.T @ y.data)

    # returns partial derivative of applyFwd(x0,x1,y) in x0
    def applyFwdDeriv(self, x0, dx0, x1, y):
        dxs = 2*x0.dot(dx0)/((1+x0.norm()**2)**2)
        mat = dxs*self.mat1
        np.copyto(y.data,mat @ x1.data)

    # adjoint partial derivative of applyFwd(x0,x1,y) in x0
    def applyAdjDeriv(self, x0, y, x1, dx0):
        dxs = (2.0/((1+x0.norm()**2)**2))
        np.copyto(dx0.data,dxs*np.dot(y.data.T,self.mat1 @ x1.data)[0][0]*x0.data)

    def myNameIs(self):
        print('sepexpl1')



        



            

    
