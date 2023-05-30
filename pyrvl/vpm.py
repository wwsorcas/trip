from abc import ABC, abstractmethod
import vcl
import vcalg

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
        



            

    
