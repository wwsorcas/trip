import vcl
import vcalg

class vpmcgjet(vcl.ScalarFunction):

    '''
    Implementation of variable projection reduction for nonlinear
    separable least squares problems, using conjugate gradient
    approximate solution of the linear least squares subproblem

    

    Parameters:
        x0 (vcl.Vector): evaluation point in nonlinear subspace 
            of domain
        F (vcl.SepFunction): separable function object
        b (vcl.Vector): rhs vector in LS objective function
        kmax (int): max number of CG iterations
        rho (float): residual tolerance
        verbose (int):: verbosity level for CG
        
    '''

    def __init__(self, x0, F, b, kmax, rho, verbose=0):
        try:
            if not isinstance(x0, vcl.Vector):
                raise Exception('first input not vcl.Vector')
            if not isinstance(F,vcl.SepFunction):
                raise Exception('first input not vcl.SepFunction')
            if b.space != F.getRange():
                raise Exception('second arg not vector in range')
        except Exception as ex:
            print(ex)
            raise Exception('called from vpm.vpmcg constructor')
        else:
            self.x = x0
            self.F = F
            self.b = b
            self.kmax = kmax
            self.rho = rho
            self.verbose = verbose
            # storage for residual 
            self.r = None
            self.w = None
            self.g = None
            self.v = None

    def getDomain(self):
        return self.F.getDomain()[0]

    def getRange(self):
        return self.F.getRange()

    def value(self):
        try:
            if self.r is None:
                A = vcl.LinearRestriction(self.F,self.x)
                self.r = vcl.Vector(A.getRange())
                self.w = vcl.Vector(A.getDomain())
                vcalg.bcg(self.w, self.b, A, self.kmax, sefl.rho,
                              self.verbose, r=self.r)
            if self.v is None:
                self.v = 0.5*(self.r.norm()**2)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmcg.value')
        else:
            return self.v

    def gradient(self):
        try:
            if self.r is None:
                A = vcl.LinearRestriction(self.F,self.x)
                self.r = vcl.Vector(A.getRange())
                self.w = vcl.Vector(A.getDomain())
                vcalg.bcg(self.w, self.b, A, self.kmax, self.rho,
                              self.verbose, r=self.r)
            if self.g is None:
                G = vcl.DerivNonlinearRestriction(self.F,self.x,self.w)
                # r = b-Ax, grad = DA^T(Ax-b)
                self.g = transp(G)*self.r
                self.g.scale(-1.0)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmcg.value')
        else:
            return self.g

    # place holder
    def Hessian(self):
        return None
            
            

    
