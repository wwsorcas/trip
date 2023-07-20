import linalg
import vcl
import vcalg
import vpm
import segyvc
import rsfvc
import os
import tempfile

TRIP = os.getenv('TRIP')

def setrms(u, urms):
    '''
    Calculate and store L2 norms of each trace in a SEGY file.
    Parameters:
    u (string): filename of input data traces
    urms (string): filename output data traces of length 1, 
    each containing L2 norm of corresponding trace in input file.
    '''
    try:
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
        if linalg.sanity(u,'su') and linalg.sanity(urms,'su'):
            cmd = os.path.join(TRIP,'iwave/trace/main/rms.x')
            ret = os.system(cmd + ' in=' + u + ' out=' + urms)
            if ret != 0:
                raise Exception('iwave/trace/main/rms.x failed')
        else:
            raise Exception('at least one input does not have .su suffix')
    except Exception as ex:
        print(ex)
        raise Exception('called from setrms')

### as of 2023.07.19, this class is deprecated in favor of compawipensol

def compawipen(u, p, u0rms, precond, alpha):
    '''
    Compute action of AWI penalty operator on a set of traces. Multiplies 
    trace samples by time and penalty weight, optionally scales by 
    reciprocal of trace-dependent scale factor, intended to be L2 norms 
    of traces for alpha=0 AWI source.
    Parameters:
    u (string): filename of input SEGY traces
    p (string): filename of output SEGY traces, same geometry as
    input traces (u)
    u0rms (string): filename of auxiliary length-1 traces containing
    scale factors for preconditioning. Same number of traces as input
    (u) and output (p)
    precond (int): preconditioning flag. 0 = no scaling, else scaling
    applied
    alpha (float): penalty weight
    '''
    try:
        #print('compawipen: u=' + u + ' p=' + p + ' u0rms=' + u0rms)
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
        if linalg.sanity(u,'su') and linalg.sanity(p,'su'):
            cmd = os.path.join(TRIP,'iwave/trace/main/awipen.x')
            if precond !=0 and linalg.sanity(u0rms,'su'):
                ret = os.system(cmd + ' in=' + u + ' rms=' + u0rms
                                    + ' out=' + p + ' precond=' + str(precond)
                                    + ' alpha=' + str(alpha))
            else:
                #print(cmd + ' in=' + u  
                #            + ' out=' + p + ' precond=0 alpha='
                #            + str(alpha))
                ret = os.system(cmd + ' in=' + u  
                                    + ' out=' + p + ' precond=0 alpha='
                                    + str(alpha))
            if ret != 0:
                raise Exception('iwave/trace/main/awipen.x failed')
        else:
            raise Exception('at least one input does not have .su suffix')
    except Exception as ex:
        print(ex)
        raise Exception('called from compawipen')

def compawipensol(u, p, alpha,  u0rms=None):
    '''
    Compute action of AWI penalty operator on a set of traces. Multiplies 
    trace samples by time and penalty weight, optionally scales by 
    reciprocal of trace-dependent scale factor, intended to be L2 norms 
    of traces for alpha=0 AWI source.
    Parameters:
    u (string): filename of input SEGY traces
    p (string): filename of output SEGY traces, same geometry as
    input traces (u)
    u0rms (string): filename of auxiliary length-1 traces containing
    scale factors for preconditioning. Same number of traces as input
    (u) and output (p)
    precond (int): preconditioning flag. 0 = no scaling, else scaling
    applied
    alpha (float): penalty weight
    '''
    try:
        #print('compawipen: u=' + u + ' p=' + p + ' u0rms=' + u0rms)
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
        if linalg.sanity(u,'su') and linalg.sanity(p,'su'):
            cmd = os.path.join(TRIP,'iwave/trace/main/awipen.x')
            if (u0rms is not None):
                if linalg.sanity(u0rms,'su'):
                    ret = os.system(cmd + ' in=' + u + ' rms=' + u0rms
                                        + ' out=' + p + ' precond=1'
                                        + ' alpha=' + str(alpha))
                else:
                    raise Exception('u0rms does not have suffix .su')
            else:
                #print(cmd + ' in=' + u  
                #            + ' out=' + p + ' precond=0 alpha='
                #            + str(alpha))
                ret = os.system(cmd + ' in=' + u  
                                    + ' out=' + p + ' precond=0'
                                    + ' alpha=' + str(alpha))
            if ret != 0:
                raise Exception('iwave/trace/main/awipen.x failed')
        else:
            raise Exception('at least one input does not have .su suffix')
    except Exception as ex:
        print(ex)
        raise Exception('called from compawipensol')

### as of 2023.07.19 this class is deprecated in favor of awipensol

class awipen(vcl.LinearOperator):
    '''
    AWI penalty operator, including AWI preconditioning and penalty weight.
    
    Constructor parameters:
        dom (segyvc:Space): domain, extended source space
        d (vcl.Vector): observed data vector in segyvc.Space
        p (vcl.Vector): predicted data vector in same sevyvc.Space as d
        precond (int): preconditioning flag - 0 for none, nonzero for yes
        alpha (float): penalty parameter
        kmax (int): max number of CG iterations for estimating u0
        eps (float): residual tolerance for CG
        rho (float): normal residual tolerance fof CG
        verbose (int): CG verbosity flag
    '''

    def __init__(self, dom, alpha=0.0, precond=0, d=None, p=None, 
                     kmax=0, eps=0.0, rho=0.0, verbose=0):
        try:
            self.dom = dom
            self.d = d
            self.p = p
            self.precond = precond
            self.alpha = alpha
            self.kmax = kmax
            self.eps = eps
            self.rho = rho
            self.verbose = verbose
            self.u0rms = None
            ### preconditioning setup
            if precond != 0:
                datapath = os.getenv('DATAPATH')
                if not os.path.exists(datapath):
                    raise Exception('Error: datapath = ' + datapath + ' not valid path')
                temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.su')
                temp.close()
                self.u0rms = temp.name
                # print('in awipen constructor: u0rms=' + self.u0rms)
                os.system('touch ' + self.u0rms)

                if not isinstance(d,vcl.Vector):
                    raise Exception('input observed data not vector')
                if not isinstance(d.space,segyvc.Space):
                    raise Exception('input observed data not SEGY')

                if not isinstance(p,vcl.Vector):
                    raise Exception('input predicted data not vector')
                if not isinstance(p.space,segyvc.Space):
                    raise Exception('input predicted data not SEGY')

                if d.space != p.space:
                    raise Exception('observed, predicted data spaces differ')
                
                Sm0 = segyvc.ConvolutionOperator(dom=self.dom,
                                            rng=self.d.space,
                                            green=self.p.data)
                # print('in awipen constructor: compute u0')
                u0=vcl.Vector(self.dom)
                vcalg.conjgrad(x=u0, b=self.d, A=Sm0, kmax=self.kmax, \
                eps=self.eps, rho=self.rho, verbose=self.verbose)
                setrms(u0.data,self.u0rms)

            self._unlink = os.unlink
            
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen constructor')

    def __del__(self):
        if self.u0rms is not None:
            self._unlink(self.u0rms)

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.dom

    def applyFwd(self,dx,dy):
        try:
            compawipen(dx.data, dy.data, self.u0rms, self.precond, self.alpha)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyFwd')
        else:
            return dy
        
    def applyAdj(self,dx,dy):
        try:
            compawipen(dx.data, dy.data, self.u0rms, self.precond, self.alpha)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyAdj')
        else:
            return dy

    def myNameIs(self):
        print('AWI penalty operator')
        print('    domain = range:')
        self.dom.myNameIs()
        print('    observed data:')
        self.d.myNameIs()
        print('    predicted data:')
        self.p.myNameIs()
        print('    preconditioning flag = ' + str(self.precond))
        print('    penalty weight = ' + str(self.alpha))
        print('    max CG iterations = ' + str(self.kmax))
        print('    relative residual tolerance = ' + str(self.eps))
        print('    relative normal residual tolerance = ' + str(self.rho))
        print('    CG verbosity flag = ' + str(self.verbose))
        print('    filename for filter trace rms = ' + self.u0rms)

class awipensol(vcl.LinearOperator):
    '''
    AWI penalty operator, including AWI preconditioning and penalty weight.
    
    Constructor parameters:
        dom (segyvc:Space): domain, extended source (awi adaptive kernel) space
        p (vcl.Vector): predicted data vector in same sevyvc.Space as d
        alpha (float): penalty parameter
        solver (vcl.LSSolver): LS solver
        d (vcl.Vector): observed data vector in segyvc.Space

        solver and d needed only for precond branch
    '''

    def __init__(self, dom, alpha=0.0, solver=None, d=None, p=None):
        try:
            # d and solver retained only for self-doc
            self.dom = dom
            self.solver = solver
            self.alpha = alpha
            self.d = d
            self.p = p
            self.u0rms = None
            self._unlink = os.unlink
                        
            ### preconditioning setup
            if (self.d is not None) and (self.p is not None) and (self.solver is not None):
                datapath = os.getenv('DATAPATH')
                if not os.path.exists(datapath):
                    raise Exception('Error: datapath = ' + datapath + ' not valid path')
                temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.su')
                temp.close()
                self.u0rms = temp.name
                # print('in awipen constructor: u0rms=' + self.u0rms)
                os.system('touch ' + self.u0rms)

                if not isinstance(d,vcl.Vector):
                    raise Exception('input observed data not vector')
                if not isinstance(d.space,segyvc.Space):
                    raise Exception('input observed data not SEGY')

                if not isinstance(p,vcl.Vector):
                    raise Exception('input predicted data not vector')
                if not isinstance(p.space,segyvc.Space):
                    raise Exception('input predicted data not SEGY')

                if d.space != p.space:
                    raise Exception('observed, predicted data spaces differ')
                
                Sm0 = segyvc.ConvolutionOperator(dom=self.dom,
                                            rng=self.p.space,
                                            green=self.p.data)
                # print('in awipen constructor: compute u0')
                u0=vcl.Vector(self.dom)
                [u0, e] = solver.solve(Sm0,d)
                setrms(u0.data,self.u0rms)
            
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen constructor')

    def __del__(self):
        if self.u0rms is not None:
            self._unlink(self.u0rms)

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.dom

    def applyFwd(self,dx,dy):
        try:
            compawipensol(dx.data, dy.data, self.alpha, self.u0rms)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyFwd')
        else:
            return dy
        
    def applyAdj(self,dx,dy):
        try:
            compawipensol(dx.data, dy.data, self.alpha, self.u0rms)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyAdj')
        else:
            return dy

    def myNameIs(self):
        print('AWI penalty operator')
        print('    domain = range:')
        self.dom.myNameIs()
        print('    predicted data:')
        self.p.myNameIs()
        print('    penalty weight = ' + str(self.alpha))
        if self.d is not None:
            print('    observed data:')
            self.d.myNameIs()
        if self.solver is not None:
            print('    LS solver:')
            self.solver.myNameIs()
        if self.u0rms is not None:
            print('    filename for filter trace rms = ' + self.u0rms)

class awiop(vcl.LinearOperator):
    '''
    Penalty-AWI operator based on FD solution of wave equation and
    scale-by-t penalty with optional preconditioning by zero-penalty
    solution norms.

    Parameters:
    F (vcl.Function): simulation operator
    d (vcl.Vector): trace data
    awisp (vcl.Space): AWI adaptive kernel space
    awisol (vcl.LSSolver): computes adaptive kernel
    data (vcl.Vector): optional data trace, flags preconditioning

    if data is provided, then it is used to compute the preconditioned
    verstion of the AWI penalty operator. Note that in the application,
    this vector should be the same as the data vector of the inverse
    problem.
    '''

    def __init__(self, dom, rng, predicted, alpha, awisol=None, observed=None):
        try:
            # range is product (data space x adaptive kernel space)
            if not isinstance(rng, vcl.ProductSpace):
                raise Exception('range not product space')
            if dom != rng[1]:
                raise Exception('domain not same as range[1]')
            if predicted.space != rng[0]:
                raise Exception('p not in 0th component of range')
            if observed is not None:
                if observed.space != rng[0]:
                    raise Exception('data not in 0th component of range')
            self.rng = rng
            # p is predicted data
            self.p = predicted
            self.alpha = alpha
            self.sol = awisol
            self.d = observed
            self.Sm0 = segyvc.ConvolutionOperator(self.rng[1], self.rng[0], self.p.data)
            self.alphaT = awipensol(self.rng[1], self.alpha,
                                            solver=self.sol, d=self.d, p=self.p)
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop constructor')

    def getDomain(self):
        return self.rng[1]

    def getRange(self):
        return self.rng

    def applyFwd(self, dx, dy):
        try:
            self.Sm0.applyFwd(dx,dy[0])
            self.alphaT.applyFwd(dx,dy[1])
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop.applyFwd')

    def applyAdj(self, dx, dy):
        try:
            cy = vcl.Vector(self.getDomain())
            self.Sm0.applyAdj(dx[0],cy)
            self.alphaT.applyAdj(dx[1],dy)
            dy.linComb(1.0,cy)
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop.applyAdj')

    def myNameIs(self):
        print('AWI operator')
        print('domain = ')
        self.getDomain().myNameIs()
        print('range = data space OPLUS domain')
        if self.d is not None:
            print('observed data:')
            self.d.myNameIs()
        print('predicted data:')
        self.p.myNameIs()
        print('penalty weight = ' + str(self.alpha))
        print('internal solver:')
        self.sol.myNameIs()

class awisep(vpm.SepFunction):
    '''
    separable function expressing AWI configuration

    domain = bulk moduli space x adaptive kernel space
    range  = data space x adaptive kernel space

    '''
    
    def __init__(self, dom, rng, F, alpha, awisol=None, observed=None):
        try:
            # range is product (data space x adaptive kernel space)
            if not isinstance(dom, vcl.ProductSpace):
                raise Exception('domain not product space')
            if dom[0] != F.getDomain():
                raise Exception('0th comp of domain not = simulator domain')
            if not isinstance(rng, vcl.ProductSpace):
                raise Exception('range not product space')
            if rng[0] != F.getRange():
                raise Exception('0th comp of range not = simulator range')
            if dom[1] != rng[1]:
                raise Exception('domain[1] not same as range[1]')
            self.dom = dom
            self.rng = rng
            self.F = F
            self.alpha = alpha
            self.awisol = awisol
            self.observed = observed

        except Exception as ex:
            print(ex)
            raise Exception('called from awisep constructor')
        
    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng
    
    def opfcn(self, m):
        try:
            return awiop(self.dom[1], self.rng, self.F(m), self.alpha,
                             awisol=self.awisol, observed=self.observed)
        except Exception as ex:
            print(ex)
            raise Exception('called from awisep.opfcn')
        
    def derfcn(self, x0, x1):
        raise Exception('awisep.defcn not defined yet (2023.07.19)')

    def myNameIs(self):
        print('AWI separable function object')
        
    





