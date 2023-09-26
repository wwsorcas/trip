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


def compawipensol(u, p, alpha=1.0, u0rms=None):
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
#        print('compawipen')
#        print('compawipen: u=' + u + ' p=' + p + ' u0rms=' + u0rms + ' alpha=' + str(alpha))
#        print('next')
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
#        print('sanity')
        if linalg.sanity(u,'su') and linalg.sanity(p,'su'):
            cmd = os.path.join(TRIP,'iwave/trace/main/awipen.x') \
                + ' in=' + u + ' rms=' + u0rms \
                + ' out=' + p + ' precond=1' \
                + ' alpha=' + str(alpha)
            if (u0rms is not None):
                if linalg.sanity(u0rms,'su'):
#                    print('alpha=' + str(alpha))
                    ret = os.system(cmd)
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

class awipensol(vcl.LinearOperator):
    '''
    AWI penalty operator, including AWI preconditioning and penalty weight.
    
    Constructor parameters:
        dom (segyvc:Space): domain, extended source (awi adaptive kernel) space
        p (vcl.Vector): predicted data vector in same sevyvc.Space as d
        sigma (float): reg parameter
        solver (vcl.LSSolver): LS solver
        d (vcl.Vector): observed data vector in segyvc.Space

        solver and d needed only for precond branch
    '''

    def __init__(self, dom, filt=None, sigma=None, solver=None, d=None):
        try:
            # d and solver retained only for self-doc
            self.dom = dom
            self.filt = filt
            self.solver = solver
            self.sigma = sigma
            self.d = d
            self.u0 = None
            self.e0 = None
            self.u0rms = None
            self._unlink = os.unlink
                        
            ### preconditioning setup
            if (self.d is not None) and (self.solver is not None) and (self.filt is not None):
                if dom != filt.getDomain():
                    raise Exception('domain not = domain of filter op')
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

                if d.space != filt.getRange():
                    raise Exception('observed data space != filter range')

                print('\nawi.awipensol: compute adaptive filter and preconditioner data')
                print('sigma = ' + str(self.sigma))
                if self.sigma is None:
                    rfilt = filt
                    dp = d
                else:
                    rfilt = vcl.LinearOpReg(filt, self.sigma)
                    dp = vcl.Vector(rfilt.getRange())
                    dp[0].copy(d)
                # print('in awipen constructor: compute u0')
                self.u0=vcl.Vector(rfilt.getDomain())
                self.e0=vcl.Vector(rfilt.getRange())
                [tmpu0, tmpe0] = solver.solve(rfilt,dp)
                self.u0.copy(tmpu0)
                self.e0.copy(tmpe0)
                setrms(self.u0.data,self.u0rms)
            
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

    def innersol(self):
        try:
            if self.u0 is None:
                raise Exception('no precond, innersol not computed')
            return self.u0
        except Exception as ex:
            print(ex)
            raise Exception('called from awi.awipensol.innersol')
        
    def innerrms(self):
        try:
            if self.u0rms is None:
                raise Exception('no precond, innerrms not computed')
            return self.u0rms
        except Exception as ex:
            print(ex)
            raise Exception('called from awi.awipensol.innerrms')

    def innererr(self):
        try:
            if self.e0 is None:
                raise Exception('no precond, innererr not computed')
            return self.e0
        except Exception as ex:
            print(ex)
            raise Exception('called from awi.awipensol.innererr')
        
    def applyFwd(self,dx,dy):
        try:
            compawipensol(dx.data, dy.data, u0rms=self.u0rms)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyFwd')
        else:
            return dy
        
    def applyAdj(self,dx,dy):
        try:
            compawipensol(dx.data, dy.data, u0rms=self.u0rms)
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyAdj')
        else:
            return dy

    def myNameIs(self):
        print('AWI penalty operator')
        print('    domain = range:')
        self.dom.myNameIs()
        print('    filter:')
        self.p.myNameIs()
        print('    reg weight = ' + str(self.sigma))
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
    dom (vcl.Space): adaptive filter space
    rng (vcl.ProductSpace): 
        comp 0 = data space, comp 1 = comp 2 = adaptive filter space
    predicted (vcl.Vector): predicted data
    alpha (float): penalty weight
    sigma (float): reg weight
    awisol (vcl.LSSolver): computes adaptive kernel
    observed (vcl.Vector): optional data trace, flags preconditioning

    if data is provided, then it is used to compute the preconditioned
    verstion of the AWI penalty operator. Note that in the application,
    this vector should be the same as the data vector of the inverse
    problem.
    '''

    def __init__(self, dom, rng, awifilt, awipen, alpha=None, sigma=None, observed=None):
        try:
            # range is product (data space x adaptive kernel space)
            if not isinstance(rng, vcl.ProductSpace):
                raise Exception('range not product space')
            if dom != rng[1]:
                raise Exception('domain not same as range[1]')
            if dom != rng[2]:
                raise Exception('domain not same as range[2]')
            if awifilt.getRange() != rng[0]:
                raise Exception('p not in 0th component of range')
            if observed is not None:
                if observed.space != rng[0]:
                    raise Exception('data not in 0th component of range')
            self.rng = rng
            # p is predicted data
            self.filt = awifilt
            self.pen = awipen
            self.alpha = alpha
            self.sigma = sigma                
            self.d = observed
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop constructor')

    def getDomain(self):
        return self.rng[1]

    def getRange(self):
        return self.rng

    def innersol(self):
        return self.pen.innersol()

    def innerrms(self):
        return self.pen.innerrms()

    def applyFwd(self, dx, dy):
        try:
            self.filt.applyFwd(dx,dy[0])
            if self.alpha is None:
                dy[1].scale(0.0)
            else:
                self.pen.applyFwd(dx,dy[1])
                dy[1].scale(self.alpha)
            if self.sigma is None:
                dy[2].scale(0.0)
            else:
                dy[2].copy(dx)
                dy[2].scale(self.sigma)
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop.applyFwd')

    def applyAdj(self, dx, dy):
        try:
            cy = vcl.Vector(self.getDomain())
            self.filt.applyAdj(dx[0],dy)
            if self.alpha is not None:
                self.pen.applyAdj(dx[1],cy)
                dy.linComb(self.alpha,cy)
            if self.sigma is not None:
                dy.linComb(self.sigma,dx[2])
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

class awiwg(vcl.ScalarJet):
    '''
    Original Warner-Guasch AWI function jet, including preconditioning.
    
    Constructor parameters:
        dom (segyvc:Space): domain, extended source (awi adaptive kernel) space
        sim (vcl.Function): simulator / modeling operator
        mod (vcl.Vector): model, in domain of sim
        solver (vcl.LSSolver): LS solver
        d (vcl.Vector): observed data vector in range of sim
        precond (int): if not = 0, then precondition per AWI
=    '''

    def __init__(self, dom, sim, mod, solver, d, precond=1):
        try:
            
            # sanity checks
            if not isinstance(d,vcl.Vector):
                raise Exception('input observed data not vector')
            if not isinstance(d.space,segyvc.Space):
                raise Exception('input observed data not SEGY')
            if not isinstance(mod,vcl.Vector):
                raise Exception('input model not vector')
            if mod.space != sim.getDomain():
                raise Exception('mod not in domain of sim')
            if d.space != sim.getRange():
                raise Exception('data vector d not in range of sim')
            # store instance data
            self.dom = dom
            self.sim = sim
            self.mod = mod
            self.solver = solver
            self.d = d
            self.precond = precond
            # minimizer of |Su-d|
            self.u0 = vcl.Vector(self.dom)
            # after application of AWI penalty
            self.au0 = vcl.Vector(self.dom)
            # file of u0 trace rms values 
            self.u0rms = None
            self._unlink = os.unlink
            
            # set up tmp filename for u0rms
            if precond != 0:
                datapath = os.getenv('DATAPATH')
                if not os.path.exists(datapath):
                    raise Exception('Error: datapath = ' + datapath + ' not valid path')
                temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.su')
                temp.close()
                self.u0rms = temp.name
                # print('in awipen constructor: u0rms=' + self.u0rms)
                os.system('touch ' + self.u0rms)
                
            # predicted data
            self.p = sim(mod)

            # residual vector Su-d
            self.e = vcl.Vector(d.space)

            # convolution by predicted data
            Sm0 = segyvc.ConvolutionOperator(dom=self.dom,
                                                rng=self.p.space,
                                                green=self.p.data)
            # decon observed by predicted
            [self.u0, self.e] = solver.solve(Sm0,d)
            if self.u0rms is not None:
                setrms(self.u0.data,self.u0rms)

            # apply scale-by-t and optionally reciprocal trace norms
            compawipensol(self.u0.data, self.au0.data, 1.0, self.u0rms)

            # compute value
            n = self.au0.norm()
            self.val = 0.5*n*n
            
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen constructor')

    def __del__(self):
        if self.u0rms is not None:
            self._unlink(self.u0rms)

    def point(self):
        return self.mod

    def value(self):
        return self.val

    def innersol(self):
        return self.u0

    def innerrms(self):
        try:
            if self.u0rms is None:
                raise Exception('no u0rms')
            return self.u0rms
        except Exception as ex:
            print(ex)
            raise Exception('called from awipensol:innerrms')

    def gradient(self):
        try:
            raise Exception('gradient not available')
        except Exception as ex:
            print(ex)
            raise Exception('called from awiwg.gradient')
        
    def Hessian(self):
        try:
            raise Exception('Hessian not available')
        except Exception as ex:
            print(ex)
            raise Exception('called from awiwg.Hessian')

    def myNameIs(self):
        print('Warner-Guasch function')
        print('    adaptive kernel space:')
        self.dom.myNameIs()
        print('    simulator:')
        self.sim.myNameIs()
        print('    model:')
        self.mod.myNameIs()
        print('    observed data:')
        self.d.myNameIs()
        print('    LS solver:')
        self.solver.myNameIs()
        if self.u0rms is not None:
            print('    preconditioned by kernel rms')
            print('    filename for kernel trace rms = ' + self.u0rms)
        else:
            print('    no preconditioning applied')

        
    





