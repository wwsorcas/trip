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


#def compawipensol(u, p, alpha=1.0, u0rms=None):
def applytxgain(u, p, u0rms=None, spower=0, tpower=1):
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
    '''
    try:
#        print('compawipen')
#        print('next')
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
# 23.11.15 old version        
#        print('sanity')
#        if linalg.sanity(u,'su') and linalg.sanity(p,'su'):
#            cmd = os.path.join(TRIP,'iwave/trace/main/awipen.x') \
#                + ' in=' + u + ' rms=' + u0rms \
#                + ' out=' + p + ' precond=1' \
#                + ' alpha=' + str(alpha)
#            if (u0rms is not None):
#                if linalg.sanity(u0rms,'su'):
#                    print('alpha=' + str(alpha))
#                    ret = os.system(cmd)
#                else:
#                    raise Exception('u0rms does not have suffix .su')
#            else:
#                #print(cmd + ' in=' + u  
#                #            + ' out=' + p + ' precond=0 alpha='
#                #            + str(alpha))
#                ret = os.system(cmd + ' in=' + u  
#                                    + ' out=' + p + ' precond=0'
#                                    + ' alpha=' + str(alpha))
#
# 23.11.15 new version
        if linalg.sanity(u,'su') and linalg.sanity(p,'su'):
            if u0rms is None:
                cmd = os.path.join(TRIP,'iwave/trace/main/txgain.x') \
                + ' in=' + u + ' out=' + p + ' spower=0 tpower=' \
                + str(tpower)
            else:
                if not linalg.sanity(u0rms,'su'):
                    raise Exception('u0rms does not have suffix .su')
                cmd = os.path.join(TRIP,'iwave/trace/main/txgain.x') \
                + ' in=' + u + ' out=' + p + ' sfile=' + u0rms \
                + ' spower=' + str(spower) + ' tpower=' + str(tpower)
            ret = os.system(cmd)
            if ret != 0:
                raise Exception('iwave/trace/main/txgain.x failed')
        else:
            raise Exception('at least one input does not have .su suffix')
    except Exception as ex:
        print(ex)
        raise Exception('called from applytxgain')

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
#            compawipensol(dx.data, dy.data, u0rms=self.u0rms)
            if self.u0rms is None:
                applytxgain(dx.data, dy.data, spower=0, tpower=1)
            else:
                applytxgain(dx.data, dy.data, u0rms=self.u0rms,
                            spower=-1, tpower=1)

        except Exception as ex:
            print(ex)
            raise Exception('called from awipen.applyFwd')
        else:
            return dy
        
    def applyAdj(self,dx,dy):
        try:
#            compawipensol(dx.data, dy.data, u0rms=self.u0rms)
            if self.u0rms is None:
                applytxgain(dx.data, dy.data, spower=0, tpower=1)
            else:
                applytxgain(dx.data, dy.data, u0rms=self.u0rms,
                            spower=-1, tpower=1)
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
    rng (vcl.ProductSpace): 
        comp 0 = data space, comp 1 = comp 2 = adaptive filter space
    awifilt (vcl.LinearOperator): predicted data filter, dom = adaptive filter 
        space, rng = data space
    awipen (vcl.LinearOperator): penalty operator, domain = range = adaptive
        filter space
    alpha (float): penalty weight
    sigma (float): regularization weight

    if data is provided, then it is used to compute the preconditioned
    verstion of the AWI penalty operator. Note that in the application,
    this vector should be the same as the data vector of the inverse
    problem.
    '''

    def __init__(self, awifilt, awipen, alpha=None, sigma=None):
        try:
            # range is product (data space x adaptive kernel space)
            if not isinstance(rng, vcl.ProductSpace):
                raise Exception('range not product space')
            if awifilt.getDomain() != rng[1]:
                raise Exception('filter domain not same as range[1]')
            if awipen.getDomain() != rng[1]:
                raise Exception('penalty domain not same as range[1]')
            if awipen.getRange() != rng[1]:
                raise Exception('penalty range not same as range[1]')          
            if awifilt.getDomain() != rng[2]:
                raise Exception('filter domain not same as range[2]')
            if awifilt.getRange() != rng[0]:
                raise Exception('filter range not same as range[0]')
            self.rng = rng
            # p is predicted data
            self.filt = awifilt
            self.pen = awipen
            self.alpha = alpha
            self.sigma = sigma                
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop constructor')

    def getDomain(self):
        return self.rng[1]

    def getRange(self):
        return self.rng

    def filter(self):
        return self.filt

    def penalty(self):
        return self.pen

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
        print('filter = ')
        self.filt.myNameIs()
        print('penalty op = ')
        self.pen.myNameIs()
        print('range = filter range OPLUS domain OPLUS domain')
        print('penalty weight = ' + str(self.alpha))
        print('regularization weight = ' + str(self.sigma))

class awisep(vpm.SepFunction):
    '''
    separable function expressing AWI configuration

    domain = bulk moduli space x adaptive kernel space
    range  = data space x adaptive kernel space

    '''
    
    def __init__(self, dom, rng, F, alpha, sigma, awisol=None, observed=None):
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
            if dom[1] != rng[2]:
                raise Exception('domain[1] not same as range[2]')
            self.dom = dom
            self.rng = rng
            self.F = F
            self.alpha = alpha
            self.sigma = sigma
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
            ### convolution with predicted data
            ### note that the third arg in the constructor of the
            ### segyvc conv op is a segy data filename, not a vcl.Vector
            filt = segyvc.ConvolutionOperator(dom=self.dom[1],
             rng=self.rng[0], green=self.F(m).data)
            ### penalty operator, owns alpha=0 soln & trace norms
            pen = awipensol(self.dom[1], filt,
                                         sigma=self.sigma,
                                         solver=self.awisol,
                                         d=self.observed)
            op = awiop(self.dom[1], self.rng, filt, pen,
                             alpha=self.alpha, sigma=self.sigma)
            return op
        except Exception as ex:
            print(ex)
            raise Exception('called from awisep.opfcn')
        
    def derfcn(self, x0, x1):
        raise Exception('awisep.defcn not defined yet (2023.07.19)')

    def compj0(self):
        if self.filt is None:
            self.filt = segyvc.ConvolutionOperator(dom=self.dom[1],
                                    rng=self.rng[0],
                                    green=self.F(m).data)
        ### penalty operator, owns alpha=0 soln & trace norms
        if self.pen is None:
            self.pen = awipensol(self.dom[1], filt=self.filt,
                                     sigma=self.sigma,
                                     solver=self.awisol,
                                     d=self.observed)
        x = pen.innererr().norm()
        x0 = pen.innererr()[0].norm()
        x1 = pen.innererr()[1].norm()
        j0 = 0.5*x*x
        print('j0 = ' + str(j0))
        print('components of j0:')
        print('0: ' + str(0.5*x0*x0))
        print('1: ' + str(0.5*x1*x1))
        ta = pen*pen.innersol()
        wg = 0.5*ta.dot(ta)
        print('wg = ' + str(wg))
        
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
            self.txu0 = vcl.Vector(self.dom)
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
#            compawipensol(self.u0.data, self.txu0.data, self.u0rms)
            if self.u0rms is None:
                applytxgain(self.u0.data, self.txu0.data, spower=0, tpower=1)
            else:
                applytxgain(self.u0.data, self.txu0.data, u0rms=self.u0rms,
                            spower=-1, tpower=1)

            # compute value
            n = self.txu0.norm()
            self.val = 0.5*n*n

            # placeholders
            self.grad = None
            self.hess = None
            
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
           # if self.grad is None:
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

        
    





