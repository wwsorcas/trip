import linalg
import vcl
import vcalg
import segyvc
import rsfvc
import os
import tempfile

TRIP = os.getenv('TRIP')

def setrms(u, urms):
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

def compawipen(u, p, u0rms, precond, alpha):
    try:
        #print('compawipen: u=' + u + ' p=' + p + ' u0rms=' + u0rms)
        if not isinstance(TRIP,str):
            raise Exception('string TRIP not defined')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
        if linalg.sanity(u,'su') and linalg.sanity(u0rms,'su') and linalg.sanity(p,'su'):
            cmd = os.path.join(TRIP,'iwave/trace/main/awipen.x')
            ret = os.system(cmd + ' in=' + u + ' rms=' + u0rms + ' out=' + p +
                                ' precond=' + str(precond) + ' alpha=' + str(alpha))
            if ret != 0:
                raise Exception('iwave/trace/main/awipen.x failed')
        else:
            raise Exception('at least one input does not have .su suffix')
    except Exception as ex:
        print(ex)
        raise Exception('called from compawipen')

class awipen(vcl.LinearOperator):
    '''
    AWI penalty operator, including AWI preconditioning and penalty weight.
    
    Constructor parameters:
        dom (segyvc:Space): domain, extended source space
        d (vcl.Vector): observed data vector in segyvc.Space
        p (vcl.Vector): predicted data vector in same sevyvc.Space as d
        precond (int): preconditioning flag - 0 for none, nonzero for yes
        alpha (float): penalty parametere 
        kmax (int): max number of CG iterations for estimating u0
        eps (float): residual tolerance for CG
        rho (float): normal residual tolerance fof CG
        verbose (int): CG verbosity flag
    '''

    def __init__(self, dom, d, p, precond, alpha, kmax, eps, rho, verbose):
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
            datapath = os.getenv('DATAPATH')
            if not os.path.exists(datapath):
                raise Exception('Error: datapath = ' + datapath + ' not valid path')
            temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.su')
            temp.close()
            self.u0rms = temp.name
            print('in awipen constructor: u0rms=' + self.u0rms)
            os.system('touch ' + self.u0rms)
            
            if d.space != p.space:
                raise Exception('observed, predicted data spaces differ')
            Sm0 = segyvc.ConvolutionOperator(dom=self.dom, rng=self.d.space, \
                                green=self.p.data)
            print('in awipen constructor: compute u0')
            u0=vcl.Vector(self.dom)
            vcalg.conjgrad(x=u0, b=self.d, A=Sm0, kmax=self.kmax, \
               eps=self.eps, rho=self.rho, verbose=self.verbose)
            setrms(u0.data,self.u0rms)

            self._unlink = os.unlink
            
        except Exception as ex:
            print(ex)
            raise Exception('called from awipen constructor')

    def __del__(self):
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

class awiop(vcl.LinearOperator):

    def __init__(self, dom, d, p, precond, alpha, kmax, eps, rho, verbose):
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

            self.rng = vcl.ProductSpace([self.d.space, self.dom])

            self.Sm0 = segyvc.ConvolutionOperator(self.dom, self.d.space, self.p.data)
            self.alphaT = awipen(dom, d, p, precond, alpha, kmax, eps, rho, verbose)

        except Exception as ex:
            print(ex)
            raise Exception('called from awiop constructor')

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def applyFwd(self,dx,dy):
        try:
            #dy[0]=self.Sm0*dx
            #dy[1]=self.alphaT*dx
            self.Sm0.applyFwd(dx,dy[0])
            self.alphaT.applyFwd(dx,dy[1])
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop.applyFwd')
        else:
            return dy
        
    def applyAdj(self,dx,dy):
        try:
#            print('AWI Op: applyAdj')
#            print('transp comp 0')
#            print('|dx[0]|=' + str(dx[0].norm()))
            # cy = vcl.transp(self.Sm0)*dx[0]
            cy = vcl.Vector(self.getDomain())
            self.Sm0.applyAdj(dx[0],cy)
#            print('|cy|=' + str(cy.norm()))
#            print('transp comp 1')
#            print('|dx[1]|=' + str(dx[1].norm()))
#            #dy = vcl.transp(self.alphaT)*dx[1]
            self.alphaT.applyAdj(dx[1],dy)
#            print('|dy|=' + str(dy.norm()))
#            print('lincomb')
            dy.linComb(1.0,cy)
#            print('|dy|=' + str(dy.norm()))
#            dy.myNameIs()
#            print('exit AWI Op: applyAdj')
        except Exception as ex:
            print(ex)
            raise Exception('called from awiop.applyAdj')
        else:
            return dy

    def myNameIs(self):
        print('AWI operator')
        print('domain = ')
        self.dom.myNameIs()
        print('range = data space OPLUS domain')
        print('observed data:')
        self.d.myNameIs()
        print('predicted data:')
        self.p.myNameIs()
        print('preconditioning flag = ' + str(self.precond))
        print('penalty weight = ' + str(self.alpha))
        print('max CG iterations = ' + str(self.kmax))
        print('relative residual tolerance = ' + str(self.eps))
        print('relative normal residual tolerance = ' + str(self.rho))
        print('CG verbosity flag = ' + str(self.verbose))

