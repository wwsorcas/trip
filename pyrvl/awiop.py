class awiop(vcl.LinearOperator):

    def __init__(self, dom, rng, m, buoy=None, w=None, alpha=0.0,
                     precond=0, d=None,
                     kmax=None, eps=None, rho=None, verbose=None,
                     order=2, sampord=1, nsnaps=20,
                     cfl=0.5, cmin=1.0, cmax=3.0, dmin=0.8, dmax=3.0,
                     nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0):
        '''
        AWI linear operator based on point-source acoustic simulation with
        known and fixed buoyancy field. Preconditioning optional.
        Parameters:
        dom (vcl.Space): extended source space
        rng (vcl.ProductSpace): first factor = data space, second factor=
        extended source space
        m (vcl.Vector): bulk modulus
        buoy (string): rsf filename for buoyancy
        w (string): su filename for source wavelet 
        d (vcl.Vector): data, member of rng[0]
        alpha (float): penalty weight
        precond (int): preconditioning flag. if set, next 5 args must have 
        sensible values        
        
        '''
        try:
            if not isinstance(rng, vcl.ProductSpace):
                raise Exception('range not product space')
            if rng[1] != dom:
                raise Exception('range[1] not same as domain')
            if buoy is None:
                raise Exception('buoyancy field must be specified')
            if not isinstance(buoy,string):
                raise Exception('bouyancy filename must be specified')
            if w is None:
                raise Exception('source wavelet must be specified')
            if not isinstance(buoy,string):
                raise Exception('source waveelt filename must be specified')
            asgop = asg.fsbop(m.space, rng[0], buoyancy=buoy, source_p=w,
                            order=order, sampord=sampord, nsnaps=nsnaps,
                            cfl=cfl, cmin=cmin, cmax=cmax, dmin=dmin, dmax=dmax,
                            nl1=nl1, nr1=nr1, nl2=nl2, nr2=nr2, pmlampl=pmlampl)
            self.dom = dom
            self.rng = rng
            self.p = asgop(m)
            self.d = d
            self.p = 
            self.precond = precond
            self.alpha = alpha
            self.kmax = kmax
            self.eps = eps
            self.rho = rho
            self.verbose = verbose

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
m
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

