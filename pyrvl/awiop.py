import vcl
import segyvc
import awi

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

    def __init__(self, p, awisp, alpha, awisol=None, data=None):
        # p is predicted data
        self.p = p
        # range is product (range of F x awisp)
        self.rng = vcl.ProductSpace([p.space, awisp])
        self.sp = awisp
        self.alpha = alpha
        self.sol = awisol
        self.d = data
        self.Sm0 = segyvc.ConvolutionOperator(self.sp, self.p.space, self.p.data)
        self.alphaT = awi.awipensol(self.sp, self.alpha, solver=self.sol, d=data, p=self.p)

    def getDomain(self):
        return self.sp

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
        self.dom.myNameIs()
        print('range = data space OPLUS domain')
        if self.d is not None:
            print('observed data:')
            self.d.myNameIs()
        print('predicted data:')
        self.p.myNameIs()
        print('penalty weight = ' + str(self.alpha))
        print('internal solver:')
        self.sol.myNameIs()

