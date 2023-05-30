import vcl
import npvc
import numpy as np
import vpm

class expl1(vcl.ScalarFunction):

    def __init__(self,sp):
        try:
            if not isinstance(sp,npvc.Space):
                raise Exception('input not npvc.Space')
            self.n = sp.dim
            self.dom = sp
            self.mat0 = np.zeros((2*n,n))
            for i in range(n):
                self.mat0[i][i] = 1.0
            self.rhs = np.zeros((2*n,1))
            for i in range(n):
                self.rhs[i]=(-i-1)**(i+1)
        except Exception as ex:
            print(ex)
            raise Exception('called from vpmtest.expl1 constructor')
        
    def getDomain(self):
        return self.dom

    def value(self,x):
        xs = x.norm()**2/(1 + x.norm()**2)
        mat = np.copy(self.mat0)
        for i in range(self.n):
            mat[n+i][i]=(i+1)*xs
        [w, res, rk, s] = np.linalg.lstsq(mat,self.rhs, rcond=None)
        #print('w = ')
        #print(w)
        r = mat @ w - self.rhs
        #print('r = ')
        #print(r)
        #print('res = ')
        #print(res)
        return 0.5*np.dot(r.T,r)[0][0]

    def raw_gradient(self,x):
        xs = x.norm()**2/(1+x.norm()**2)
        mat = np.copy(self.mat0)
        for i in range(self.n):
            mat[n+i][i]=(i+1)*xs
        [w, res, rk, s] = np.linalg.lstsq(mat,self.rhs, rcond=None)
        r = mat @ w - self.rhs
        g = vcl.Vector(self.getDomain())
        p = 0.0
        for i in range(self.n):
            p = p + (i+1)*w[i]*r[n+i]
        np.copyto(g.data, 2.0*x.data*p/((1+x.norm()**2)**2))
        return g

    def raw_Hessian(self,x):
        return None

    def myNameIs(self):
        print('vpm expl 1, n = ' + str(n))

###################

class sepexpl1(vpm.SepFunction):

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
            self.mat0 = np.zeros((2*n,n))
            self.mat1 = np.zeros((2*n,n))
            for i in range(n):
                self.mat0[i][i] = 1.0
                self.mat1[n+i][i] = i+1
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


######################

n=4
sp0 = npvc.Space(n)
sp1 = npvc.Space(n)
sp2 = npvc.Space(2*n)

sp = vcl.ProductSpace([sp0, sp1])

f = expl1(sp0)

x = vcl.Vector(f.getDomain())

print('value at origin = ' + str(f(x)))

x.data[0] = 1.0
x.data[1] = 0.2
x.data[2] =-0.5
x.data[3] =-0.4

fx = vcl.StandardJet(f,x)

print('value at [1.0,0.2,-0.5,-0.4] = ' + str(fx.value()))
print('gradient at [1.0,0.2,-0.5,-0.4] = ')
print(fx.gradient().data)

x.data[0] = 1.0
x.data[1] = 0.0
x.data[2] = 0.0
x.data[3] = 0.0

fx = vcl.StandardJet(f,x)

print('value at [1.0,0.0,0.0,0.0] = ' + str(fx.value()))
print('gradient at [1.0,0.0,0.0,0.0] = ')
print(fx.gradient().data)

for i in range(11):
    x.data[0]=i*0.1
    fx = vcl.StandardJet(f,x)
    print('value at x[0]=' + str(x.data[0]) + ' = ' + str(fx.value()))

#############

print('\nvpmcgjet evaluation:')

g = sepexpl1(sp,sp2)

kmax = 10
eps  = 1.e-3
rho  = 1.e-3
b = vcl.Vector(g.getRange())
# nonzeros in b only in top half
for i in range(n):
    b.data[i]=(-i-1)**(i+1)

gx = vpm.vpmcgjet(x, g, b, kmax, eps, rho, verbose=2)

print('\nvpmcgjet value = ' + str(gx.value()))
print('\ngradient = ')
print(gx.gradient().data)
