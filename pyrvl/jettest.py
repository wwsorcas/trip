import vcl
import npvc
import vcalg

class RosieJet(vcl.ScalarJet):

    def __init__(self,x,rhs):
        self.x = vcl.Vector(x.space)
        self.x.copy(x)
        self.f = npvc.DoubleRosie(x.space)
        self.b = rhs
        self.ls = vcl.LeastSquares(self.f,self.b)
        
    def point(self):
        return self.x

    def value(self):
        return self.ls(self.x)

    def gradient(self):
        return self.ls.gradient(self.x)

    def Hessian(self):
        return self.ls.Hessian(self.x)

    def myNameIs(self):
        print('Double Rosie Jet at x = ')
        self.x.myNameIs()
        print('and b = ')
        self.b.myNameIs()
        print('and value = ' + str(self.value()))
        
        
sp = npvc.Space(4)
x = vcl.Vector(sp)
x.data[0]=0
x.data[1]=0
x.data[2]=0
x.data[3]=0

b = vcl.Vector(sp)
b.data[0]=0
b.data[1]=-1
b.data[2]=0
b.data[3]=-1

j = RosieJet(x,rhs=b)

j.myNameIs()

print('\n')

class foo:
    def __init__(self,J,rhs):
        self.J = J
        self.rhs = rhs

    def bar(self,x):
        Jx = self.J(x,rhs = self.rhs)
        print('\njet at x:')
        Jx.myNameIs()

args = dict(rhs=b)
woo = foo(RosieJet,**args)
woo.bar(x)

f = npvc.DoubleRosie(sp)
ls = vcl.LeastSquares(f,b)

args1 = dict(fcn=ls)
jx = vcl.StandardJet(x,**args1)
print('\nStandardJet at x')
jx.myNameIs()

jy = vcl.StandardJet(x,fcn=ls)
print('\nStarndardJet at x, alt call')
jy.myNameIs()

print('\n****************** ')
print('\n Newton test')
print('\n****************** ')

x.data[0]=-1.2
x.data[1]=1.0
x.data[2]=-1.2
x.data[3]=1.0
b.data[0]=0.0
b.data[1]=-1.0
b.data[2]=0.0
b.data[3]=-1.0

vcalg.trcgnewt(x, vcl.StandardJet, newtmax=40, newteps=0.001,
             cgmax=10, cgeps=1.e-6, Delta=10.0, \
             mured=0.5, muinc=1.8, \
             gammared=0.1, gammainc=0.9, nverbose=2, cgverbose=0, \
             maxreds=0, gnorm0=None, jetargs=args1)

print('\nSolution')
x.myNameIs()

x.data[0]=-1.2
x.data[1]=1.0
x.data[2]=-1.2
x.data[3]=1.0

print('\nAnother way to go:')
vcalg.trcgnewt(x, RosieJet, newtmax=40, newteps=0.001,
             cgmax=10, cgeps=1.e-6, Delta=10.0, \
             mured=0.5, muinc=1.8, \
             gammared=0.1, gammainc=0.9, nverbose=2, cgverbose=0, \
             maxreds=0, gnorm0=None, jetargs=args)

print('\nSolution')
x.myNameIs()



    

    
