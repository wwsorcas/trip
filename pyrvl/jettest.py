import vcl
import npvc

class RosieJet(vcl.ScalarJet):

    def __init__(self,x):
        self.x = vcl.Vector(x.space)
        self.x.copy(x)
        self.b = vcl.Vector(self.x.space)
        self.b.data[0]=0
        self.b.data[1]=-1
        self.b.data[2]=0
        self.b.data[3]=-1
        self.f = npvc.DoubleRosie(x.space)
        self.ls = vcl.LeastSquares(self.f,self.b)
        
    def point(self):
        return self.x

    def value(self):
        return self.ls(x)

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

j = RosieJet(x)

j.myNameIs()

print('\n')

class foo:
    def __init__(self,J,x):
        self.J = J
        self.x = x

    def bar(self):
        Jx = self.J(self.x)
        print('jet at x:')
        Jx.myNameIs()

woo = foo(RosieJet,x)
woo.bar()
    

    
