import vcl
import npvc

def f0(xxx,yyy):
    try:
        yyy.linComb(1.0,xxx)
    except Exception as ex:
        print(ex)
        raise Exception('called from f1')
    
def f1(x,y):
    try:
        x.linComb(1.0,y)
    except Exception as ex:
        print(ex)
        raise Exception('called from f1')

def f2(x,y):
    try:
        print(x.dot(y))
    except Exception as ex:
        print(ex)
        raise Exception('called from f2')

def f3(dom,y):
    try:
        dom.printData(y.data)
    except Exception as ex:
        print(ex)
        raise Exception('called from f3')

def f4(x,d):
    try:
        x.link(d)
    except Exception as ex:
        print(ex)
        raise Exception('called from f4')

def f5(sp1,sp2):
    spl=[]
    spl.append(sp1)
    spl.append(sp2)
    ps=vcl.ProductSpace(spl)
    ps.myNameIs()

def f6(sp1,sp2):
    spl=[]
    spl.append(sp1)
    spl.append(sp2)
    ps=vcl.ProductSpace(spl)
    dl=[]
    dl.append(sp1.getData())
    dl.append(sp2.getData())
    print('return f6')
    print(ps.isData(dl))

def f7(sp1,sp2):
    spl=[]
    spl.append(sp1)
    spl.append(sp2)
    ps=vcl.ProductSpace(spl)
    dl=[]
    dl.append(sp1.getData())
    dl.append(sp1.getData())
    print('return f7')
    print(ps.isData(dl))

def f8(sp1,sp2):
    f = npvc.OpExpl1(sp1,sp2)
    x = vcl.Vector(sp1)
    y = vcl.Vector(sp2)
    x.data[0]=1
    x.data[1]=-2
    print('input vector:')
    x.myNameIs()
    f.apply(x,y)
    print('output of apply method:')
    y.myNameIs()
    dfx = f.deriv(x)
    print('output of deriv method:')
    dfx.myNameIs()
    dx=vcl.Vector(sp1)
    dx.data[0]=2
    dx.data[1]=-3
    dy=vcl.Vector(sp2)
    print('input to deriv.applyFwd')
    dx.myNameIs()
    dfx.applyFwd(dx,dy)
    print('output of deriv.applyFwd')
    dy.myNameIs()

def f9(sp1,sp2):
    f = npvc.OpExpl2(sp1,sp2)
    x = vcl.Vector(sp1)
    yy = vcl.Vector(sp2)
    x.data[0][0]=1
    x.data[1][0]=-2
    print('input vector:')
    x.myNameIs()
    f.apply(x,yy)
    print('output of apply method:')
    yy.myNameIs()
    dfx = f.deriv(x)
    print('output of deriv method:')
    dfx.myNameIs()
    dx=vcl.Vector(sp1)
    dx.data[0][0]=2
    dx.data[1][0]=-3
    dy=vcl.Vector(sp2)
    print('input to deriv.applyFwd')
    dx.myNameIs()
    dfx.applyFwd(dx,dy)
    print('output of deriv.applyFwd')
    dy.myNameIs()
    dfx0 = f.partialDeriv(x,0)
    cdx0 = dx.component(0)
    dy0=vcl.Vector(sp2)
    dfx0.applyFwd(cdx0,dy0)
    print('input of partial deriv 0')
    cdx0.myNameIs()
    print('output of partial deriv 0')
    dy0.myNameIs()
    dfx1 = f.partialDeriv(x,1)
    cdx1 = dx.component(1)
    dy1=vcl.Vector(sp2)
    dfx1.applyFwd(cdx1,dy1)
    print('input of partial deriv 1')
    cdx1.myNameIs()
    print('output of partial deriv 1')
    dy1.myNameIs()
    
dom=npvc.Space(2)
dom1=npvc.Space(1)
dom2=npvc.Space(1)
doml=[]
doml.append(dom1)
doml.append(dom2)
pdom=vcl.ProductSpace(doml)
rng=npvc.Space(3)
x=vcl.Vector(dom)
x.data[0]=1
x.data[1]=2
z=vcl.Vector(dom)
z.data[0]=3
z.data[1]=4
w=vcl.Vector(dom)
y=vcl.Vector(rng)

print()
try:
    f0(x,z)
except Exception as ex:
    print(ex)
    
print()
try:
    f1(x,y)
except Exception as ex:
    print(ex)

print()
try:
    f2(x,y)
except Exception as ex:
    print(ex)
    
print()
try:
    f3(dom,y)
except Exception as ex:
    print(ex)

print()
try:
    f4(x,z.data)
except Exception as ex:
    print(ex)
else:
    print('successful link')

print()
try:
    f4(x,y.data)
except Exception as ex:
    print(ex)
else:
    print('successful link')

print()
try:
    f4(x,w.data)
except Exception as ex:
    print(ex)
else:
    print('successful link')

print()
try:
    f5(dom,rng)
except Exception as ex:
    print(ex)

print()
try:
    f6(dom,rng)
except Exception as ex:
    print(ex)

print()
try:
    f7(dom,rng)
except Exception as ex:
    print(ex)

print()
try:
    f8(dom,rng)
except Exception as ex:
    print(ex)

print()
try:
    f9(pdom,rng)
except Exception as ex:
    print(ex)
