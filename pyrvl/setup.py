import vcl
import npvc

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
    
dom=npvc.Space(2)
rng=npvc.Space(3)
x=vcl.Vector(dom)
z=vcl.Vector(dom)
w=vcl.Vector(dom)
y=vcl.Vector(rng)

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
