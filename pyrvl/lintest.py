import vcl
import npvc

def f(x,y):
    y=x+1

d=0
e=0
f(d,e)
print('x=' + str(e))

sp=npvc.Space(3)
xx=vcl.Vector(sp)
yy=vcl.Vector(sp)
xx.data[0]=1
xx.data[1]=2
xx.data[2]=3
yy.data[0]=4
yy.data[1]=5
yy.data[2]=6
yy.data=sp.linComb(1.0,xx.data,yy.data,1.0)
print('at end')
print(yy.data)

print('Vector version')
xx.data[0]=1
xx.data[1]=2
xx.data[2]=3
yy.data[0]=4
yy.data[1]=5
yy.data[2]=6
yy.linComb(1.0,xx)
yy.myNameIs()
