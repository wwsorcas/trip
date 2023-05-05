import vcl
import npvc

dom0 = npvc.Space(1)
dom1 = npvc.Space(1)
dom = vcl.ProductSpace([dom0,dom1])
rng  = npvc.Space(3)

x = vcl.Vector(dom)
x.data[0][0] = 1.0
x.data[1][0] = 2.0

print('x[1].data=' + str(x[1].data) + ' id=' + str(id(x[1].data)))
print('x.data[1]=' + str(x.data[1]) + ' id=' + str(id(x.data[1])))

y = vcl.Vector(dom1)
y.data[0]=3.0
print('y.data = ' + str(y.data))
x[1].copy(y)
print('x[1].data=' + str(x[1].data) + ' id=' + str(id(x[1].data)))
print('x.data[1]=' + str(x.data[1]) + ' id=' + str(id(x.data[1])))

y.data[0]=4.0
dom1.copy(y.data,x[1].data)

#x.compcopy(y,1)
print('y.data = ' + str(y.data))
print('x[1].data=' + str(x[1].data) + ' id=' + str(id(x[1].data)))
print('x.data[1]=' + str(x.data[1]) + ' id=' + str(id(x.data[1])))

z = vcl.Vector(dom)
z.copy(x)
print('z after = ' + str(z.data))


