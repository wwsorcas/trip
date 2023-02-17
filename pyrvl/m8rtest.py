import m8r
import numpy as np
import data
import vcl
import asg
import rsfvc
import os

try:
    # bulk modulus with less focussing lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)
    inp = m8r.Input('m.rsf')
    mdata=inp.read()
    
    outp = m8r.Output('cx.rsf')
    outp.put('in','cx.rsf@')
    outp.put('label','')
    outp.put('unit','')
    
    outp.write(mdata)
    
    bulksp = rsfvc.Space('m.rsf')
    cxsp   = rsfvc.Space('cx.rsf')
    #cxsp.myNameIs()
    
    m = vcl.Vector(bulksp,'m.rsf')
    
    cmax = 6.0
    cmin = 1.0
    dmin = 0.5
    dmax = 3.0
    
    inv = asg.invrntobulkfb(bulksp, cxsp, 'bym.rsf', cmin, cmax, dmin, dmax)
    
    cxv = inv(m)
    
    cx = vcl.Vector(cxsp,'cx.rsf')
    
    cx.copy(cxv)

    fwd = asg.rntobulkfb(cxsp, bulksp, 'bym.rsf', cmin, cmax, dmin, dmax)

    mx = fwd(cxv)

    os.system('sfcp ' + mx.data + ' mx.rsf')

    dmx = vcl.Vector(bulksp)

    dmx.copy(m)
    dmx.linComb(-1.0,mx)

    print('rel error = ' + str(dmx.norm()/m.norm()))

except Exception as ex:
    print(ex)



