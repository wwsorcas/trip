import m8r
import numpy as np
import data
import vcl
import asg
import rsfvc
import os
import linalg

try:
    print('test 1: transform bulk mod to rn rep, transform back, compute error')
    
    # bulk modulus with less focussing lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)
    inp = m8r.Input('m.rsf')
    mdata=inp.read()

    # create expanded velo space
    os.system('sfput < m.rsf label="Real" unit="None" > cx.rsf')
    
    bulksp = rsfvc.Space('m.rsf')
    cxsp   = rsfvc.Space('cx.rsf')

    # bulk mod vector
    m = vcl.Vector(bulksp,'m.rsf')

    # set initial bounds
    cmax = 6.0
    cmin = 1.0
    dmin = 0.5
    dmax = 3.0

    # inverse expanded c
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

    print('\ntest 2: check linearization error convergence')
    print('  random pert uniform in [-0.1,0.1], step *= 0.5')
    print('  error should decrease by approx 0.5 per step')
              

    # perturbation, perturbed
    dcx = vcl.Vector(cxsp)
    pcx = vcl.Vector(cxsp)

    # random pert of every sample, uniform in [-0.1,0.1]
    linalg.rand(dcx.data)
    dcx.scale(0.1)

    # derivative at m0
    der = asg.drntobulkfb(cxsp, bulksp, cxv, 'bym.rsf', cmin, cmax, dmin, dmax)

    step = 1.0
    print('step         norm of lin error')
    for i in range(0,5):
        dm = der*dcx
        pcx.copy(cxv)
        pcx.linComb(1.0,dcx)
        pm = fwd(pcx)
        pm.linComb(-1.0,m)
        pm.linComb(-1.0,dm)
        print('%10.6e  %10.6e' % (step,pm.norm()))
        step *= 0.5
        dcx.scale(0.5)          
        
    print('\ntest 3: inverse mapping withwrong bounds, cmax too small')
    print('  should throw exception')
    cmax = 1.5
    winv = asg.invrntobulkfb(bulksp, cxsp, 'bym.rsf', cmin, cmax, dmin, dmax)
    wxv = winv(m)

except Exception as ex:
    print(ex)



