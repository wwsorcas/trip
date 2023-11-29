import linalg
import data
import vcl
import vcalg
import vpm
import segyvc
import rsfvc
import asg
import os
import awi

    ################## domain, range spaces ###################

try:
    wsc=1.0e+3
    mysigma=1.e-5

    mykmax=1000
    myrho=0.01
    
    data.rechdr(file='baru.su',nt=251,dt=8.0,\
                ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=-1000)

    usp=segyvc.Space('baru.su')
    u0 = vcl.Vector(usp,'baru.su')

    # bulk modulus with less focussing lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201,
                   dx=20, dz=20, lensfac=0.7)

    # homogeneous bulk modulus 
    data.model(bulkfile='m0.rsf', bulk=4.0, nx=401, nz=201,
                   dx=20, dz=20, lensfac=1.0)

    # bandpass filter source at sx=4200, sz=3000 (single trace)
    data.bpfilt(file='wstar.su',nt=251,dt=8.0,s=1.0,
                    f1=1.0,f2=2.5,f3=7.5,f4=12,sx=4200,sz=3000)
    linalg.scale('wstar.su',wsc)
    
    # create zero data file with same source position, rz=500, rx=[2000,6000]
    data.rechdr(file='g.su',nt=626,dt=8.0,ntr=201,
                    rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0)
    
    bulksp = rsfvc.Space('m.rsf')
    datasp = segyvc.Space('g.su')

    # target model
    m  = vcl.Vector(bulksp,'m.rsf')
    # homog model
    m0 = vcl.Vector(bulksp,'m0.rsf')

    F = asg.fsbop(dom=bulksp, rng=datasp, \
                buoyancy='bym.rsf', source_p='wstar.su', \
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

    print('wavelet scale = ' + str(wsc))

    print('compute observed data')
    d = F(m)

    vals = []
    errs = []
    dnorm = d.norm()
    for i in range(13):
        mm = vcl.Vector(bulksp)
        mm.copy(m0)
        mm.scale(0.1*(10-i))
        mm.linComb(0.1*i,m)
        wg = awi.awiwg(dom=usp, sim=F, mod=mm, data=d, sigma=mysigma, kmax=mykmax, rho=myrho, precond=1, verbose=2)
        vals.append(wg.value())
        errs.append(wg.innererr().norm()/dnorm)
    print('\nvals, rel errors at mtest = (1-t)*m0 + t*m')
    print('  t       val       relerr')
    for i in range(13):
        print('%2.1f  %10.4e  %10.4e' % (0.1*i, vals[i], errs[i]))
                                            
# construct jet at m0

    #wg = awi.awiwg(dom=usp, sim=F, mod=m0, data=d, sigma=mysigma, kmax=mykmax, rho=myrho, precond=1, verbose=2)

    # homogeneous bulk modulus 
    #data.model(bulkfile='g.rsf', bulk=4.0, nx=401, nz=201,
#                   dx=20, dz=20, lensfac=1.0)

    #g = vcl.Vector(bulksp,'g.rsf')
    #g.copy(wg.gradient())

    #u0.copy(wg.innersol())
    
except Exception as ex:
    print(ex)

        

