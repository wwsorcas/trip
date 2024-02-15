import os
import sys

import linalg
import data
import vcl
import vcalg
import vpm
import segyvc
import rsfvc
import asg
import awi

    ################## domain, range spaces ###################

try:

    # arg1 = job choice: 'vals', 'grad0', 'grad', or gradtest (default)
    # arg2 = frequency index
    # arg3 = lens index
    # arg4 = alpha
    # arg5 = sigma
    # arg6 = rho
    
    choice = sys.argv[1]
    idx = int(sys.argv[2])
    lens = int(sys.argv[3])
    myalpha = float(sys.argv[4])
    mysigma = float(sys.argv[5])
    myrho = float(sys.argv[6])

    print('mswitest:')
    print('choice = ' + choice)
    print('freq index = ' + str(idx))
    print('lens index = ' + str(lens))
    print('alpha = ' + str(myalpha))
    print('sigma = ' + str(mysigma))
    print('rho   = ' + str(myrho))
    print(' ')
    
    wsc=1.0e+3
    mykmax=1000

    NTU = [251,501]
    NTD = [626,1251]
    NX  = [401, 801]
    NZ  = [201, 401]
    DELRTW = [-1000, -500]
    DT  = [8, 4]
    DX  = [20, 10]
    F1  = [1.0, 2.0]
    F2  = [2.5, 5.0]
    F3  = [7.5, 15.0]
    F4  = [12.5, 25.0]

    LENSRADD = [0.2, 0.4]
    LENSFAC = [0.7, 0.5]

    print('build spaces, models')
    
    data.rechdr(file='baru.su',nt=NTU[idx],dt=DT[idx],\
                ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=-1000)

    usp=segyvc.Space('baru.su')
    u0 = vcl.Vector(usp,'baru.su')

    # bulk modulus with less focussing lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=NX[idx], nz=NZ[idx],
                   dx=DX[idx], dz=DX[idx], lensfac=LENSFAC[lens], lensradd=LENSRADD[lens])

    # homogeneous bulk modulus 
    data.model(bulkfile='m0.rsf', bulk=4.0, nx=NX[idx], nz=NZ[idx],
                   dx=DX[idx], dz=DX[idx], lensfac=1.0)

    # bandpass filter source at sx=4200, sz=3000 (single trace)
    data.bpfilt(file='wstar.su',nt=NTU[idx],dt=DT[idx],s=1.0,
                    f1=F1[idx],f2=F2[idx],f3=F3[idx],f4=F4[idx],sx=4200,sz=3000)
    linalg.scale('wstar.su',wsc)
    
    # create zero data file with same source position, rz=500, rx=[2000,6000]
    data.rechdr(file='g.su',nt=NTD[idx],dt=DT[idx],ntr=201,
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

    print('compute observed data at m, predicted data at m0')
    os.system('/bin/cp g.su d.su')
    os.system('/bin/cp g.su d0.su')
    d = vcl.Vector(datasp,'d.su')
    d.copy(F(m))
    
    if choice == 'vals':
        print('compute values')        
        vals = []
        errs = []
        dnorm = d.norm()

        for i in range(13):
            print('t = ' + str(0.1*(10-i)))
            mm = vcl.Vector(bulksp)
            mm.copy(m0)
            mm.scale(0.1*(10-i))
            mm.linComb(0.1*i,m)
            j = awi.mswi(dom=usp, sim=F, mod=mm, data=d, alpha=myalpha, sigma=mysigma, kmax=mykmax, rho=myrho, verbose=2)
            vals.append(j.value())
            errs.append(j.dataerr().norm()/dnorm)
        print('\nvals, rel errors at mtest = (1-t)*m0 + t*m')
        print('alpha = %10.4e sigma = %10.4e rho = %10.4e'
                      % (myalpha, mysigma, myrho))
        print('  t       val       relerr')
        for i in range(13):
            print('%2.1f  %10.4e  %10.4e' % (0.1*i, vals[i], errs[i]))
                                            
# construct jet at m0

    else:

        if choice == 'grad0':
        # homogeneous bulk modulus 
            j0 = awi.mswi(dom=usp, sim=F, mod=m0, data=d, alpha=myalpha, sigma=mysigma, kmax=mykmax, rho=myrho, verbose=2)
            data.model(bulkfile='grad0.rsf', bulk=4.0, nx=NX[idx], nz=NZ[idx],
                    dx=DX[idx], dz=DX[idx], lensfac=1.0)
            g = vcl.Vector(bulksp,'grad0.rsf')
            g.copy(j0.gradient())

        elif choice == 'grad':
            # target bulk modulus
            j = awi.mswi(dom=usp, sim=F, mod=m, data=d, alpha=myalpha, sigma=mysigma, kmax=mykmax, rho=myrho, verbose=2)
            data.model(bulkfile='grad.rsf', bulk=4.0, nx=NX[idx], nz=NZ[idx],
                    dx=DX[idx], dz=DX[idx], lensfac=1.0)
            g = vcl.Vector(bulksp,'grad.rsf')
            g.copy(j.gradient())

        elif choice == 'gradtest':

            data.model(bulkfile='dm.rsf', bulk=4.0, nx=NX[idx], nz=NZ[idx],
                    dx=DX[idx], dz=DX[idx], lensfac=1.0)
            dm = vcl.Vector(bulksp,'dm.rsf')
            dm.copy(m)
            dm.linComb(-1.0,m0)
            j0 = awi.mswi(dom=usp, sim=F, mod=m0, data=d, alpha=myalpha, sigma=mysigma, kmax=mykmax, rho=myrho, verbose=2)

            mtest = vcl.Vector(bulksp)
            g = vcl.Vector(bulksp)
            g.copy(j0.gradient())
            v0 = j0.value()

            h = 0.1
            steps = 5
            secs = []
            for ih in range(steps):
                mtest.copy(m0)
                mtest.linComb(h,dm)
                jtest = awi.mswi(dom=usp, sim=F, mod=mtest, data=d,
                            alpha=myalpha, sigma=mysigma,
                            kmax=mykmax, rho=myrho, verbose=0)
                secs.append((jtest.value()-v0)/h)
                h *= 0.5

            print('\n\nmswi gradient test')
            print('alpha = %10.4e sigma = %10.4e rho = %10.4e'
                      % (myalpha, mysigma, myrho))
            print('  v at m0  = %10.4e' % (v0))
            print('  slope    = %10.4e' % (g.dot(dm)))
            h = 0.1
            print('  h           sec         err/h')

            for ih in range(steps):
                print('%6.4f     %10.4e  %10.4e'
                          % (h, secs[ih], abs(secs[ih]-g.dot(dm))/h))                
                h *= 0.5
        else:
            print('no choice at all!')
    
except Exception as ex:
    print(ex)

        

