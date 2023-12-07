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
    mysigma=1.e-6
    myalphalist = [1.0e-8, 5.0e-9, 2.0e-9, 1.0e-9, 5.0e-10, 2.0e-10, 1.0e-10]

    mykmax=1000
    myeps=0.00001
    myrho=0.0001
    

    data.rechdr(file='baru.su',nt=251,dt=8.0,\
                ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=-1000)

    usp=segyvc.Space('baru.su')

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
    
    print('compute predicted data')
    p = F(m0)
    
    # convolution by predicted data
    Sm0 = segyvc.ConvolutionOperator(dom=usp,
                                         rng=datasp,
                                         green=p.data)

    # construct CG solver for internal AWI use
    sol = vcalg.cgne(kmax=mykmax, eps=myeps, rho=myrho, verbose=1)

    # for awi-pen calc
    # construct product domain, range for awi sep
    dom = vcl.ProductSpace([bulksp, usp])
    rng = vcl.ProductSpace([datasp, usp, usp])

    # augmented model
    m0p = vcl.Vector(dom)
    m0p[0].copy(m0)

    # augmented data
    dp = vcl.Vector(rng)
    dp[0].copy(d)

    for i in range(len(myalphalist)):
        print('\nAWI JET COMP, alpha=' + str(myalpha))
#        op1 = awi.awiop(usp, rng, awifilt=Sm0, awipen=pen, alpha=myalpha, sigma=mysigma, observed=d)
#        [u1,ep1] = sol.solve(op1,dp)
#        x = ep1.norm()
#        j1 = 0.5*x*x
        sep = awisep(dom, rng, F, myalpha, mysigma, awisol=sol, observed=d)
        if i==0:
            sep.opfcn(m0).compj0()
        jet = vpm.vpmjet(m0, sep, dp, sol)
        j1 = jet.value()
        print('j1 = ' + str(j1))
        print('\nalpha^-2 * (j1-j0) = ' + str((j1-j0)/(myalpha*myalpha)))
        
except Exception as ex:
    print(ex)

        

