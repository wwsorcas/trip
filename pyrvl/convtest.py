import vcl
import vcalg
import segyvc
import rsfvc
import data
import asg
import os

def runtest():
    # bulk modulus with lens
    data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)

    # bandpass filter source at sx=4200, sz=3000 (single trace)
    data.bpfilt(file='wstar.su',nt=251,dt=8.0,s=1.0,f1=1.0,f2=2.5,f3=7.5,f4=12,sx=4200,sz=3000)
    
    # create zero data file with same source position, rz=500, rx=[2000,6000]
    data.rechdr(file='g.su',nt=626,dt=8.0,ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20)
    
    bulksp = rsfvc.Space('m.rsf')
    datasp = segyvc.Space('g.su')
    
    m=vcl.Vector(bulksp,'m.rsf')
    
    ######################$ operator ############################
    F = asg.fsbop(dom=bulksp, rng=datasp, \
                buoyancy='bym.rsf', source_p='wstar.su', \
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)
    ########################## modeling ############################
    
    # evaluate F[m]
    Fm = F(m)

    data.rechdr(file='baru.su',nt=251,dt=8.0,\
            ntr=201,rx=2000,rz=1000,sx=4200,sz=3000,drx=20,delrt=-1000)
    usp=segyvc.Space('baru.su')
    baru=vcl.Vector(usp,'baru.su')

    barSm = segyvc.ConvolutionOperator(dom=usp, rng=datasp, green=Fm.data)

    # max allowed iterations
    kmax = 50
# residual norm reduction for termination
    eps  = 0.01
# normal residual norm reduction for termination
    rho  = 0.01

    # vcalg.conjgrad(x=baru, b=Fm, A=barSm, \
    #                kmax=kmax, eps=eps, rho=rho, verbose=2)

    vcalg.trconjgrad(x=baru, b=Fm, A=barSm, \
                    kmax=kmax, rho=rho, Delta=10, verbose=2)
    
runtest()    


