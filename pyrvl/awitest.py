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
    linalg.scale('wstar.su',1.0e+3)
    
    # create zero data file with same source position, rz=500, rx=[2000,6000]
    data.rechdr(file='g.su',nt=626,dt=8.0,ntr=201,
                    rx=2000.0,rz=1000.0,sx=4200,sz=3000,drx=20.0)
    
    bulksp = rsfvc.Space('m.rsf')
    datasp = segyvc.Space('g.su')
    
    F = asg.fsbop(dom=bulksp, rng=datasp, \
                buoyancy='bym.rsf', source_p='wstar.su', \
                order=2, sampord=1, nsnaps=20,\
                cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0,\
                nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

    print('construct product domain and range')
    dom = vcl.ProductSpace([bulksp,usp])
    rng = vcl.ProductSpace([datasp,usp])
                
    print('compute lens data')
    m = vcl.Vector(bulksp,'m.rsf')
    d = vcl.Vector(rng)
    d[0].copy(F(m))
    barSm = segyvc.ConvolutionOperator(dom=usp, rng=datasp, 
                                        green=d[0].data)
    print('compute hom data')
    m0 = vcl.Vector(bulksp,'m0.rsf')
    d0 = vcl.Vector(rng)
    d0[0].copy(F(m0))
    barSm0 = segyvc.ConvolutionOperator(dom=usp, rng=datasp, 
                                        green=d0[0].data)

    #print('construct AWI penalty operator, precond=0')
    #op0 = awi.awipensol(dom=usp, alpha=1.e-4)
    #print('construct AWI operator, precond=0')
    #op0 = awiop.awiop(Fm, usp, alpha=1.0e-4)


    print('construct augmented rhs')
    

    print('construct CG solver for internal AWI use')
    sol = vcalg.cgne(kmax=20, eps=0.01, rho=0.01, verbose=1)
    
    print('construct awisep, precond=0')
    sep0 = awi.awisep(dom, rng, F, alpha=1.0e-4)

    print('construct solver for use in jet (use same)')
    print('construct jet at m')
    j01 = vpm.vpmjet(m, sep0, d, sol)
    print('value at m  = ' + str(j01.value()))
    j00 = vpm.vpmjet(m0, sep0, d, sol)
    print('value at m0 = ' + str(j00.value()))
    print('extract awiop, precond=0')
    op0 = sep0.opfcn(m)

    barull = vcl.Vector(usp)
    baruhl = vcl.Vector(usp)
    pent0ll = vcl.Vector(usp)
    pent1ll = vcl.Vector(usp)
    pent0hl = vcl.Vector(usp)
    pent1hl = vcl.Vector(usp)

    #print('construct AWI penalty operator, precond=1, case lens-lens')
    #opll1 = awi.awipensol(dom=usp, alpha=1.0e-4, solver=sol, d=Fm, p=Fm)

    print('construct AWI SEPFunction, precond=1, case lens-lens')
    #opll1 = awi.awiop(dom[1], rng, predicted=Fm, alpha=1.0e-4, awisol=sol, observed=Fm)
    sepll1 = awi.awisep(dom, rng, F, alpha=1.0e-4, awisol=sol, observed=d[0])
    #print('construct AWI penalty operator, precond=1, case hom-lens')
    opll1 = sepll1.opfcn(m)
    
    print('construct AWI SEPFunction, precond=1, case hom-lens')
    #opll1 = awi.awiop(dom[1], rng, predicted=Fm0, alpha=1.0e-4, awisol=sol, observed=Fm)
    sephl1 = awi.awisep(dom, rng, F, alpha=1.0e-4, awisol=sol, observed=d[0])
    #print('construct AWI penalty operator, precond=1, case hom-lens')
    ophl1 = sephl1.opfcn(m0)
    
    dir = {
        'lens-lens': [d[0], barSm, barull, pent0ll, pent1ll, opll1],
        'lens-hom': [d0[0], barSm0, baruhl, pent0hl, pent1hl, ophl1]
        }

    for po in ['lens-lens', 'hom-lens']:
        print('compute awi kernel in case ' + po)
        dir[po][2]=vcl.Vector(usp)
        vcalg.conjgrad(x=dir[po][2], b=dir[po][0], A=dir[po][1], \
                    kmax=20, eps=0.01, rho=0.01, verbose=0)

        linalg.simplot(dir[po][2].data)

        print('apply AWI penalty operator, precond=1, case ' + po)
        dir[po][4] = dir[po][5]*dir[po][2]
        linalg.simplot(dir[po][4][1].data)

        print('apply AWI penalty operator, precond=0, case ' + po)
        dir[po][3] = op0*dir[po][2]
        linalg.simplot(dir[po][3][1].data)

    
    

    

except Exception as ex:
    print(ex)

        


#rhs[0].copy(Fm)
#e = vcl.Vector(op.getRange())

#vcalg.conjgrad(x=baru, b=rhs, A=op, kmax=50, eps=0.01, rho=0.01, verbose=2)



#x1 = vcl.transp(op)*rhs
#y1 = op*x1
#x2 = vcl.transp(op)*y1
#y2 = op*x2
#x3 = vcl.transp(op)*y2

#print('x1 dot op^Ty2 = ' + str(x1.dot(x3)))
#print('op x1 dot y2 = ' + str(y1.dot(y2)))
