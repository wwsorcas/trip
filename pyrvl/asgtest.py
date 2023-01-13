import linalg
import data
import vcl
import segyvc
import rsfvc
import asg
import os

################### input data #####################

# bulk modulus with less focussing lens
data.model(bulkfile='m.rsf', bulk=4.0, nx=401, nz=201, dx=20, dz=20, lensfac=0.7)

# bandpass filter source at sx=4200, sz=3000 (single trace)
data.bpfilt(file='wstar.su',nt=251,dt=8.0,s=1.0,f1=1.0,f2=2.5,f3=7.5,f4=12,sx=4200,sz=3000)

# create zero data file with same source position, rz=500, rx=[2000,6000]
data.rechdr(file='g.su',nt=626,dt=8.0,rxmin=2000,rxmax=6000,ntr=201,rz=1000,sx=4200,sz=3000)

rsfroot = os.getenv('RSFROOT')
sfcp = os.path.join(rsfroot,'bin/sfcp')
sfrm = os.path.join(rsfroot,'bin/sfrm')
sfnoise = os.path.join(rsfroot,'bin/sfnoise')
sfscale = os.path.join(rsfroot,'bin/sfscale')
sfboxsmooth = os.path.join(rsfroot,'bin/sfboxsmooth')
cmd = sfcp + ' m.rsf dm.rsf; ' + sfnoise + ' < dm.rsf > ndm.rsf; ' + sfboxsmooth + ' < ndm.rsf > dm.rsf rect1=10 rect2=10 rept=5'
os.system(cmd)

def runtests():
    
    ################## domain, range spaces ###################

    bulksp = rsfvc.Space('m.rsf')
    srcsp = segyvc.Space('wstar.su')
    domlist = [bulksp,srcsp]
    dom = vcl.ProductSpace(domlist)
    rng = segyvc.Space('g.su')

    ################## construct FB operator ###################

    F = asg.fbop(dom,rng,buoyancy='bym.rsf', order=2, sampord=1, nsnaps=20, cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0, nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

    print('\n*** operator')
    F.myNameIs()

    ################### initialize input vector ##################

    print('\n*** domain space')
    dom.myNameIs()
    
    x = vcl.Vector(dom,['m.rsf', 'wstar.su'])
    
    print('\n*** WHOLE VECTOR')
    x.myNameIs()
    
    print('\n*** COMPONENT 0')
    x[0].myNameIs()
    
    print('\n*** COMPONENT 1')
    x[1].myNameIs()
    
    xx = vcl.Vector(dom)
    
    xx.copy(x)
    
    ##################### apply operator ########################
    
    #y = vcl.Vector(rng)
    #F.apply(x,y)
    y = F(xx)
    
    print('\n*** OUTPUT VECTOR ***')
    y.myNameIs()

    ##################### dot product test #################

    print('\n*** DOT PRODUCT TEST - FIXED BUOYANCY')

    dm = vcl.Vector(dom.spl[0],'dm.rsf')

    DFx = F.deriv(xx)

    dy = DFx[0]*dm
    dydotdy =  dy.dot(dy)
    dmdotdfstardy = dm.dot(vcl.transp(DFx[0])*dy)
    relerr = abs(dydotdy - dmdotdfstardy)/dydotdy

    print('\n***********************************************')
    print('DFx*dm dot DFx*dm         = %12.6e' % (dydotdy))
    print('dm dot transp(DFx)*DFx*dm = %12.6e' % (dmdotdfstardy))
    print('relative error            = %12.6e' % (relerr))

    print('\n*** DOT PRODUCT TEST - FIXED SOURCE AND BUOYANCY')

    G = asg.fsbop(dom[0],rng,buoyancy='bym.rsf', source_p='wstar.su', order=2, sampord=1, nsnaps=20, cfl=0.5, cmin=1.0, cmax=3.0,dmin=0.8, dmax=3.0, nl1=250, nr1=250, nl2=250, nr2=250, pmlampl=1.0)

    print('\n*** operator')
    G.myNameIs()    

    DGx = G.deriv(xx[0])

    dy = DGx*dm

    dydotdy =  dy.dot(dy)
    dmdotdfstardy = dm.dot(vcl.transp(DGx)*dy)
    relerr = abs(dydotdy - dmdotdfstardy)/dydotdy

    print('\n***********************************************')
    print('DFx*dm dot DFx*dm         = %12.6e' % (dydotdy))
    print('dm dot transp(DFx)*DFx*dm = %12.6e' % (dmdotdfstardy))
    print('relative error            = %12.6e' % (relerr))    
    
runtests()

####################### clean up globals ########################

cmd1 = sfrm + ' m.rsf dm.rsf ndm.rsf bym.rsf'
os.system(cmd1)
os.unlink('g.su')
os.unlink('wstar.su')
os.unlink('cout0.txt')
os.unlink('jnk')














