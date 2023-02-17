import os
import vcl
import segyvc
import rsfvc
import m8r
import numpy as np
import linalg

class rntobulkfb(vcl.Function):

    '''
    maps Rn dimensionless parameter vector representing velocity to
    bulk modulus for given buoyancy array.
    '''

    def __init__(self, dom, rng,\
                 buoyancy, cmin, cmax, dmin, dmax):
        ''' args:
        dom = domain (dimless velo repn)
        rng = range (bulk modulus space, same geometry)
        buoyancy = rsf filename for buoyancy
        cmin, cmax = vel bounds
        dmin, dmax = den bounds
        '''

        try:
            if not isinstance(dom, rsfvc.Space):
                raise Exception('Error: input domain not rsf space')
            if not isinstance(rng, rsfvc.Space):
                raise Exception('Error: input range not segy space')
            if not linalg.rsfcomp(dom.filename, rng.filename, checkunit=False):
                raise Exception('Error: grid geometries of dom, rng differ')
            # ignore unit of domain - should be None
            inprng = m8r.Input(rng.filename)
            if linalg.m8rstr(inprng.get('label')) != 'Bulk_modulus':
                raise Exception('Error: incorrect label for range rsf header\n \
                should be "Bulk_modulus"')
            if linalg.m8rstr(inprng.get('unit')) != 'GPa':
                raise Exception('Error: incorrect unit for range rsf header\n \
                should be "GPa"')
                
            self.dom=dom
            self.rng=rng
            
            inpbuoy = m8r.Input(buoyancy)
            if linalg.m8rstr(inpbuoy.get('label')) != 'Buoyancy':
                raise Exception('Error: incorrect label for buoyancy header\n \
                should be "Buoyancy"')
            if linalg.m8rstr(inpbuoy.get('unit')) != 'cc/g':
                raise Exception('Error: incorrect unit for buoyancy header\n \
                should be "cc/g"')            
            # cmax, cmin interpreted as m/ms or km/s
            # dmax, dmin interpreted as g/cc
            
            self.buoydata = inpbuoy.read()
            
            if np.min(self.buoydata)<=1.0/dmax or np.max(self.buoydata)>=1.0/dmin:
                raise Exception('Error: buoyancy exceeds stated bounds')
        except Exception as ex:
            print(ex)
            raise Exception('called from rntobulkfb constructor')
        else:
            self.cmax=cmax
            self.cmin=cmin
            self.dmax=dmax
            self.dmin=dmin

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def apply(self,x,y):
        inx = m8r.Input(x.data)
        xdata = inx.read()
        cdata = 0.5*(self.cmax+self.cmin) + 0.5*(self.cmax-self.cmin)*xdata\
          *(np.power(1.0 + xdata*xdata, -0.5))
        old_settings = np.seterr(divide='raise')
        ydata = cdata*cdata/self.buoydata
        iny = m8r.Input(y.data)
        datafile = linalg.m8rstr(iny.get('in'))
        ydata.tofile(datafile)
        return y

    def raw_deriv(self,x):
        pass
#        return drntobulkfb(self.dom, self.rng, self.cmax, self.cmin, x)

    def myNameIs(self):
        print('ulbounds function')


class invrntobulkfb(vcl.Function):

    '''
    maps bulk modulus to Rn dimensionless parameter vector representing
    velocity for given buoyancy array. Inverse to rntobulkfb.
    
    '''

    def __init__(self, dom, rng,\
                 buoyancy, cmin, cmax, dmin, dmax):
        try:
            if not isinstance(dom, rsfvc.Space):
                raise Exception('Error: input domain not rsf space')
            if not isinstance(rng, rsfvc.Space):
                raise Exception('Error: input range not segy space')
            if not linalg.rsfcomp(dom.filename, rng.filename, checkunit=False):
                raise Exception('Error: grid geometries of dom, rng differ')
            # ignore unit of range - should be None
            inpdom = m8r.Input(dom.filename)
            if linalg.m8rstr(inpdom.get('label')) != 'Bulk_modulus':
                raise Exception('Error: incorrect label for dom rsf header\n \
                should be "Bulk_modulus"')
            if linalg.m8rstr(inpdom.get('unit')) != 'GPa':
                raise Exception('Error: incorrect unit for dom rsf header\n \
                should be "GPa"')

            self.dom=dom
            self.rng=rng
            
            inpbuoy = m8r.Input(buoyancy)
            if linalg.m8rstr(inpbuoy.get('label')) != 'Buoyancy':
                raise Exception('Error: incorrect label for buoyancy header\n \
                should be "Buoyancy" instead of ' + m8rstr(inpbuoy.get('label')))
            if linalg.m8rstr(inpbuoy.get('unit')) != 'cc/g':
                raise Exception('Error: incorrect unit for buoyancy header\n \
                should be "cc/g"')            
            # cmax, cmin interpreted as m/ms or km/s
            # dmax, dmin interpreted as g/cc
            
            self.buoydata = inpbuoy.read()
            
            if np.min(self.buoydata)<=1.0/dmax or np.max(self.buoydata)>=1.0/dmin:
                raise Exception('Error: buoyancy exceeds stated bounds')
        except Exception as ex:
            print(ex)
            raise Exception('called from invrntobulkfb constructor')
        else:
            self.cmax=cmax
            self.cmin=cmin
            self.dmax=dmax
            self.dmin=dmin

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def apply(self,x,y):
        try:
            inx = m8r.Input(x.data)
            xdata = inx.read()
            cdata = np.power(xdata*self.buoydata,0.5)
            if np.min(cdata) <= self.cmin or np.max(cdata) >= self.cmax:
                raise Exception('Error: vel from bulk and buoy exceeds bounds')
            zdata = (cdata -    import linalg
                0.5*(self.cmax+self.cmin))/(0.5*(self.cmax-self.cmin))
            ydata = zdata*np.power(1.0-zdata*zdata,-0.5)
            iny = m8r.Input(y.data)
            datafile = linalg.m8rstr(iny.get('in'))
            ydata.tofile(datafile)
            return y
        except Exception as ex:
            print(ex)
            raise Exception('called from invrntobulkfb.apply')

    def raw_deriv(self,x):
        pass

    def myNameIs(self):
        print('invulbounds function')

# Fixed source and buoyancy op
class fsbop(vcl.Function):
    
    def __init__(self, dom, rng, \
                     buoyancy, source_p, order, sampord, nsnaps,\
                     cfl, cmin, cmax, dmin, dmax, \
                     nl1, nr1, nl2, nr2, pmlampl):
                     
        try:
            if not isinstance(dom, rsfvc.Space):
                raise Exception('Error: input domain not rsf space')
            if not isinstance(rng, segyvc.Space):
                raise Exception('Error: input range not segy space')
        except Exception as ex:
            print(ex)
            raise Exception('called from asg.fsbop constructor')
        else:
            self.dom = dom
            self.rng = rng
            # next two are pathnames to appropriate files
            self.buoyancy = buoyancy
            self.source_p = source_p
            self.order = order
            self.sampord = sampord
            self.nsnaps = nsnaps
            self.cfl = cfl
            self.cmin = cmin
            self.cmax = cmax
            self.dmin = dmin
            self.dmax = dmax
            self.nl1 = nl1
            self.nr1 = nr1
            self.nl2 = nl2
            self.nr2 = nr2
            self.pmlampl = pmlampl
            self.therest = ' order=' + str(self.order) + \
              ' sampord=' + str(self.sampord) + \
              ' nsnaps=' + str(self.nsnaps) + \
              ' cfl=' + str(self.cfl) + \
              ' cmin=' + str(self.cmin) + \
              ' cmax=' + str(self.cmax) + \
              ' dmin=' + str(self.dmin) + \
              ' dmax=' + str(self.dmax) + \
              ' nl1=' + str(self.nl1) + \
              ' nr1=' + str(self.nr1) + \
              ' nl2=' + str(self.nl2) + \
              ' nr2=' + str(self.nr2) + \
              ' pmlampl=' + str(self.pmlampl) 

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def apply(self,x,y):
        try:
            #os.system('ls /var/tmp')
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            args = ' bulkmod=' + x.data + \
                   ' buoyancy=' + self.buoyancy + \
                   ' source_p=' + self.source_p + \
                   ' data_p=' + y.data + \
                   ' deriv=0 adjoint=0' + self.therest
                   
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: return ' + str(ret) + ' from sim.x')
            #os.system('ls /var/tmp')
        except Exception as ex:
            print(ex)
            raise Exception('called from asg.apply')

    def raw_deriv(self,x):
        try:
            return fsbderiv(self.dom,self.rng,x,self.buoyancy,self.source_p,self.therest)
        except Exception as ex:
            print(ex)
            raise Exception('called from asg.fsbop.raw_deriv')

    def myNameIs(self):
        print('2D acoustic simulator: iwave/asg/main/sim.x')
        print('parameters:')
        print(self.therest)
        
class fsbderiv(vcl.LinearOperator):

    def __init__(self,dom,rng,x,buoyancy,source_p,therest):
        try:
            if not isinstance(dom, rsfvc.Space):
                raise Exception('Error: input domain not rsf space')
            if not isinstance(rng, segyvc.Space):
                raise Exception('Error: input range not segy space')
            if x.space != dom:
                raise Exception('Error: input vector not in damain')
        except Exception as ex:
            print(ex)
            raise Exception('called from asg.fsbderiv constructor')
        else:
            self.dom = dom
            self.rng = rng
            self.x = x
            self.buoyancy = buoyancy
            self.source_p = source_p
            self.therest = therest

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def applyFwd(self,dx,dy):
        try:
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            bsp = rsfvc.Space(self.buoyancy)
            dbuoyancy = vcl.Vector(bsp)
            dbuoyancy.scale(0.0)
            args = ' bulkmod=' + self.x.data + \
              ' bulkmod_d1=' + dx.data + \
              ' buoyancy=' + self.buoyancy + \
              ' buoyancy_d1=' + dbuoyancy.data + \
              ' source_p=' + self.source_p + \
              ' data_p=' + dy.data + \
              ' deriv=1 adjoint=0' + self.therest
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: return ' + str(ret) + ' from sim.x')
        except Exception as ex:
            print(ex)
            raise Exception('called from fdderiv.raw_applyFwd')

    def applyAdj(self,dx,dy):
        try:
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            bsp = rsfvc.Space(self.buoyancy)
            dbuoyancy = vcl.Vector(bsp)
            args = ' bulkmod=' + self.x.data + \
              ' bulkmod_b1=' + dy.data + \
              ' buoyancy=' + self.buoyancy + \
              ' buoyancy_b1=' + dbuoyancy.data + \
              ' source_p=' + self.source_p + \
              ' data_p=' + dx.data + \
              ' deriv=1 adjoint=1' + self.therest
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: return ' + str(ret) + ' from sim.x')
            #print('asg.fsbderiv.applyAdj ret=' + str(ret))
        except Exception as ex:
            print(ex)
            raise Exception('called from fdderiv.raw_applyAdj')

    def myNameIs(self):
        print('asg.fsbderiv: derivative of asg.fsbop simulator')

#############################
    
# Fixed buoyancy op
# example of therest:
# ' buoyancy =xxxx order=2 cfl=0.5 cmin=1.0 cmax=3.0' + \
# ' dmin=0.8 dmax=3.0 nl1=250 nr1=250 nl2=250 nr2=250 pmlampl=1.0' + \
# ' sampord=1 '

class fbop(vcl.Function):

    def __init__(self, dom, rng, \
                     buoyancy, order, sampord, nsnaps,\
                     cfl, cmin, cmax, dmin, dmax, \
                     nl1, nr1, nl2, nr2, pmlampl):
        try:
            if not isinstance(dom, vcl.ProductSpace):
                raise Exception('Error: input domain not product space')
            if len(dom.spl) != 2:
                raise Exception('Error: intput domain #factors ne 2')
            if not isinstance(dom.spl[0], rsfvc.Space):
                raise Exception('Error: input domain comp 0 not rsf space')
            if not isinstance(dom.spl[1], segyvc.Space):
                raise Exception('Error: input domain comp 1 not segy space')
        except Exception as ex:
            print(ex)
            raise Exception('called from asg constructor')
        else:
            self.dom = dom
            self.rng = rng
            ### buoyancy is pathname
            self.buoyancy = buoyancy
            self.order = order
            self.sampord = sampord
            self.nsnaps = nsnaps
            self.cfl = cfl
            self.cmin = cmin
            self.cmax = cmax
            self.dmin = dmin
            self.dmax = dmax
            self.nl1 = nl1
            self.nr1 = nr1
            self.nl2 = nl2
            self.nr2 = nr2
            self.pmlampl = pmlampl
            self.therest = ' order=' + str(self.order) + \
              ' sampord=' + str(self.sampord) + \
              ' nsnaps=' + str(self.nsnaps) + \
              ' cfl=' + str(self.cfl) + \
              ' cmin=' + str(self.cmin) + \
              ' cmax=' + str(self.cmax) + \
              ' dmin=' + str(self.dmin) + \
              ' dmax=' + str(self.dmax) + \
              ' nl1=' + str(self.nl1) + \
              ' nr1=' + str(self.nr1) + \
              ' nl2=' + str(self.nl2) + \
              ' nr2=' + str(self.nr2) + \
              ' pmlampl=' + str(self.pmlampl) 

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def apply(self,x,y):
        try:
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            args = ' bulkmod=' + x[0].data + \
                   ' buoyancy=' + self.buoyancy + \
                   ' source_p=' + x[1].data + \
                   ' data_p=' + y.data + \
                   ' deriv=0 adjoint=0' + self.therest
                   
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: failure return from sim.x')
        except Exception as ex:
            print(ex)
            raise Exception('called from asg.raw_apply')

    def raw_deriv(self,x):
        oplist = [fbpartial(self.dom,self.rng,x,0,self.buoyancy,self.therest),
            fbpartial(self.dom,self.rng,x,1,self.buoyancy,self.therest)]
        return vcl.RowLinearOperator(self.dom,self.rng,oplist)

    def myNameIs(self):
        print('2D acoustic simulator: iwave/asg/main/sim.x')
        print('parameters:')
        print(self.therest)

class fbpartial(vcl.LinearOperator):

    def __init__(self,dom,rng,x,i,buoyancy,therest):
        self.dom = dom.spl[i]
        self.rng = rng
        self.x = x
        self.i = i
        self.buoyancy = buoyancy
        self.therest = therest

    def getDomain(self):
        return self.dom

    def getRange(self):
        return self.rng

    def applyFwd(self,dx,dy):
        try:
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            if self.i==0:
                bsp = rsfvc.Space(self.buoyancy)
                dbuoyancy = vcl.Vector(bsp)
                dbuoyancy.scale(0.0)
                args = ' bulkmod=' + self.x[0].data + \
                  ' bulkmod_d1=' + dx.data + \
                  ' buoyancy=' + self.buoyancy + \
                  ' buoyancy_d1=' + dbuoyancy.data + \
                  ' source_p=' + self.x[1].data + \
                  ' data_p=' + dy.data + \
                  ' deriv=1 adjoint=0' + self.therest
            if self.i==1:
                args = ' bulkmod=' + self.x[0].data + \
                  ' buoyancy=' + self.buoyancy + \
                  ' source_p=' + dx.data + \
                  ' data_p=' + dy.data + \
                  ' deriv=0 adjoint=0' + self.therest                  
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: return ' + str(ret) + ' from sim.x')
        except Exception as ex:
            print(ex)
            raise Exception('called from fdderiv.raw_applyFwd')

    def applyAdj(self,dx,dy):
        try:
            TRIP = os.getenv('TRIP')
            cmd = os.path.join(TRIP,'iwave/asg/main/sim.x')
            if self.i==0:
                bsp = rsfvc.Space(self.buoyancy)
                dbuoyancy = vcl.Vector(bsp)
                args = ' bulkmod=' + self.x[0].data + \
                  ' bulkmod_b1=' + dy.data + \
                  ' buoyancy=' + self.buoyancy + \
                  ' buoyancy_b1=' + dbuoyancy.data + \
                  ' source_p=' + self.x[1].data + \
                  ' data_p=' + dx.data + \
                  ' deriv=1 adjoint=1' + self.therest
            if self.i==1:
                args = ' bulkmod=' + self.x[0].data + \
                  ' buoyancy=' + self.buoyancy + \
                  ' source_p=' + dy.data + \
                  ' data_p=' + dx.data + \
                  ' deriv=0 adjoint=1' + self.therest                  
            ret = os.system(cmd + args)
            if ret != 0:
                raise Exception('Error: return ' + str(ret) + ' from sim.x')
        except Exception as ex:
            print(ex)
            raise Exception('called from fdderiv.raw_applyAdj')

    def myNameIs(self):
        print('fdpartial: partial derivative ' + str(self.i))
        print('  of asg simulator')
        
