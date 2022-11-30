import os
import tempfile

def bpfilt(file,nt,dt,s,f1,f2,f3,f4,sx,sz):

    CWPROOT = os.getenv('CWPROOT')
    suspike = os.path.join(CWPROOT,'bin/suspike')
    sufilter = os.path.join(CWPROOT,'bin/sufilter')
    sugain = os.path.join(CWPROOT,'bin/sugain')
    sushw = os.path.join(CWPROOT,'bin/sushw')

    cmd = \
    suspike + ' nt=' + str(nt) + ' ntr=1 offset=0 ix1=1 nspk=1 it1=' + \
    str((1000*s/dt)+1) + ' dt=' + str(0.001*dt) + ' | ' + \
    sufilter + ' f=' + str(f1) + ',' + str(f2) + ',' + \
    str(f3) + ',' + str(f4) + ' | ' + \
    sugain + ' scale=' + str(1.0/dt) + ' | ' + \
    sushw + ' key=delrt,gelev,selev,gx,sx a=0.0,' + str(-sz) + ',' + \
    str(-sz) + ',' + str(sx) + ',' + str(sx) + \
    ' > ' + file

#    print(cmd)

    os.system(cmd)

# zero-phase delta
def delta(file,nt,dt,sx,sz):

    CWPROOT = os.getenv('CWPROOT')
    suspike = os.path.join(CWPROOT,'bin/suspike')
    sugain = os.path.join(CWPROOT,'bin/sugain')
    sushw = os.path.join(CWPROOT,'bin/sushw')

    cmd = \
    suspike + ' nt=' + str(nt) + ' ntr=1 offset=0 ix1=1 nspk=1 it1=' +\
    str((nt/2)+1) + ' dt=' + str(0.001*dt) + ' | ' + \
    sugain + ' scale=' + str(1.0/dt) + ' | ' + \
    sushw + ' key=delrt,gelev,selev,gx,sx a=' + str(-dt*(nt-1)/2) + ',' + \
    str(-sz) + ',' + str(-sz) + ',' + str(sx) + ',' + str(sx) + \
    ' > ' + file

#    print(cmd)

    os.system(cmd)

# gather of zero-phase deltas - iteration over delta
def deltagather(file,nt,dt,ntr,sxstart,szstart,dsx,dsz):

    delta(file,nt,dt,sxstart,szstart)

    temp = tempfile.NamedTemporaryFile(delete=False)

    for i in range(1,ntr):
        delta(temp.name,nt,dt,sxstart+i*dsx,szstart+i*dsz)
        os.system('/bin/cat ' + temp.name + ' >> ' + file)
                      
    temp.close()
    os.unlink(temp.name)
    
# trace header
# compute header files for receiver (rechdr) and source (srchdr) lines, at gelev
# = -1000 and -3000 resp.

def rechdr(file,nt,dt,rxmin,rxmax,ntr,rz,sx,sz):

    CWPROOT = os.getenv('CWPROOT')
    sunull = os.path.join(CWPROOT,'bin/sunull')
    sugain = os.path.join(CWPROOT,'bin/suscale')
    sushw = os.path.join(CWPROOT,'bin/sushw')
    suchw = os.path.join(CWPROOT,'bin/suchw')

    cmd = sunull + ' nt=' + str(nt) + \
         ' ntr=' + str(ntr) + \
         ' dt=' + str(0.001*dt) + ' | ' + \
         sushw + ' key=gx a=' + str(rxmin) + \
         ' b=' + str((rxmax-rxmin)/(ntr-1)) + \
         ' j=' + str(ntr) + ' | ' + \
         sushw + ' key=gelev,selev,sx' + \
         ' a=' + str(-rz) + ',' + str(-sz) + ',' + str(sx) + ' | ' + \
         suchw + ' key1=offset key2=gx key3=sx c=-1 > ' + file         

#    print(cmd)

    os.system(cmd) 

# spatial model - only homog buoy, optional lens in bulk
# value at ctr of lens = lensfac*bulk
def model(bulkfile, bulk, nx, nz, dx, dz, lensfac, buoy=1.0):
    RSFROOT = os.getenv('RSFROOT')
    makevel = os.path.join(RSFROOT,'bin/sfmakevel')
    put   = os.path.join(RSFROOT,'bin/sfput')
    cmd = makevel + \
         ' n1=' + str(nz) + \
         ' n2=' + str(nx) + \
         ' d1=' + str(dz) + \
         ' d2=' + str(dx) + \
         ' label1=Depth label2=Distance' + \
         ' unit1=m unit2=m' + \
         ' label=Bulk_modulus unit=GPa' +\
         ' v000=' + str(bulk) + \
         ' x1lens=' + str(nz*dz/2.0) + \
         ' x2lens=' + str(nx*dx/2.0) + \
         ' dlens=' + str(nz*dz/5.0) + \
         ' tlens=' + str(nx*dx/5.0) + \
         ' vlens=' + str(bulk*(lensfac-1.0)) + ' | ' +\
         ' sfput dim=2 gdim=2 id1=0 id2=1 > ' + bulkfile + ' && ' + \
         makevel + \
         ' n1=' + str(nz) + \
         ' n2=' + str(nx) + \
         ' d1=' + str(dz) + \
         ' d2=' + str(dx) + \
         ' label1=Depth label2=Distance' + \
         ' unit1=m unit2=m' + \
         ' label=Buoyancy unit=cc/g' +\
         ' v000=' + str(buoy) + ' | ' +\
         ' sfput dim=2 gdim=2 id1=0 id2=1 > by' + bulkfile 

#    print(cmd)

    os.system(cmd)

