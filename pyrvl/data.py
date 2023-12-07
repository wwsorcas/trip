import os
import tempfile

def hello():
    print('hello from data')

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

# gather of zero phase deltas - iteration over delta
def deltagather(file,nt,dt,ntr,sxstart,szstart,dsx,dsz):

    delta(file,nt,dt,sxstart,szstart)

    temp = tempfile.NamedTemporaryFile(delete=False)

    for i in range(1,ntr):
        delta(temp.name,nt,dt,sxstart+i*dsx,szstart+i*dsz)
        os.system('/bin/cat ' + temp.name + ' >> ' + file)
                      
    temp.close()
    os.unlink(temp.name)

# gather of bandpass filters - iteration over bpfilt
def bpfiltgather(file,nt,dt,s,f1,f2,f3,f4,ntr,sxstart,szstart,dsx,dsz):

    bpfilt(file,nt,dt,s,f1,f2,f3,f4,sxstart,szstart)

    temp = tempfile.NamedTemporaryFile(delete=False)

    for i in range(1,ntr):
        bpfilt(temp.name,nt,dt,s,f1,f2,f3,f4,sxstart+i*dsx,szstart+i*dsz)
        os.system('/bin/cat ' + temp.name + ' >> ' + file)
                      
    temp.close()
    os.unlink(temp.name)
    
def rechdr(file,nt,dt,ntr,rx,rz,sx,sz,drx,
               delrt=0,nshot=1,dsx=0,fixed=True):
    '''
    Creates a SU (SEGY without reel header) file with zero 
    data values and correct header (metadata) values describing 
    the geometry of a 2D seismic survey. Produces headers for either 
    fixed spread (same receiver positions for every shot) or towed
    streamer (same receiver offsets for every shot) geometry. Geometry
    is regular, that is, same step between each pair of adjacent 
    reeivers and each pair of adjacent shot locations. All shot 
    records are assumed to have the same number of receivers. Also, all 
    receivers, respectively all sources, are assumed to have the same
    depth (z), and differ only in horizontal position (x). These 
    properties are idealizations, suitable for synthetic data, and
    unlikely to hold for field data.

    Parameter units are m for distances/lengths and ms for times.

    Parameters:
        file (str): filename for output, typically with suffix '.su'
        nt (int): number if time samples in each receiver trace
        dt (float): time sample interval
        ntr (int): number of receivers per shot
        rx (float): first receiver horizontal position (fixed spread)
            or offset (towed streamer)
        rz (float): receiver depth
        sx (float): first shot horizontal position
        sz (float): shot depth
        drx (float): interval between receiver positions within
            shot record (all belonging to same shot)
        delrt (float): time of first sample
        nshot (int): number of shots in survey
        dsx (float): interval between shot positions
        fixed (boolean): fixed spread if True, else towed streamer
        return value (int): 0 if successful
    '''
        
    CWPROOT = os.getenv('CWPROOT')
    sunull = os.path.join(CWPROOT,'bin/sunull')
    sugain = os.path.join(CWPROOT,'bin/suscale')
    sushw = os.path.join(CWPROOT,'bin/sushw')
    suchw = os.path.join(CWPROOT,'bin/suchw')

    if fixed:
        cmd = sunull + ' nt=' + str(nt) + \
          ' ntr=' + str(ntr*nshot) + \
          ' dt=' + str(0.001*dt) + ' | ' + \
          sushw + ' key=gx a=' + str(rx) + \
          ' b=' + str(drx) + \
          ' j=' + str(ntr) + ' | ' + \
          sushw + ' key=sx a=' + str(sx) + \
          ' c=' + str(dsx) + \
          ' j=' + str(ntr) + ' | ' + \
          sushw + ' key=gelev,selev,delrt' + \
          ' a=' + str(-rz) + ',' + str(-sz) + ',' + str(delrt) + ' | ' + \
          suchw + ' key1=offset key2=gx key3=sx c=-1 > ' + file         
    else:
        cmd = sunull + ' nt=' + str(nt) + \
          ' ntr=' + str(ntr*nshot) + \
          ' dt=' + str(0.001*dt) + ' | ' + \
          sushw + ' key=gx a=' + str(rx + sx) + \
          ' b=' + str(drx) + \
          ' c=' + str(dsx) + \
          ' j=' + str(ntr) + ' | ' + \
          sushw + ' key=sx a=' + str(sx) + \
          ' c=' + str(dsx) + \
          ' j=' + str(ntr) + ' | ' + \
          sushw + ' key=gelev,selev,delrt' + \
          ' a=' + str(-rz) + ',' + str(-sz) + ',' + str(delrt) + ' | ' + \
          suchw + ' key1=offset key2=gx key3=sx c=-1 > ' + file              
          
    ret = os.system(cmd)

    if ret != 0:
       print('\nattempt to create shot record file ' + file + \
                 ' via function rechdr failed')

    return ret

def rsffile(file, datatype, unit, n1, n2, d1, d2, val=1.0):
    '''
    creates 2D RSF data file pair, spatially homogeneous data,
    suitable for further processing via NumPy and m8r.

    For all parameters, unit of length is m.

    Parameters:
        file (str): name of header (rsf) file - data file will be DATAPATH/file@
        datatype (str): data type, eg. velocity, density,... - no embedded blanks!
        unit (str): data unit - no embedded blanks
        n1 (int): number of points on axis 1 (depth)
        n2 (int): number of points on axis 2 (horizontal)
        d1 (float): increment on axis 1 (depth)
        d2 (float): increment on axis 2 (horizontal)
        val (float): value assigned to all data points
        return value (int): return from os.system, = 0 for success

    Also initializes header words dim, gdim, id1, id2, needed by
    some IWAVE applications but ignored by M8R
    '''

    RSFROOT = os.getenv('RSFROOT')
    makevel = os.path.join(RSFROOT,'bin/sfmakevel')
    put   = os.path.join(RSFROOT,'bin/sfput')
    cmd = makevel + \
         ' n1=' + str(n1) + \
         ' n2=' + str(n2) + \
         ' d1=' + str(d1) + \
         ' d2=' + str(d2) + \
         ' v000=' + str(val) + ' | ' +\
         put + ' dim=2 gdim=2 id1=0 id2=1 | ' +\
         put + ' unit1=m unit2=m | ' +\
         put + ' label1=Depth label2=Distance | ' +\
         put + ' label=' + datatype + ' unit=' + unit + ' > ' + file

    ret = os.system(cmd)

    if ret != 0:
       print('\nattempt to create file ' + file + \
                 ' via function rsffile failed')

    return ret

def model(bulkfile, bulk, nx, nz, dx, dz, lensfac, buoy=1.0, lensradd=0.2, lensradt=0.2):
    ''' 
    creates (bulk modulus, buoyancy) pair of rsf pairs for input
    to iwave simulation code. bulk modulus has optional circular 
    lens in center, with gaussian profile.

    For all length parameters, unit = m.

    Parameters:

        bulkfile (str): name of bulk modulus header file. data file is 
            DATAPATH/bulkfile@, buoyancy header is "by" + bulkfile
        bulk (float): background bulk modulus, unit = GPa
        nx (int): number of points on axis 2 (horizontal)
        nz (int): number of points on axis 1 (depth)
        dx (float): increment on axis 2 (horizontal)
        dz (float): increment on axis 1 (depth)
        lensfac (float): relative bulkmod change from background
            to center of lens
        buoy (float): buoyancy (homogeneous), unit = cc/g
        lensradd (float): lens diameter, as proportion of total range
        lensradt (float): lens thickness, as proportion of total range
        return value (int): sum of returns from os.system, = 0 for success
    '''
    
    RSFROOT = os.getenv('RSFROOT')
    makevel = os.path.join(RSFROOT,'bin/sfmakevel')
    put   = os.path.join(RSFROOT,'bin/sfput')
    cmd = makevel + \
         ' n1=' + str(nz) + \
         ' n2=' + str(nx) + \
         ' d1=' + str(dz) + \
         ' d2=' + str(dx) + \
         ' v000=' + str(bulk) + \
         ' x1lens=' + str(nz*dz/2.0) + \
         ' x2lens=' + str(nx*dx/2.0) + \
         ' dlens=' + str(nz*dz*lensradd) + \
         ' tlens=' + str(nx*dx*lensradt) + \
         ' vlens=' + str(bulk*(lensfac-1.0)) + ' | ' +\
         put + ' dim=2 gdim=2 id1=0 id2=1 ' + ' | ' +\
         put + ' unit1=m unit2=m' + ' | ' +\
         put + ' label1=Depth label2=Distance ' + ' | ' +\
         put + ' label=Bulk_modulus unit=GPa > ' + bulkfile

    ret = os.system(cmd)

    if ret != 0:
       print('\nattempt to crete bulk modulus file ' + bulkfile + \
                 ' via function model failed')

    retby = rsffile('by' + bulkfile, 'Buoyancy', 'cc/g',
                            nz, nx, dz, dx, val=buoy)

    if retby != 0:
       print('\nattempt to crete buoyancy file by' + bulkfile + \
                 ' via function model failed')

    return ret + retby


    
