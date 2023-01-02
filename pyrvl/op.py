import os
import linalg

def fdop(m,w,d,adj=0):

    # sanity-check argsd
    if not linalg.sanity(m,'rsf'):
        return False
    if not linalg.sanity(w,'su'):
        return False
    if not linalg.sanity(d,'su'):
        return False
    
    # check that top-level TRIP directory exists
    TRIP = os.getenv('TRIP')
    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # build command including many variables not updated
    sim = os.path.join(TRIP,'iwave/asg/main/sim.x')
    cmd = sim + \
        ' bulkmod=' + m + \
        ' buoyancy=by' + m + \
        ' data_p=' + d + \
        ' source_p=' + w + \
        ' deriv=0 adjoint=' + str(adj) + ' order=2 cfl=0.5 cmin=1.0 cmax=3.0' + \
        ' dmin=0.8 dmax=3.0 nl1=250 nr1=250 nl2=250 nr2=250 pmlampl=1.0' + \
        ' sampord=1 '
    #
    # execute
#    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print('command failed with return value ' + str(ret))
        return False

    return True

def convop(g,w,d,adj=0):

    # sanity-check argsd
    if not linalg.sanity(g,'su'):
        return False
    if not linalg.sanity(w,'su'):
        return False
    if not linalg.sanity(d,'su'):
        return False
    
    # check that top-level TRIP directory exists
    TRIP = os.getenv('TRIP')
    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # build command including many variables not updated
    sim = os.path.join(TRIP,'iwave/trace/main/SEGYConv.x')
    if (adj==0):
        cmd = sim + \
        ' in=' + w + \
        ' out=' + d + \
        ' ker=' + g + \
        ' adj=0 >& crud'
    else:
        cmd = sim + \
        ' in=' + d + \
        ' out=' + w + \
        ' ker=' + g + \
        ' adj=1 >& crud'        
    # execute
    #    print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        print('command failed with return value ' + str(ret))
        return False

    return True

def cleanup():
    RSFROOT = os.getenv('RSFROOT')
    sfrm = os.path.join(RSFROOT,'bin/sfrm')
    dirlist = os.listdir()
    for file in dirlist:
        if file[len(file)-3:len(file)]=='.su':
            os.system('/bin/rm ' + file)
        if file[len(file)-4:len(file)]=='.rsf':
            os.system(sfrm + ' ' + file)
 

