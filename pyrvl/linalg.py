import os
import m8r

TRIP = os.getenv('TRIP')
    
def sanity(name,suffix):

    # check that name is a string
    if not isinstance(name,str):
        print('linalg.sanity: first arg not a string')
        return False

    # check that suffix is a string
    if not isinstance(suffix,str):
        print('linalg.sanity: second arg not a string')
        return False    

    # check that name is valid path
    if not os.path.exists(name):
        print('linalg.sanity: path ' + name + ' not found\n')
        return False

    # check that name is long enough
    if len(name) < len(suffix)+2:
        print('linalg.sanity: filename ' + name + ' too short\n')
        return False
    
    # check that name has correct suffix
    if name[len(name)-len(suffix)-1:len(name)] != '.' + suffix:
        # print('linalg.sanity: filename ' + name + ' does not have correct suffix .' + suffix)
        return False

    return True

# copy vec1 over vec2
def copy(vec1,vec2):

    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        print(TRIP)
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYNorm.x
    if sanity(vec1,'su'):
        ret = os.system('/bin/cp ' + vec1 + ' ' + vec2)
        if ret != 0:
            print('SU copy failed with return value ' + str(ret))
            return False
        return True
    # if file is RSF, call grid/main/GridDot.x
    elif sanity(vec1,'rsf'): # and sanity(vec2,'rsf'):
        #rsfroot = os.getenv('RSFROOT')
        cmd = os.path.join(TRIP,'iwave/grid/main/GridCopy.x')
        ret = os.system(cmd + ' in=' + vec1 + ' out=' + vec2)
        #cmd = os.path.join(rsfroot,'bin/sfcp')
        #ret = os.system(cmd + ' ' + vec1 + ' ' + vec2)
        if ret != 0:
            print('RSF copy failed with return value ' + str(ret))
            return False
        return True
    else:
        print(vec1 + ', ' + vec2 + ' not su or rsf file')
        print('or incompatible')
        return False
    
def norm(vec):
    
    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYNorm.x
    if sanity(vec,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYNorm.x')
        ret = os.system(cmd + ' in=' + vec + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
        else:
            with open('jnk','r') as f:
                res = f.readline()
                return float(res[0:len(res)-1])
    # if file is RSF, call grid/main/GridNormt.x
    elif sanity(vec,'rsf'):
        cmd = os.path.join(TRIP,'iwave/grid/main/GridNorm.x')
        ret = os.system(cmd + ' in=' + vec + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
        else:
            with open('jnk','r') as f:
                res = f.readline()
                return float(res[0:len(res)-1])

    else:
        print(vec + ' not name of legal su or rsf file')
        return False

# returns dot product
def dot(vec1, vec2):

    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYDot.x
    if sanity(vec1,'su') and sanity(vec2,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYDot.x')
        ret = os.system(cmd + ' in1=' + vec1 + ' in2=' + vec2 + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in1=' + vec1 + ' in2=' + vec2 + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
        else:
            with open('jnk','r') as f:
                res = f.readline()
                return float(res[0:len(res)-1])
            
    # if files are RSF, call grid/main/GridDot.x
    elif sanity(vec1,'rsf') and sanity(vec2,'rsf'):
        cmd = os.path.join(TRIP,'iwave/grid/main/GridDot.x')
        ret = os.system(cmd + ' in1=' + vec1 + ' in2=' + vec2 + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in1=' + vec1 + ' in2=' + vec2 + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
        else:
            with open('jnk','r') as f:
                res = f.readline()
                return float(res[0:len(res)-1])

    else:
        print(vec1 + ' not name of legal su or rsf file, or')
        print(vec2 + ' not name of legal su or rsf file, or')
        print('different data structures')
        return False

# vec2 = a*vec1 + b*vec2
# with the default value of b, the classic axpy
def lincomb(a, vec1, vec2, b=1.0):

    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYDot.x
    if sanity(vec1,'su') and sanity(vec2,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYLinComb.x')
        ret = os.system(cmd + ' a=' + str(a) + ' in1=' + vec1 + \
                        ' b=' + str(b) + ' in2=' + vec2)
        if ret != 0:
            print('command:')
            print(cmd + ' in1=' + vec1 + ' in2=' + vec2 + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
            
    # if files are RSF, call grid/main/GridDot.x
    elif sanity(vec1,'rsf') and sanity(vec2,'rsf'):
        cmd = os.path.join(TRIP,'iwave/grid/main/GridLinComb.x')
        ret = os.system(cmd + ' a=' + str(a) + ' in1=' + vec1 + \
                        ' b=' + str(b) + ' in2=' + vec2)
        if ret != 0:
            print('command:')
            print(cmd + 'a=' + str(a) + ' in1=' + vec1 + \
                        ' b=' + str(b) + ' in2=' + vec2)
            print('failed with return value ' + str(ret))
            return False

    else:
        print(vec1 + ' not name of legal su or rsf file, or')
        print(vec2 + ' not name of legal su or rsf file, or')
        print('different data structures')
        return False

    return True

# a special case that occurs often enough to be worth
# singling out
def scale(vec, a):
    
    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYNorm.x
    if sanity(vec,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYScale.x')
        ret = os.system(cmd + ' in=' + vec + ' scale=' + str(a))
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec + ' scale=' + str(a))
            print('failed with return value ' + str(ret))
            return False
    # if file is RSF, call grid/main/GridNormt.x
    elif sanity(vec,'rsf'):
        cmd = os.path.join(TRIP,'iwave/grid/main/GridScale.x')
        ret = os.system(cmd + ' in=' + vec + ' scale=' + str(a))
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec + ' > jnk')
            print('failed with return value ' + str(ret))
            return False
    else:
        print(vec + ' not name of legal su or rsf file')
        return False
    return True

def rand(vec):

    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYRand.x
    if sanity(vec,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYRand.x')
        ret = os.system(cmd + ' in=' + vec + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec)
            print('failed with return value ' + str(ret))
            return False
    # if file is RSF, call grid/main/GridRand.x
    elif sanity(vec,'rsf'):
        cmd = os.path.join(TRIP,'iwave/grid/main/GridRand.x')
        ret = os.system(cmd + ' in=' + vec + ' > jnk')
        if ret != 0:
            print('command:')
            print(cmd + ' in=' + vec)
            print('failed with return value ' + str(ret))
            return False
    else:
        print(vec + ' not name of legal su or rsf file')
        return False

def hdrcomp(vec1, vec2):
    
    if not isinstance(TRIP,str):
        print('path to TRIP package not defined')
        return False

    if not os.path.exists(TRIP):
        print('TRIP package = ' + TRIP + ' not found')
        return False
    
    # if file is SU, call trace/main/SEGYDot.x
    if sanity(vec1,'su') and sanity(vec2,'su'):
        cmd = os.path.join(TRIP,'iwave/trace/main/SEGYCmpHdrs.x')
        ret = os.system(cmd + ' in1=' + vec1 + ' in2=' + vec2)
        if ret != 0:
            #print('command:')
            #print(cmd + ' in1=' + vec1 + ' in2=' + vec2)
            #print('failed with return value ' + str(ret))
            return False
    else:
        return False
    return True

def m8rint(x,y):
    if str(x) != "b''":
#        print('x = ' + str(x))
        return int(x)
    else:
        return y

def m8rfloat(x,y):
    if str(x) != "b''":
        return float(x)
    else:
        return y

def m8rstr(x):
    return str(x)[2:len(str(x))-3]

def rsfcomp1(vec1, vec2):
    return True

def rsfcomp(vec1, vec2, checkunit=True):
    ''' 
    returns true if headers of two rsf files are compatible
    else false. If checkunit=False, ignores data unit and label,
    so tests only grid geometry.
    '''

    if sanity(vec1,'rsf') and sanity(vec2,'rsf'):

        dtol=0.001

        n1=[1,1,1]
        d1=[1.0,1.0,1.0]
        o1=[0.0,0.0,0.0]
        unit1=['','','']
        n2=[1,1,1]
        d2=[1.0,1.0,1.0]
        o2=[0.0,0.0,0.0]
        unit2=['','','']
        

        inp1=m8r.Input(vec1)
        inp2=m8r.Input(vec2)

        n1[0]=m8rint(inp1.get('n1'),1)
        n1[1]=m8rint(inp1.get('n2'),1)
        n1[2]=m8rint(inp1.get('n3'),1)
        n2[0]=m8rint(inp2.get('n1'),1)
        n2[1]=m8rint(inp2.get('n2'),1)
        n2[2]=m8rint(inp2.get('n3'),1)

        d1[0]=m8rfloat(inp1.get('d1'),1.0)
        d1[1]=m8rfloat(inp1.get('d2'),1.0)
        d1[2]=m8rfloat(inp1.get('d3'),1.0)
        d2[0]=m8rfloat(inp2.get('d1'),1.0)
        d2[1]=m8rfloat(inp2.get('d2'),1.0)
        d2[2]=m8rfloat(inp2.get('d3'),1.0)
 
        o1[0]=m8rfloat(inp1.get('o1'),0.0)
        o1[1]=m8rfloat(inp1.get('o2'),0.0)
        o1[2]=m8rfloat(inp1.get('o3'),0.0)
        o2[0]=m8rfloat(inp2.get('o1'),0.0)
        o2[1]=m8rfloat(inp2.get('o2'),0.0)
        o2[2]=m8rfloat(inp2.get('o3'),0.0)

        unit1[0]=str(inp1.get('unit1'))
        unit1[1]=str(inp1.get('unit2'))
        unit1[2]=str(inp1.get('unit3'))
        unit2[0]=str(inp2.get('unit1'))
        unit2[1]=str(inp2.get('unit2'))
        unit2[2]=str(inp2.get('unit3'))
        
        for i in range(0,3):
#            print('n1[' + str(i) + '] = ' + str(n1[i]))
#            print('n2[' + str(i) + '] = ' + str(n2[i]))            
            if n1[i] != n2[i] or abs(d1[i]-d2[i])>dtol*d1[i] or abs(o1[i]-o2[i])>dtol*d1[i] or unit1[i] != unit2[i]:
                return False

        dataunit1=str(inp1.get('unit'))
        dataunit2=str(inp2.get('unit'))
        datalabel1=str(inp1.get('label'))
        datalabel2=str(inp2.get('label'))
        
        if checkunit and (dataunit1 != dataunit2 or datalabel1 != datalabel2):
            return False

        return True

    else:

        return False

def rsfboundstest(vec, u, l):

    if not sanity(vec,'rsf'):
        raise Exception('Error: rsfboundstest - input not rsf')

    inp=m8r.Input(vec)

    arr=inp.read()

    if np.min(arr) <= l:
        print('rsfboundstest: file = ' + vec.data)
        print('  lower bound violated')
        return False
    
    if np.max(arr) >= u:
        print('rsfboundstest: file = ' + vec.data)
        print('  upper bound violated')
        return False

    return True

def simplot(f, addcb=False, clip=None, minval=None, maxval=None, width=7, asprat=-1):

    '''
    simple plotting tool using matplotlib to plot su and rsf data.
    color scale controls use syntax from SU (ximage, psimage) and M8R (sfgrey)
    f = data filename
    addcb = colorbar switch
    clip = abs max color scale value for oscillatory data
    minval, maxval = min and max color scale values for signed data 
    '''

    RSFROOT = os.getenv('RSFROOT')
    if not os.path.exists(RSFROOT):
        print('simplot')
        print('root M8R directory RSFROOT = ' + RSFROOT + ' not found')
        print('M8R not correcly installed')
        return False
    
    if sanity(f,'su'):
        ff = f[0:len(f)-3]+'.rsf'
        cmd = os.path.join(RSFROOT,'bin/sfsuread')
        ret=os.system(cmd + ' read=data endian=0 < ' + f + ' > ' + ff)
        if ret != 0:
            print('simplot')            
            print('command:')
            print(cmd + ' read=data endian=0 < ' + f + ' > ' + ff) 
            print('failed with return value ' + str(ret))
            return False
        asprat=1
     
    elif sanity(f,'rsf'):
        ff=f
        
    else:
        print('simplot: arg does not name file will either .su or .rsf suffix')
        return False

    # plot

    import m8r
    import numpy as np
    import matplotlib.pyplot as plt

    inp=m8r.Input(ff)
    n1=int(inp.get('n1'))
    o1=float(inp.get('o1'))
    d1=float(inp.get('d1'))
    n2=int(inp.get('n2'))
    o2=float(inp.get('o2'))
    d2=float(inp.get('d2'))
    datafile=str(inp.get('in'))
    # print('datafile=' + datafile[2:len(datafile)-3])
    data=inp.read()

    # start with the clip
    #if clip is not None:
    #    np.clip(a=data, a_min=-clip, a_max=clip, out=data)

    # line plot
    if n2==1:
        
        Z = np.linspace(o1, o1+(n1-1)*d1, n1)
        fig, ax = plt.subplots(figsize=(width,width))
        ax.set_title(f)
        plt.plot(Z,data)
        plt.show()

    # color scale
    else:

        # define grid
        X, Y = np.meshgrid(np.linspace(o2, o2+(n2-1)*d2, n2), np.linspace(o1, o1+(n1-1)*d1, n1))

        # if nonsense asprat has not been overridden, override it now
        if asprat < 0:
           asprat = n1*d1/(n2*d2)
       
        fig, ax = plt.subplots(figsize=(width,asprat*width))

        # if data is of one sign, set color extremes to min and max. Otherwise, 
        # use symmetric interval about zero

        if np.min(data) < 0 and np.max(data) > 0:
            if clip is not None:
                colmax = clip
            else:
                colmax = max(abs(np.min(data)),np.max(data))
            colmin = -colmax
        else:
            if minval is not None:
                colmin = minval
            else:
                colmin = np.min(data)
            if maxval is not None:
                colmax = maxval
            else:
                colmax = np.max(data)

        print('colmax=' + str(colmax) + ' colmin=' + str(colmin))
            
        pc = ax.pcolormesh(X, Y, (data.T)[0:n1-1,0:n2-1], vmin=colmin, vmax=colmax, cmap='RdBu_r')
        if addcb:
            fig.colorbar(pc, ax=ax)
              
        ax.set_title(f)

        # y axis increases downwards, seismic-style
        plt.ylim(max(plt.ylim()), min(plt.ylim()))
        plt.show()
        print('simplot: data min = %10.4e, data max = %10.4e' % (np.min(data), np.max(data)))
        if clip is not None:
            print('simplot: clip = %10.4e' %(clip))

    # cleanup - remove temp rsf file if necessary
        # print('ff=' + ff)
        if sanity(f,'su'):
            os.unlink(ff)
            os.unlink(datafile[2:len(datafile)-3])

import m8r
import numpy as np
import array

def ndarraytorsfdata(x, rsfname):
    '''
    Overwrite NumPy ndarray onto data file of rsf pair. ndarray size
    must be same as data file size, but axis lengths are not checked.
    Note that original rsf data is destroyed. NumPy data is implicitly 
    dimensionless, so data inherits unit of rsf. Data written as 32 bit 
    IEEE floats, native endian.

    Parameters:
        x (ndarray): input NumPy data
        rsfname (string): name of output rsf header file
        returns None
    '''
       
    try:
        if not isinstance(rsfname,str):
            raise Exception('Error: input second argument not string')
        if not sanity(rsfname,'rsf'):
            raise Exception('Error: input file name does not have .rsf suffix')
        if not isinstance(x,np.ndarray):
            raise Exception('Error: input first argument not ndarray')
        inp=m8r.Input(rsfname)
        n1=int(inp.get('n1'))
        n2=int(inp.get('n2'))
        n3=int(inp.get('n3'))
        if x.size != n1*n2*n3:
            raise Exception('Error: ndarray of size = ' + str(x.size) +
                                ', rsf array size = ' + str(n1*n2*n3))
        dataname = m8rstr(inp.get('in'))
        datafile = open(dataname,'wb')
        # allocate single precision array
        dataarray = array.array('f')
        # read in list data, converting to single precision
        dataarray.fromlist((np.ravel(x)).tolist())
        # output signle precision data to file
        dataarray.tofile(datafile)
        datafile.close()

    except Exception as ex:
        print(ex)
        raise Exception('called from function ndarraytorsfdata')

def ndarraytorsf(x, rsfname, n=None, d=None, o=None, unit=None, units=None):
    '''
    Write NumPy ndarray and aux data onto rsf pair. ndarray size
    must be same as data file size, but axis lengths are not checked.
    Data written as 32 bit IEEE floats, native endian.

    Parameters:
        x (ndarray): input NumPy data
        rsfname (string): name of output rsf header file
        n (int array or None): axis lengths
        d (float array or None): axis increments
        o (float array or None): axis origins
        unit (string): quantity unit
        units... (string array): axis units
        returns None
    '''
       
    try:
        dim = 1
        nout = []
        oout = []
        dout = []
        if n is None:
            nout.append(x.size)
        else:
            dim = len(n)
            tot = 1
            for i in range(dim):
                nout.append(n[i])
                tot *= n[i]
            if tot != x.size:
                raise Exception('product of axis lengths != x.size')
                
        if d is None:
            for i in range(dim):
                dout.append(1.0)
        else:
            if len(d) != dim:
                raise Exception('d.size != n.size')
            for i in range(dim):
                dout.append(d[i])
        
        if o is None:
            for i in range(dim):
                oout.append(0.0)
        else:
            if len(o) != dim:
                raise Exception('o.size != n.size')
            for i in range(dim):
                oout.append(o[i])

        # write header file
        hdr = open(rsfname, 'w')
        for i in range(dim):
            hdr.write('n' + str(i+1) + '=' + str(nout[i]) +
                          ' d' + str(i+1) + '=' + str(dout[i]) +
                          ' o' + str(i+1) + '=' + str(oout[i]) + '\n')
        if unit is not None:
            hdr.write('unit=' + unit)
        if units is not None:
            if len(units) != dim:
                raise Exception('units.size != n.size')
            for i in range(dim):
                hdr.write(' unit' + str(i+1) + '=' + units[i])
        if unit is not None or units is not None:
            hdr.write('\n')
        hdr.write('data_format="native_float" esize=4 \n')
        hdr.write('in=' + rsfname + '@')
        hdr.close()

        # write data file
        datafile = open(rsfname + '@', 'wb')
        # allocate single precision array
        dataarray = array.array('f')
        # read in list data, converting to single precision
        dataarray.fromlist((np.ravel(x)).tolist())
        # output signle precision data to file
        dataarray.tofile(datafile)
        datafile.close()
            
    except Exception as ex:
        print(ex)
        raise Exception('called from function ndarraytorsf for file ' + rsfname)


def ndarrayfromrsfdata(rsfname):
    '''
    reads rsf data into NumPy ndarray. Returned ndarray has axis 
    lengths specified by rsf. Note that NumPy data is implicitly 
    dimensionless, so unit and axis unit info is lost.
   
    Parameters:
        rsfname (string): name of rsf header file
        returns NumPy ndarray
    '''
    try:
        if not sanity(rsfname,'rsf'):
            raise Exception('Error: input file name does not have .rsf suffix')
        inp=m8r.Input(rsfname)
        data = inp.read()
        return data
    except Exception as ex:
        print(ex)
        raise Exception('called from function ndarrayfromrsfdata')

