import vcl
import os
import linalg
import tempfile
import op
import sys

class Space(vcl.Space):

    def __init__(self,f):
        self.filename = f
        self._unlink = os.unlink

    def getData(self):
        try: 
            datapath = os.getenv('DATAPATH')
            if not os.path.exists(datapath):
                raise Exception('Error: datapath = ' + datapath + ' not valid path')
            temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.su')
            temp.close()
            if not linalg.copy(self.filename,temp.name):
                raise Exception('Error: linalg.copy')
            if not linalg.scale(temp.name, 0.0):
                raise Exception('Error: linalg.scale')
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.Space.getData')
        else:
            return temp.name

    def isData(self,x):
        return linalg.hdrcomp(self.filename,x)

    def saveto(self,x,y):
        try:
            if not self.isData(x):
                raise Exception('input data ' + str(x) + ' not rsf')
            if not linalg.copy(x,y):
                raise Exception('output data ' + str(y) + ' not legal path')
        except Exception as ex:
            print(ex)
            raise Exception('called from rsfvc.Space.saveto')
    
    # operates on data = filenames
    def raw_linComb(self,a,x,y,b=1.0):
        linalg.lincomb(a,x,y,b)
        return y

    # operates on data = filenames        
    def raw_dot(self,x,y):
        return linalg.dot(x,y)
     
    # operates on data = filenames            
    def raw_scale(self,x,c):
        linalg.scale(x,c)
        return x

    # convenience functions
    def raw_copy(self,x,y):
        linalg.copy(x,y)
        return y
    
    # for use in vector destructor - x is data 
    def cleanup(self,x):
        if self.filename != x:
            # print('SEGY CLEANUP file = ' + x)
            self._unlink(x)

    def raw_printData(self,x):
        print('file with name = ' + x)
       
    def myNameIs(self):
        print('SEGYSpace based on SEGY file ' + self.filename)
            
class ConvolutionOperator(vcl.LinearOperator):

    # dom and rng are SEGYSpaces, green is a filename
    def __init__(self,dom,rng,green):        
        self.dom = dom
        self.rng = rng
        self.g   = green 
    
    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.rng
    
    def applyFwd(self,x, y):
        try:
            if not op.convop(self.g,x.data,y.data,adj=0):
                raise Exception('Error: op.convop call failed')
        except Exception as ex:
            print(ex)
            print('\nkernel = ' + self.g)
            print('\nInput vector =')
            x.myNameIs()
            print('\nOutput vector =')
            y.myNameIs()
            raise Exception('called from segyvc.OperatorConvolution.applyFwd')
            
    def applyAdj(self,x, y):
        try:
            if not op.convop(self.g,y.data,x.data,adj=1):
                raise Exception('Error: op.convop call failed')
        except Exception as ex:
            print(ex)
            print('\nkernel = ' + self.g)
            print('\nInput vector =')
            x.myNameIs()
            print('\nOutput vector =')
            y.myNameIs()
            raise Exception('called from segyvc.OperatorConvolution.applyAdj')
        
    def myNameIs(self):
        print('Convolution operator with kernel ' + self.g)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()

        
class MultiTraceConvOp(vcl.LinearOperator):

    '''
    trace-by-trace convolution of two SEGY data sets (input and kernel), 
    output to another SEGY data set. Input and output represented as vcl.Vectors
    in SEGYSpaces, kernel is given as a filename (can be data member of another
    SEGYSpace vcl.Vector). Input and output must have at least as many traces as 
    kernel. Only the trace range between trmin and trmax is convolved and written
    to the corresponding traces in output. Both convolution and adjoint convolution
    (cross-correlation) are represented as fwd and adj options of this linear op. 
    Note that since possibly only part of the output data is altered, this is really
    a AXPY, rather than a linear op, unless all traces are written.

    Constructor Parameters:
    dom (SEGYSpace) - domain (inputs)
    rng (SEGYSpace) - range (outputs)
    ker (string) - filename for kernel traces
    trmin (int) - min trace index (tracl)
    trmax (int) - max trace index
    '''
    
    def __init__(self,dom,rng,ker,trmin=0,trmax=sys.maxsize):        
        self.dom = dom
        self.rng = rng
        if not linalg.sanity(ker,'su'):
            raise Exception('Error: convolution kernel file not SEGY') 
        self.ker   = ker
        TRIP = os.getenv('TRIP')
        if not os.path.exists(TRIP):
            raise Exception('TRIP package = ' + TRIP + ' not found')
        sim = os.path.join(TRIP,'iwave/trace/main/SEGYConv.x')
        self.trmin = trmin
        self.trmax = trmax
        self.cmd = sim + ' min=' + str(trmin) + ' max=' + str(trmax)
        # print('\ncmd = ' + self.cmd)
        
    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.rng
    
    def applyFwd(self,x, y):
        try:
            loccmd = self.cmd + ' in=' + x.data + \
              ' out=' + y.data + ' ker=' + self.ker + \
              ' adj=0'
            print('fwd: loccmd = ' + loccmd)
            ret= os.system(loccmd)
            if ret != 0: 
                msg = 'Failed to execute ' + loccmd + \
                  '\ncalled from segyvc.MultiTraceConvOp.applyFwd'
                raise Exception(msg)
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.MultiTraceConvOp.applyFwd')
        
    def applyAdj(self,x, y):
        try:
            loccmd = self.cmd + ' in=' + x.data + \
              ' out=' + y.data + ' ker=' + self.ker + \
              ' adj=1'
            print('adj: loccmd = ' + loccmd)
            ret= os.system(loccmd)
            if ret != 0: 
                msg = 'Failed to execute ' + loccmd + \
                  '\ncalled from segyvc.MultiTraceConvOp.applyAdj'
                raise Exception(msg)
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.MultiTraceConvOp.applyFwd')
        
    def myNameIs(self):
        print('Multi Trace Convolution operator with kernel ' + self.ker)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()
        print('trace range [' + str(self.trmin) + ', ' + str(self.trmax-1) + ']')
        
class TScaleOperator(vcl.LinearOperator):

    def __init__(dom,scale):
        self.dom = dom
        self.scale = scale

    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.dom   

    def applyFwd(self,x, y):
        try:
            cwproot = os.getenv('CWPROOT')
            if not isinstance(cwproot,str):
                raise Exception('Error: path CWPROOT not defined')
            if not os.path.exists(cwproot):
                raise Exception('Error: CWPROOT = ' + cwproot + ' not found')
            cmd = os.path.join(cwproot,'bin/sugain')
            ret = os.system(cmd + ' scale=' + str(self.scale) + ' tpow=1.0 ' +
                                '< ' + x.data + ' > '  + y.data)
            if ret != 0:
                raise Exception('Error: call to sugain failed')
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.TScaleOperator.applyFwd')
        
    def applyAdj(self,x,y):
        try:
            self.applyFwd(x,y)
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.TScaleOperator.applyAdj')

    def myNameIs(self):
        print('segyvc.TScaleOperator')
        print('  uses sugain to multiply traces by scale*t')
        ptinr('  scale = ' + str(self.xcale))
        print('  domain = rainge = ')
        self.dom.myNameIs()

################################# parconv ###############################

# parallel convolution and supporting functions

import multiprocessing

def getntr(name):
    cwproot = os.getenv('CWPROOT')
    cmd = os.path.join(cwproot,'bin/sucountkey')
    ret = os.system(cmd + ' < ' + name + ' key=tracl verbose=0 | sed s/"tracl"// > jnk')
    if ret != 0:
        raise Exception('call to sucountkey failed')
        
    with open('jnk','r') as fn:
        res = fn.readline()
        return int(res[0:len(res)-1])


def divvy(ntr=2,partask=1):
    '''
    returns array of length partask listing task intervals for
    divistion of ntr tasks into partask chunks.
    '''
    try:
        if not isinstance(ntr,int) or not isinstance(partask,int):
            raise Exception('inputs must be integers')
        if not ntr > 0 or not partask > 0:
            raise Exception('inputs must be positive')
        if partask > ntr:
            raise Exception('number of jobs (ntr) must be at lease number of tasks (partask)')
        
        chunk=int(ntr/partask)
        if ntr/partask > chunk:
            chunk += 1
        bot = 0
        top = chunk-1
        work = []
        while bot < ntr:
            top = min(bot + chunk - 1, ntr-1)
            work.append([bot, top])
            bot += chunk
        return work
    except Exception as ex:
        print(ex)
        raise Exceeption('called from segyvc:divvy')

import time
import psutil

def doit(cmd):
    ''' 
    ad-hoc version of exec
    '''
    try:
        print('doit enter pid=' + str(os.getpid()))
        rpt = 100
        for i in range(rpt):
            ret = os.system(cmd)
        if ret != 0: 
            raise Exception('Failed to execute ' + cmd + '\nos.system=' + str(ret))
        time.sleep(5)
        print ('doit exit pid=' + str(os.getpid()))
    except Exception as ex:
        print(ex)
        raise Exception('called from segyvc.doit')
    

class ParConvOp(vcl.LinearOperator):

    '''
    trace-by-trace convolution of two SEGY data sets (input and kernel), 
    output to another SEGY data set. Input and output represented as vcl.Vectors
    in SEGYSpaces, kernel is given as a filename (can be data member of another
    SEGYSpace vcl.Vector). Input and output must have at least as many traces as 
    kernel. Only the trace range between trmin and trmax is convolved and written
    to the corresponding traces in output. Both convolution and adjoint convolution
    (cross-correlation) are represented as fwd and adj options of this linear op. 
    Note that since possibly only part of the output data is altered, this is really
    a AXPY, rather than a linear op, unless all traces are written.

    Constructor Parameters:
    dom (SEGYSpace) - domain (inputs)
    rng (SEGYSpace) - range (outputs)
    ker (string) - filename for kernel traces
    partask (int) - number of processes - 0 = serial (same as 1)
    '''
    
    def __init__(self,dom,rng,ker,partask=0):
        try:
            self.dom = dom
            self.rng = rng
            self.ker = ker
            self.partask = partask

            # sanity checks
            if not linalg.sanity(ker,'su'):
                raise Exception('convolution kernel file not SEGY')
            self.ntr = getntr(ker)
            if not isinstance(dom,Space):
                raise Exception('domain not segyvc.Space')
            if self.ntr != getntr(dom.filename):
                raise Exception('domain space, kernel do not have same number of traces')
            if not isinstance(rng,Space):
                raise Exception('range not segyvc.Space')
            if self.ntr != getntr(rng.filename):
                raise Exception('range space, kernel do not have same number of traces')            
            TRIP = os.getenv('TRIP')
            if not os.path.exists(TRIP):
                raise Exception('TRIP package = ' + TRIP + ' not found')
            self.mconv = os.path.join(TRIP,'iwave/trace/main/SEGYMConv.x')
            self.tasklist = divvy(self.ntr, self.partask)
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.ParConvOp constructor')
        
    def getDomain(self):
        return self.dom
    
    def getRange(self):
        return self.rng

    def applyFwd(self,x, y):
        try:
            taskarray = []
            for k in range(len(self.tasklist)):
                taskarray.append(self.mconv + ' in=' + x.data + \
              ' out=' + y.data + ' ker=' + self.ker + \
              ' adj=0' + ' trmin=' + str(self.tasklist[k][0]) + \
              ' trmax=' + str(self.tasklist[k][1]))

            print('taskarray:')
            print(taskarray)

            pool = multiprocessing.Pool(self.partask)
            pool.map(doit, taskarray)
            
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.ParConvOp.applyFwd')
        
    def applyAdj(self,x, y):
        try:
            taskarray = []
            for k in range(len(self.tasklist)):
                taskarray.append(self.mconv + ' in=' + x.data + \
              ' out=' + y.data + ' ker=' + self.ker + \
              ' adj=1 trmin=' + str(self.tasklist[k][0]) + \
              ' trmax=' + str(self.tasklist[k][1]))

            pool = multiprocessing.Pool(self.partask)
            pool.map(doit, taskarray)
    
        except Exception as ex:
            print(ex)
            raise Exception('called from segyvc.MultiTraceConvOp.applyFwd')
        
    def myNameIs(self):
        print('Parallel Convolution operator with kernel ' + self.ker)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()
        print('number of tasks = ' + str(self.partask))
        
