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

#import multiprocessing

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
        sim = os.path.join(TRIP,'iwave/trace/main/SEGYMConv.x')
        self.trmin = trmin
        self.trmax = trmax
        self.cmd = sim + ' min=' + str(trmin) + ' max=' + str(trmax)
        print('\ncmd = ' + self.cmd)
        
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
        
