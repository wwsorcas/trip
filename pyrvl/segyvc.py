import vcl
import os
import linalg
import tempfile
import op

class SEGYSpace(vcl.Space):

    def __init__(self,f):
        self.filename = f

    def getData(self):
        temp = tempfile.NamedTemporaryFile(delete=False,dir='/var/tmp',suffix='.su')
        temp.close()
        linalg.copy(self.filename,temp.name)
        linalg.scale(temp.name, 0.0)
        return temp.name

    def isData(self,x):
        return linalg.hdrcomp(self.filename,x)
    
    # operates on data = filenames
    def linComb(self,a,x,y,b=1.0):
        if self.isData(x) and self.isData(y):
            linalg.lincomb(a,x,y,b)
            return True
        return False

    # operates on data = filenames        
    def dot(self,x,y):
        if self.isData(x) and self.isData(y):
            return linalg.dot(x,y)
        return False
    
    # operates on data = filenames            
    def zero(self,x):
        if self.isData(x):
            linalg.scale(x,0.0)

    # convenience functions
    def copy(self,x,y):
        if self.isData(x) and self.isData(y):        
            linalg.copy(x,y)
            return True
        return False

    # for use in vector destructor - x is data 
    def cleanup(self,x):
        if self.isData(x) and self.filename != x:
            os.unlink(x)

    def printData(self,x):
        if self.isData(x):        
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
        # first check that x in dom, y in range
        if x.space != self.dom:
            print('ConvOp.applyFwd: input not in domain')
            return False
        if y.space != self.rng:
            print('ConvOp.applyFwd: output not in range')
            return False
        return op.convop(self.g,x.data,y.data,adj=0)
            
    def applyAdj(self,x, y):
        # first check that x in range, y in domain
        if x.space != self.rng:
            print('ConvOp.applyAdj: input not in range')
            return False
        if y.space != self.dom:
            print('ConvOp.applyAdj: output not in domain')
            print('output:')
            y.space.myNameIs()
            print('domain:')
            self.dom.myNameIs()
            return False
        return op.convop(self.g,y.data,x.data,adj=1)

    def myNameIs(self):
        print('Convolution operator with kernel ' + g)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()

