import vcl
import os
import linalg
import tempfile
import op

class Space(vcl.Space):

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
    def raw_linComb(self,a,x,y,b=1.0):
        linalg.lincomb(a,x,y,b)

    # operates on data = filenames        
    def raw_dot(self,x,y):
        return linalg.dot(x,y)
     
    # operates on data = filenames            
    def raw_zero(self,x):
        linalg.scale(x,0.0)

    # convenience functions
    def raw_copy(self,x,y):
        linalg.copy(x,y)

    # for use in vector destructor - x is data 
    def raw_cleanup(self,x):
        if self.filename != x:
            os.unlink(x)

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
    
    def raw_applyFwd(self,x, y):
        op.convop(self.g,x.data,y.data,adj=0)
            
    def raw_applyAdj(self,x, y):
        op.convop(self.g,y.data,x.data,adj=1)

    def myNameIs(self):
        print('Convolution operator with kernel ' + self.g)
        print('domain:')
        self.dom.myNameIs()
        print('range:')
        self.rng.myNameIs()

