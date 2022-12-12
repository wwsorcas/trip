import vcl
import os
import linalg
import tempfile

class RSFSpace(vcl.Space):

    def __init__(self,f):
        self.filename = f

    def getData(self):
        temp = tempfile.NamedTemporaryFile(delete=False,dir='/var/tmp',suffix='.su')
        temp.close()
        linalg.copy(self.filename,temp.name)
        linalg.scale(temp.name, 0.0)
        return temp.name

    def isData(self,x):
        return linalg.rsfcomp(self.filename,x)
    
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
        print('RSFSpace based on file ' + self.filename)
        
    
    
