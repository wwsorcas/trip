import vcl
import os
import linalg
import tempfile

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
        return linalg.rsfcomp(self.filename,x)
    
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
    def raw_copy(self,x,y)
        linalg.copy(x,y)

    # for use in vector destructor - x is data 
    def raw_cleanup(self,x):
        if self.isData(x) and self.filename != x:
            os.unlink(x)

    def raw_printData(self,x):
            print('RSF file with name = ' + x)
       
    def myNameIs(self):
        print('RSF Space based on file ' + self.filename)
        
    
    
