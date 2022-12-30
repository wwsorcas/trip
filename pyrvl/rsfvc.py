import vcl
import os
from os import system as sys
import linalg
import tempfile

class Space(vcl.Space):

    def __init__(self,f):
        self.filename = f
        self.rsfroot = os.getenv('RSFROOT')
        self.rmcmd = os.path.join(self.rsfroot,'bin/sfrm')

    def getData(self):
        datapath = os.getenv('DATAPATH')
        temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.rsf')
        temp.close()
        linalg.copy(self.filename,temp.name)
        linalg.scale(temp.name, 0.0)
        return temp.name

    def isData(self,x):
        return linalg.rsfcomp(self.filename,x)
    
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
        try:
            # print('cleanup')
            notme = (x != self.filename)
            if not notme:
                raise Exception('cleanup called on space-defining object')
            #rrr = os.getenv('RSFROOT')
            #cmd = os.path.join(rsfroot,'bin/sfrm')
            ret = sys(self.rmcmd + ' ' + x)
            #ret = os.system(self.rmcmd + ' ' + x)
            if ret != 0:
                raise Exception('rsfvc.Space.cleanup: sfrm failed code=' \
                                    + str(ret))
        except Exception as ex:
            print(ex)
            raise Exception('called from rsfvc.Space.raw_cleanup')

    def raw_printData(self,x):
        print('Header file = ' + x)
        cmd = os.path.join(self.rsfroot,'bin/sfin')
        os.system(cmd + ' < ' + x)        
        # cmd = os.path.join(self.rsfroot,'bin/sfattr')
        # os.system(cmd + ' < ' + x)
        
    def myNameIs(self):
        print('RSF Space based on file ' + self.filename)
        cmd = os.path.join(self.rsfroot,'bin/sfin')
        os.system(cmd + ' < ' + self.filename)

        
    
    
