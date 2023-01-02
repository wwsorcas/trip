import vcl
import os
#from os import system as sys
import linalg
import tempfile

class Space(vcl.Space):

    def __init__(self,f):
        self.filename = f
        self.rsfroot = os.getenv('RSFROOT')
        self.rmcmd = os.path.join(self.rsfroot,'bin/sfrm')
        self._unlink = os.unlink

    def getData(self):
        datapath = os.getenv('DATAPATH')
        temp = tempfile.NamedTemporaryFile(delete=False,dir=datapath,suffix='.rsf')
        temp.close()
        filename = temp.name
        # unlink because copy tries to use data from existing target file
        os.unlink(temp.name)
        linalg.copy(self.filename,filename)
        linalg.scale(filename, 0.0)
        return filename

    def isData(self,x):
        return ((self.filename == x) or linalg.rsfcomp(self.filename,x))
    
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
        print('RSF CLEANUP file = ' + x)
        try:
            notme = (x != self.filename)
            if not notme:
                raise Exception('cleanup called on space-defining object')
            #rsfroot = os.getenv('RSFROOT')
            #cmd = os.path.join(rsfroot,'bin/sfrm')
            #ret = os.system(cmd + ' ' + x)
            #if ret != 0:
            #    raise Exception('rsfvc.Space.cleanup: sfrm failed code=' \
            #                        + str(ret))
            self._unlink(x)
            self._unlink(x + '@')
        except Exception as ex:
            print(ex)
            raise Exception('called from rsfvc.Space.cleanup')

    def raw_printData(self,x):
        print('RSF Header file = ' + x)
        cmd = os.path.join(self.rsfroot,'bin/sfin')
        os.system(cmd + ' < ' + x)        
        # cmd = os.path.join(self.rsfroot,'bin/sfattr')
        # os.system(cmd + ' < ' + x)
        
    def myNameIs(self):
        print('RSF Space based on file ' + self.filename)
        #cmd = os.path.join(self.rsfroot,'bin/sfin')
        #os.system(cmd + ' < ' + self.filename)

        
    
    
