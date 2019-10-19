import os
import sys
import array
import math

class InputError(Exception):
    """Exception raised for errors in the input.
    Attributes:
        inarg -- missing input argument key
        message -- explanation of the error
    """

    def __init__(self, inarg):
        self.inarg = inarg
        self.message = 'Missing input argument, '+str(inarg)+'\n'


if __name__ == '__main__':
    """ 
    This method reads in SU files and to compute approximate convergence rates
    via Richardson extrapolation.

    Arguments:
        infile -- input binary file
        it -- input time index for snapshot
        nt -- size of time traces
        ntr -- number of traces
        outfile -- name of output ASCII file, containing convergence rates 
                   computed using L2 and sup norms in time.
    """

    try:
        #reading input args
        it = []
        infile = []
        nt = []
        ntr = []
        outfile = []
        
        for tmp in sys.argv:
            if tmp.find("it=")==0:
                it = int(tmp[3:])
            if tmp.find("infile=")==0:
                infile = tmp[7:]
            if tmp.find("nt=")==0:
                nt = int(tmp[3:])
            if tmp.find("ntr=")==0:
                ntr = int(tmp[4:])
            if tmp.find("outfile=")==0:
                outfile = tmp[8:]

        #checking input args
        if it==[]:
            raise InputError("it")
        if infile==[]:
            raise InputError("infile")
        if nt==[]:
            raise InputError("nt")
        if ntr==[]:
            raise InputError("ntr")
        if outfile==[]:
            raise InputError("outfile")

        data = array.array('f')
        fp = open(infile,'rb')
        data.fromfile(fp,nt*ntr)
        fp.close()

        buff = []
        f = open(outfile,'w') 

        print('after having read file\n');

        for itr in range(ntr):
            buff = data[it+itr*nt]
            f.write(str(buff)+' ')
        f.write('\n')

        f.close()
    except InputError as err:
        print 'ERROR from snapshot.py :',err.message
        raise
