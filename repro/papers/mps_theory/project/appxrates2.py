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
        in_h -- input SU file, solution with h spatial discretization
        in_2h -- input SU file, solution with 2h spatial discretization
        in_4h -- input SU file, solution with 4h spatial discretization
        nt -- size of time traces
        ntr -- number of traces
        outfile -- name of output ASCII file, containing convergence rates 
                   computed using L2 or sup norms in time.
    """

    try:
        #reading input args
        in_h  = []
        in_2h = []
        in_4h = []
        nt = []
        ntr = []
        outfile = []
        norm = []
        
        for tmp in sys.argv:
            if tmp.find("in_h=")==0:
                in_h = tmp[5:]
            if tmp.find("in_2h=")==0:
                in_2h = tmp[6:]
            if tmp.find("in_4h=")==0:
                in_4h = tmp[6:]
            if tmp.find("nt=")==0:
                nt = int(tmp[3:])
            if tmp.find("ntr=")==0:
                ntr = int(tmp[4:])
            if tmp.find("outfile=")==0:
                outfile = tmp[8:]
            if tmp.find("norm=")==0:
                norm = tmp[5:]

        #checking input args
        if in_h==[]:
            raise InputError("in_h")
        if in_2h==[]:
            raise InputError("in_2h")
        if in_4h==[]:
            raise InputError("in_4h")
        if nt==[]:
            raise InputError("nt")
        if ntr==[]:
            raise InputError("ntr")
        if outfile==[]:
            raise InputError("outfile")
        if norm==[]:
            raise InputError("norm")

        data_h = array.array('f')
        f_h = open(in_h,'rb')
        data_h.fromfile(f_h,nt*ntr)
        f_h.close()

        data_2h = array.array('f')
        f_2h = open(in_2h,'rb')
        data_2h.fromfile(f_2h,nt*ntr)
        f_2h.close()

        data_4h = array.array('f')
        f_4h = open(in_4h,'rb')
        data_4h.fromfile(f_4h,nt*ntr)
        f_4h.close()


        norm_top = []
        norm_bot = []
        f = open(outfile,'w') 

        #calculating rates via L2 norm
        if norm=='L2':
            for itr in range(ntr):

                norm_top = 0.0
                norm_bot = 0.0

                for it in range(nt):                
                    tmp = data_4h[it+itr*nt] - data_2h[it+itr*nt]
                    tmp = tmp*tmp
                    norm_top = norm_top + tmp
                    
                    tmp = data_2h[it+itr*nt] - data_h[it+itr*nt]
                    tmp = tmp*tmp
                    norm_bot = norm_bot + tmp
                    
                if norm_bot==0:
                    f.write('nan ')
                else:
                    rate = math.log(math.sqrt(norm_top/norm_bot))/math.log(2.0)
                    f.write(str(rate)+' ')        

        #calculating rates via sup norm
        if norm=='sup':
            for itr in range(ntr):

                norm_top = math.fabs( data_4h[itr*nt] - data_2h[itr*nt] )
                norm_bot = math.fabs( data_2h[itr*nt] - data_h[itr*nt] )

                for it in range(nt):       
         
                    tmp = math.fabs( data_4h[it+itr*nt] - data_2h[it+itr*nt] )
                    if tmp>norm_top:
                        norm_top = tmp
                        
                    tmp = math.fabs( data_2h[it+itr*nt] - data_h[it+itr*nt] )
                    if tmp>norm_bot:
                        norm_bot = tmp
                        
                    if norm_bot==0:
                        f.write('nan ')
                    else:
                        rate = math.log(norm_top/norm_bot)/math.log(2.0)
                        f.write(str(rate)+' ')        

        f.close()
    except InputError as err:
        print 'ERROR from appxrates2.py :',err.message
        raise


    








