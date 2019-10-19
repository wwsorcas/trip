import os
import sys

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
    This method parses output from CG algorithm and spits out an ASCII file 
    with residual and gradient norms.

    Arguments:
        infile -- name of input CG file
        out_nres -- file name of output ASCII file with norm of residuals
        out_grad -- file name of output ASCII file with norm of gradient (normal residuals)
        out_xnrm -- file name of output ASCII file with norm of solution
        out_Lxnrm --
        out_xerr --
        out_Lxerr --
        Niter -- size of output vector containing residual and gradient norms.
                 If Niter > # iterations reported by CG file then the ASCII 
                 file will be padded by NaN.
    """

    try:
        #reading input args
        infile = []
        out_nres = []
        out_grad = []
        out_xnrm = []
        out_Lxnrm = []
        out_xerr = []
        out_Lxerr = []
        Niter = []
        
        for tmp in sys.argv:
            if tmp.find("infile=")==0:
                infile = tmp[7:]
            if tmp.find("out_nres=")==0:
                out_nres = tmp[9:]
            if tmp.find("out_grad=")==0:
                out_grad = tmp[9:]
            if tmp.find("out_xnrm=")==0:
                out_xnrm = tmp[9:]
            if tmp.find("out_Lxnrm=")==0:
                out_Lxnrm = tmp[10:]
            if tmp.find("out_xerr=")==0:
                out_xerr = tmp[9:]
            if tmp.find("out_Lxerr=")==0:
                out_Lxerr = tmp[10:]
            if tmp.find("Niter=")==0:
                Niter = int(tmp[6:])
                
        #checking input args
        if infile==[]:
            raise InputError("infile")
        if out_nres==[]:
            raise InputError("out_nres")
        if out_grad==[]:
            raise InputError("out_grad")
        if Niter==[]:
            raise InputError("Niter")

        #reading in CG file
        iter=[]
        nres=[]
        grad=[]
        xnrm=[]
        Lxnrm=[]
        xerr=[]
        Lxerr=[]        
        maxiter=[]

        with open(infile,'r') as f:

            buf = []
            
            #extracting max iterations
            for line in f:
                if line.find("* max iterations")==0:
                    buff = line.split()
                    maxiter = int(buff[4])
                    break

            if maxiter==[]:
                raise InputError("maxcount in CG file")
            
            #placing file pointer
            for line in f:
                if line.find("Iteration")==0:
                    break

            #reading in data
            for line in f:
               buf = line.split() 

               #cases where CG converged under maxiter
               if line.find("Gradient")==0:
                   break
               if line.find("Residual")==0:
                   break
               if line.find("-")==0:
                   continue

               iter = int(buf[0])
               nres.append(float(buf[1]))
               grad.append(float(buf[2]))

               if out_xnrm!=[]:
                   xnrm.append(float(buf[3]))
               if out_Lxnrm!=[]:
                   Lxnrm.append(float(buf[4]))
               if out_xerr!=[]:
                   xerr.append(float(buf[5]))
               if out_Lxerr!=[]:
                   Lxerr.append(float(buf[6]))

               if iter == maxiter:
                   break

        #writing out
        f = open(out_nres,'w')
        for i in range(0,Niter):
            if i<iter:
                f.write(str(nres[i])+' ')
            else:
                f.write('nan ')

        f = open(out_grad,'w')
        for i in range(0,Niter):
            if i<iter:
                f.write(str(grad[i])+' ')
            else:
                f.write('nan ')

        if out_xnrm!=[]:
            f = open(out_xnrm,'w')
            for i in range(0,Niter):
                if i<iter:
                    f.write(str(xnrm[i])+' ')
                else:
                    f.write('nan ')

        if out_Lxnrm!=[]:
            f = open(out_Lxnrm,'w')
            for i in range(0,Niter):
                if i<iter:
                    f.write(str(Lxnrm[i])+' ')
                else:
                    f.write('nan ')            

        if out_xerr!=[]:
            f = open(out_xerr,'w')
            for i in range(0,Niter):
                if i<iter:
                    f.write(str(xerr[i])+' ')
                else:
                    f.write('nan ')

        if out_Lxerr!=[]:
            f = open(out_Lxerr,'w')
            for i in range(0,Niter):
                if i<iter:
                    f.write(str(Lxerr[i])+' ')
                else:
                    f.write('nan ')            
        
    except InputError as err:
        print 'ERROR from parseCG.py :',err.message
        raise






