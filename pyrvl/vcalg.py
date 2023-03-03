import gc
import os
import vcl
import math
from vcl import Vector
from vcl import transp

# single steps are separated from loops as separate functions in order
# to force garbage-collection.

###################### CONJUGATE GRADIENT ITERATION ######################

# residual and normal residual vectors e and r are arguments, so
# that they can be examined after completion. Other work vectors are
# local, deleted on return.

# args:
#  x = estimated solution on return
#  b = data
#  A = operator
#  kmax = max iterations for termination
#  eps = relative residual decrease for termination
#  rho = relative normal residual decrease for termination
#  e = residual on return (optional)
#  r = normal residual on return (optional)
#  verbose = verbosity level
#    0 - no print output
#    1 - print number of its performed, residual and normal residual
#    2 - (or more) print it, res, nres at every it

def cgstep(x, A, k, p, r, e, rnorm, enorm, gamma, verbose):
#    os.system('echo top; ls /var/tmp')
    try:
        # print('k='+str(k))

        # print('#1. $q = Ap$')
        q = A*p
            
        # print('#2. $s = A^Tq$')
        s = transp(A)*q
            
        #os.system('echo q s; ls /var/tmp')
            
        # print('#3. $\alpha = gamma / \langle p, s\rangle$')
        alpha = gamma/p.dot(s)

        # print('#4. $x \leftarrow x+\alpha p$')
        x.linComb(alpha,p)

        # print('#5. $e \leftarrow e-\alpha q$')
        e.linComb(-alpha,q)
        enorm = e.norm()
        
        # print('#6. $r \leftarrow r-\alpha s$')
        r.linComb(-alpha,s)
        
        # print('#7. $\delta = \langle r, r \rangle$')  
        delta = r.dot(r)
        rnorm = math.sqrt(delta)
        
        # print('#8. $\beta = \delta / \gamma$')
        beta = delta/gamma
        
        # print('#9. $p \leftarrow r + \beta p$')
        p.linComb(1.0,r,beta)
        
        # print('#10. $\gamma \leftarrow delta$')
        gamma=delta
        
        # print('#11. $k \leftarrow k+1$')
        k=k+1
        
        # print('#12. print $k$ (iteration) $\|e\|$ (residual norm),
        # $\|r\|$ (normal residual norm)')
        if verbose > 1:
            print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))

        # without the following code, the vector destructor is not called
        # until the program exits
        # reduce ref count for q, s to 1, which is unreachable - the
        # reference created by invoking the constructor 
        del q
        del s
        # force garbage collection
        n=gc.collect()

        return [k, enorm, rnorm, gamma]

    except Exception as ex:
        print(ex)
        raise Exception('called from cgstep')
                
def conjgrad(x, b, A, kmax, eps, rho, verbose=0, e=None, r=None):

    try:

        #Initialize:
        # print('#1. $x = 0$')
        x.scale(0.0)
        
        # print('#2. $e = b$')
        if e is None:
            e = b.dup()
        else:
            e.copy(b)
        
        # print('#3. $r = A^Tb$')
        p = transp(A)*b
        
        # print('#4. $p = r$')
        if r is None:
            r = p.dup()
        else:
            r.copy(p)
        
        # print('#7. $\gamma_0 = \langle r, r \rangle$')
        gamma0 = r.dot(r)
        
        # print('#8. $\gamma = \gamma_0$')
        gamma=gamma0
        
        # print('#9. $k=0$')
        k=0
        enorm0=e.norm()
        rnorm0=r.norm()
        enorm=enorm0
        rnorm=rnorm0
    
        if verbose > 1:
            print('  k       |e|       |r|=')
            print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))
    
        #Repeat while $k<k_{\rm max}$, $\|e\|>\epsilon \|d\|$:
        while k<kmax and enorm>eps*enorm0 and rnorm>rho*rnorm0:

            [k, enorm, rnorm, gamma] = cgstep(x, A, k, p, r, e, rnorm, enorm, gamma, verbose)
            
        if verbose > 0:
            print('----------------------------------------------------')
            print('  k       |e|     |e|/|e0|      |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e  %10.4e  %10.4e' % (k, enorm, enorm/enorm0, rnorm, rnorm/rnorm0))

    except Exception as ex:
        print(ex)
        raise Exception("called from cg.conjgrad")
        

########### TRUST-RADIUS MODIFIED CONJUGATE GRADIENT ITERATION ##########
def trcgstep(x, A, k, p, r, rnorm, rnorm0, rho, gamma, Delta, cgend, verbose):

            # print('k='+str(k))

            # print('#1. $q = Ap$')
            # A.applyFwd(p,q)
            q = A*p

            # print('#2. $s = A^Tq$')
            # A.applyAdj(q,s)
            s = transp(A)*q

            # print('#3. $\alpha = gamma / \langle p, s\rangle$')
            alpha = gamma/p.dot(s)

            # print('#4. $x \leftarrow x+\alpha p$')
            x.linComb(alpha,p)
            xnorm = x.norm()
            
            if xnorm > Delta:

                x.scale(Delta/xnorm)
                cgend = True
                if verbose > 0:
                    print('trust radius exceeded - return scaled ' + \
                              'soln of length ' + str(Delta))

            else:
                
                # print('#6. $r \leftarrow r-\alpha s$')
                r.linComb(-alpha,s)
                
                # print('#7. $\delta = \langle r, r \rangle$')  
                delta = r.dot(r)
                rnorm = math.sqrt(delta)

                if rnorm>rho*rnorm0:
                    
                    # print('#8. $\beta = \delta / \gamma$')
                    beta = delta/gamma
                
                    # print('#9. $p \leftarrow r + \beta p$')
                    p.linComb(1.0,r,beta)
                
                    # print('#10. $\gamma \leftarrow delta$')
                    gamma=delta
                
                # print('#11. $k \leftarrow k+1$')
                k=k+1

            # print('#12. print $k$ (iteration) $\|e\|$ (residual norm),
            # $\|r\|$ (normal residual norm)')
                if verbose > 1:
                    print('%3d  %10.4e' % (k, rnorm))

            del q
            del s
            # force garbage collection
            n=gc.collect()

            return [x, k, p, r, rnorm, gamma, cgend]
    
# residual and normal residual vectors e and r are internal
# "normal" termination when normal residual falls below tolerance
# so no need to test residual

# args:
#  x = estimated solution on return
#  b = data
#  A = operator
#  kmax = max iterations for termination
#  rho = relative normal residual decrease for termination
#  Delta = trust radius
#  verbose = verbosity level
#    0 - no print output
#    1 - print number of its performed, residual and normal residual
#    2 - (or more) print it, res, nres at every it

def trconjgrad(x, b, A, kmax, rho, Delta, verbose=0):

    try:

        #Initialize:
        # print('#1. $x = 0$')
        x.scale(0.0)
        
        # print('#3. $r = A^Tb$')
        r = transp(A)*b
        
        # print('#4. $p = r$')
        p = r.dup()

        # print('#5. $\gamma_0 = \langle r, r \rangle$')
        gamma0 = r.dot(r)
        
        # print('#6. $\gamma = \gamma_0$')
        gamma=gamma0
        
        # print('#7. $k=0$')
        k=0
        rnorm0=r.norm()
        rnorm=rnorm0

        xnorm = x.norm()
    
        if verbose > 1:
            print('  k       |r|')
            print('%3d  %10.4e' % (k, rnorm))    

        # flag to indicate trust region truncation, indicating
        # return to gn loop
        cgend = False
        
        #Repeat while $k<k_{\rm max}$, rnorm > rho*rnorm0, |x|<Delta
        while k<kmax and rnorm>rho*rnorm0 and not cgend:
            
            [x, k, p, r, rnorm, gamma, cgend] = trcgstep(x, A, k, p, r, rnorm, rnorm0, rho, gamma, Delta, cgend, verbose)
            
        if verbose > 0:
            print('-------------------->>> CG >>>------------------------')
            print('  k       |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e' % (k, rnorm, rnorm/rnorm0))
            print('--------------------<<< CG <<<------------------------')

    except Exception as ex:
        print(ex)
        raise Exception("called from cg.conjgrad")

    else:
        return k


########################## TRUST REGION GAUSS-NEWTON #######################

# arguments:
#  x        = solution estimate, Vector in domain of F
#  b        = rhs, Vector in range of F
#  F        = Function
#  imax     = max GN iterates
#  eps      = gradient reduction - GN stopping criterion
#  kmax     = max CG iterates per GN iteration
#  rho      = normal residual reduction - CG stopping criterion
#  Delta    = trust radius
#  mured    = Delta reduction factor
#  muinc    = Delta increase factor
#  gammared = G-A reduction threshhold
#  gammainc = G-A increse threshhold

def trgn(x, b, F, imax, eps, kmax, rho, Delta, mured=0.5, muinc=1.8, \
             gammared=0.1, gammainc=0.9, gnverbose=0, cgverbose=0, \
             maxback=0, gnorm0=None):

    try:

        # total CG steps
        ktot = 0
        # total function evals
        jtot = 0
        # total gradient evals
        gtot = 0
            
        # initialize gradient
        
        # res=vcl.Vector(F.getRange())
        # grad=vcl.Vector(F.getDomain())

        # negative residual b - F(x)
        # F.apply(x,res)
        res=F(x)
        jtot+= 1
        res.linComb(1.0,b,-1.0)
        
        # negative gradient
        DFx=F.deriv(x)
        # DFx.applyAdj(res,grad)
        grad = transp(DFx)*res
        gtot+= 1
        gnorm = grad.norm()
        if gnorm0 is None:
            gnorm0 = gnorm

        # storage for step
        s=Vector(F.getDomain())

        # storage for trial update, residual
        #xp = vcl.Vector(F.getDomain())
        #resp = vcl.Vector(F.getRange())

        # current value
        Jc=0.5*res.dot(res)
        
        # GN iteration
        i=0
        if gnverbose>0:
            print('  i      J        |grad J|      Delta')
            print('%3d  %10.4e  %10.4e  %10.4e' % (i, Jc, gnorm, Delta))

        while i<imax and gnorm>eps*gnorm0:

            # compute step
            k = trconjgrad(s, res, DFx, kmax, rho, Delta, cgverbose)
            ktot += k
            print('executed ' + str(k) + ' CG steps, total so far = ' + str(ktot))
            
            actred = 0.0
            predred=0.5*s.dot(grad)

            # trial update, actred, predred
            # xp.copy(x)
            xp = x.dup()
            xp.linComb(1.0,s)
            #F.apply(xp,resp)
            resp = F(xp)
            jtot+= 1
            resp.linComb(1.0,b,-1.0)
            Jp=0.5*resp.dot(resp)
            actred=Jc-Jp
            # remember grad is NEGATIVE grad
            if gnverbose>1:
                print('actred=' + str(actred) + ' predred=' + str(predred))

            # backtracking loop
            j = 0
            while actred < gammared*predred and j < maxback:
                # trust radius reduction
                Delta *= mured
                # scale step
                s.scale(mured)
                # scale predred
                predred *= mured
                # recompute
                xp.copy(x)
                xp.linComb(1.0,s)
                resp = F(xp)
                jtot+= 1
                resp.linComb(1.0,b,-1.0)
                Jp=0.5*resp.dot(resp)
                actred=Jc-Jp
                if gnverbose>1:
                    print('backtrack step ' + str(j) + ':')
                    print('actred=' + str(actred) + ' predred=' + str(predred))
                j += 1

            # if still not there, reduce Delta and re-run CG
            if actred < gammared*predred:
                Delta *= mured

            # update
            else:
                x.copy(xp)
                res.copy(resp)
                # res = resp
                Jc = Jp
                # delete in case significant memory
                # is involved 
                # del DFx
                DFx = F.deriv(x)
                # del grad
                DFx.applyAdj(res,grad)
                gtot+= 1
                # grad = transp(DFx)*res
                gnorm = grad.norm()
                # trust radius increase
                if actred > gammainc*predred:
                    Delta *= muinc
                i=i+1
            if gnverbose > 0:
                print('%3d  %10.4e  %10.4e  %10.4e' % (i, Jc, gnorm, Delta))
                
        print('total function evals     = ' + str(jtot))
        print('total gradient evals     = ' + str(gtot))
        print('total CG steps           = ' + str(ktot))

        return [Delta, gnorm0]
        
    except Exception as ex:
        print(ex)
        raise Exception('called from trgn')
