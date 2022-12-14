import vcl
import math

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
#  e = residual on return
#  r = normal residual on return
#  verbose = verbosity level
#    0 - no print output
#    1 - print number of its performed, residual and normal residual
#    2 - (or more) print it, res, nres at every it

def conjgrad(x, b, A, kmax, eps, rho, e, r, verbose=0):

    try:

        p = vcl.Vector(A.getDomain())
        s = vcl.Vector(A.getDomain())
        q = vcl.Vector(A.getRange()) 
        
        #Initialize:
        # print('#1. $x = 0$')
        x.zero()
        
        # print('#2. $e = b$')
        e.copy(b)
        
        # print('#3. $r = A^Tb$')
        A.applyAdj(b,r)
    
        # print('#4. $p = r$')
        p.copy(r)

        # print('#5. $q = Ap$')
        A.applyFwd(p,q)

        # print('#6. $s = A^Tq$')
        A.applyAdj(q,s)

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
            # print('k='+str(k))

            # print('#1. $\alpha = gamma / \langle q, q\rangle$')
            # write it this way, as quotient of norms, to avoid
            # possible precision issues at sqrt(macheps) level
            alpha = math.sqrt(gamma)/q.norm()
            alpha = alpha*alpha

            # print('#2. $x \leftarrow x+\alpha p$')
            x.linComb(alpha,p)

            # print('#3. $e \leftarrow e-\alpha q$')
            e.linComb(-alpha,q)
            enorm = e.norm()

            # print('#4. $r \leftarrow r-\alpha s$')
            r.linComb(-alpha,s)

            # print('#5. $\delta = \langle r, r \rangle$')  
            delta = r.dot(r)
            rnorm = math.sqrt(delta)

            # print('#6. $\beta = \delta / \gamma$')
            beta = delta/gamma

            # print('#7. $p \leftarrow r + \beta p$')
            p.linComb(1.0,r,beta)

            # print('#8. $\gamma \leftarrow delta$')
            gamma=delta

            # print('#9. $q \leftarrow Ap$')
            A.applyFwd(p,q)

            # print('#10. $s \leftarrow A^Tq$')
            A.applyAdj(q,s)

            # print('#11. $k \leftarrow k+1$')
            k=k+1

            # print('#12. print $k$ (iteration) $\|e\|$ (residual norm), $\|r\|$ (normal residual norm)')
            if verbose > 1:
                print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))
                
        if verbose > 0:
            print('----------------------------------------------------')
            print('  k       |e|     |e|/|e0|      |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e  %10.4e  %10.4e' % (k, enorm, enorm/enorm0, rnorm, rnorm/rnorm0))

    except Exception as ex:
        print(ex)
        raise Exception("called from cg.conjgrad")
        
