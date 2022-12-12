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

    p = vcl.Vector(A.getDomain())
    s = vcl.Vector(A.getDomain())
    q = vcl.Vector(A.getRange()) 

    #Initialize:
    #1. $x = 0$
    x.zero()

    #2. $e = b$
    e.copy(b)

    #3. $r = F[m]^Td$
    if not A.applyAdj(b,r):
        print('cg.conjgrad: failed at initial A.applyAdj(b,r)')
        return False

    #4. $p = r$
    p.copy(r)

    #5. $q = F[m]p$
    if not A.applyFwd(p,q):
        print('cg.conjgrad: failed at initial A.applyFwd(p,q)')
        return False

    #6. $s = F[m]^Tq$
    if not A.applyAdj(q,s):
        print('cg.conjgrad: failed at initial A.applyAdj(q,s)')
        return False

    #7. $\gamma_0 = \langle r, r \rangle$
    gamma0 = r.dot(r)

    #8. $\gamma = \gamma_0$
    gamma=gamma0

    #9. $k=0$
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

        #1. $\alpha = gamma / \langle q, q\rangle$
        # write it this way, as quotient of norms, to avoid
        # possible precision issues at sqrt(macheps) level
        alpha = math.sqrt(gamma)/q.norm()
        alpha = alpha*alpha

        #2. $x \leftarrow x+\alpha p$
        x.linComb(alpha,p)

        #3. $e \leftarrow e-\alpha q$
        e.linComb(-alpha,q)
        enorm = e.norm()

        #4. $r \leftarrow r-\alpha s$
        r.linComb(-alpha,s)

        #5. $\delta = \langle r, r \rangle$  
        delta = r.dot(r)
        rnorm = math.sqrt(delta)

        #6. $\beta = \delta / \gamma$
        beta = delta/gamma

        #7. $p \leftarrow r + \beta p$
        p.linComb(1.0,r,beta)

        #8. $\gamma \leftarrow delta$
        gamma=delta

        #9. $q \leftarrow F[m]p$
        if not A.applyFwd(p,q):
            print('cg.conjgrad: failed at A.applyAdj(p,q)')
            return False
        

        #10. $s \leftarrow F[m]^Tq$
        if not A.applyAdj(q,s):
            print('cg.conjgrad: failed at A.applyAdj(q,s)')
            return False

        #11. $k \leftarrow k+1$
        k=k+1

        #12. print $k$ (iteration) $\|e\|$ (residual norm), $\|r\|$ (normal residual norm)
        if verbose > 1:
            print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))

    if verbose > 0:
        print('----------------------------------------------------')
        print('  k       |e|     |e|/|e0|      |r|        |r|/r0|')
        print('%3d  %10.4e  %10.4e  %10.4e  %10.4e' % (k, enorm, enorm/enorm0, rnorm, rnorm/rnorm0))

    return True
        
