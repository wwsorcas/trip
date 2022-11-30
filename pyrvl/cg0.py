import linalg
import op
import math
import tempfile


# implementation of CG limited to convolution operator
# op.convop, and vectors represented by names of SU files

# make residual and normal residual vectors e and r arguments, so
# that they can be examined after completion. Other work vectors are
# temp files, deleted on return.

# args:
#  w = estimated solution on return
#  d = data
#  g = green's function
#  kmax = max iterations for termination
#  eps = relative residual decrease for termination
#  rho = relative normal residual decrease for termination
#  e = residual on return
#  r = normal residual on return
#  verbose = verbosity level
#    0 - no print output
#    1 - print number of its performed, residual and normal residual
#    2 - (or more) print it, res, nres at every it

def conjgrad0(w, d, g, kmax, eps, rho, e, r, verbose=0):

    pfile = tempfile.NamedTemporaryFile(suffix='.su',delete=False)
    qfile = tempfile.NamedTemporaryFile(suffix='.su',delete=False)
    sfile = tempfile.NamedTemporaryFile(suffix='.su',delete=False)

    p = pfile.name
    q = qfile.name
    s = sfile.name

    # print('p=' + p)
    # print('q=' + q)
    # print('s=' + s)

    # initialize metadata for these files by copies
    # e and r ***should*** be initialized on call, but
    # can avoid presuming this by copies
    linalg.copy(w,p)
    linalg.copy(d,q)
    linalg.copy(w,s)

    linalg.copy(w,r)

    #Initialize:
    #1. $w = 0$
    linalg.scale(vec=w, a=0.0)

    #2. $e = d$
    linalg.copy(d,e)

    #3. $r = F[m]^Td$
    op.convop(g,r,d,adj=1)
    # print('r norm = ' + str(linalg.norm(r)))

    #4. $p = r$
    linalg.copy(r,p)

    #5. $q = F[m]p$
    op.convop(g,p,q,adj=0)
    # print('q norm = ' + str(linalg.norm(q)))

    #6. $s = F[m]^Tq$ 
    op.convop(g,s,q,adj=1)

    #7. $\gamma_0 = \langle r, r \rangle$
    gamma0 = linalg.dot(r,r)

    #8. $\gamma = \gamma_0$
    gamma=gamma0

    #9. $k=0$
    k=0
    enorm0=linalg.norm(e)
    rnorm0=linalg.norm(r)
    enorm=enorm0
    rnorm=rnorm0

    if verbose > 1:
        print('  k       |e|       |r|=')
        print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))
    

    #Repeat while $k<k_{\rm max}$, $\|e\|>\epsilon \|d\|$:
    while k<kmax and enorm>eps*enorm0 and rnorm>rho*rnorm0:

        #1. $\alpha = gamma / \langle q, q\rangle$
        alpha = math.sqrt(gamma)/linalg.norm(q)
        # print('sqrt(alpha)=' + str(alpha))
        alpha = alpha*alpha
        # print('pnorm='+str(linalg.norm(p)))
        # print('wnorm before update ='+str(linalg.norm(w)))
        #2. $w \leftarrow w+\alpha p$
        linalg.lincomb(alpha,p,w)
        # print('wnorm after update ='+str(linalg.norm(w)))
        #3. $e \leftarrow e-\alpha q$
        linalg.lincomb(-alpha,q,e)
        enorm = linalg.norm(e)

        #4. $r \leftarrow r-\alpha s$
        linalg.lincomb(-alpha,s,r)

        #5. $\delta = \langle r, r \rangle$  
        delta = linalg.dot(r,r)
        rnorm = math.sqrt(delta)

        #6. $\beta = \delta / \gamma$
        beta = delta/gamma

        #7. $p \leftarrow r + \beta p$
        linalg.lincomb(1.0,r,p,beta)

        #8. $\gamma \leftarrow delta$
        gamma=delta

        #9. $q \leftarrow F[m]p$
        op.convop(g,p,q,adj=0)

        #10. $s \leftarrow F[m]^Tq$
        op.convop(g,s,q,adj=1)

        #11. $k \leftarrow k+1$
        k=k+1

        #12. print $k$ (iteration) $\|e\|$ (residual norm), $\|r\|$ (normal residual norm)
        if verbose > 1:
            print('%3d  %10.4e  %10.4e' % (k, enorm, rnorm))

    if verbose > 0:
        print('----------------------------------------------------')
        print('  k       |e|     |e|/|e0|      |r|        |r|/r0|')
        print('%3d  %10.4e  %10.4e  %10.4e  %10.4e' % (k, enorm, enorm/enorm0, rnorm, rnorm/rnorm0))
        
