import gc
import os
import vcl
import math
from vcl import Vector
from vcl import transp
from abc import ABC, abstractmethod

# single steps are separated from loops as separate functions in order
# to force garbage-collection.

############## CONJUGATE GRADIENT ITERATION FOR THE NORMAL EQUATION ##############

def cgstep(x, A, k, p, r, e, rnorm, enorm, gamma, sigsq=None, verbose=0):
    ''' 
    Single step of CG iteration for solution of normal equations, with
    optional regularization. See doc for function conjgrad.
    '''
    
#    os.system('echo top; ls /var/tmp')
    try:
        # print('k='+str(k))

        # print('#1. $q = Ap$')
        q = A*p
            
        # print('#2. $s = A^T q + sig^2 p$')
        s = transp(A)*q
        if sigsq is not None:
            s.linComb(sigsq,p)
            
        #os.system('echo q s; ls /var/tmp')
            
        # print('#3. $\alpha = gamma / \langle p, s\rangle$')
        alpha = gamma/p.dot(s)
        #print('|p|=' + str(p.norm()) + ' |q|=' + str(q.norm()) + ' |s|=' + str(s.norm()) + ' alpha=' + str(alpha))

        # print('#4. $x \leftarrow x+\alpha p$')
        x.linComb(alpha,p)

        # print('#5. $e \leftarrow e-\alpha q$')
        #print('|e| before =' + str(e.norm()))
        e.linComb(-alpha,q)
        enorm = e.norm()
        #print('|e| after  =' + str(enorm))
        
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
                
def conjgrad(x, b, A, kmax, eps, rho, sig=None, verbose=0, e=None, r=None):
    '''
    conjugate gradient iteration for the normal equation, with optional
    Tihonov regularization. Approximate minimization of

    |Ax-b|^2 + sig^2 |x|^2

    The second term is optional: to avoid it without extra arithmetic, set
    sig = None (not sig = 0).

    Aim is to reduce residual norm |Ax-b| or normal residual norm |A^T(Ax-b)|
    below a given fraction of initial value.

    Note that the "residual" in this discussion means b-Ax, not the augmented
    residual (b-Ax, sig x)

    Iteration terminates if
        - max number of iterations reached
        - residual norm reduced below fraction of initial value
        - normal residual norm reduced below fraction if initial value

    Parameters:
    x (vcl.Vector):           solution estimate, initially = 0
    b (vcl.Vector):           right-hand side
    A (vcl.LinearOperator):   linear op defining problem
    kmax (int):               iteration limit
    eps (float):              residual reduction
    rho (float):              normal residual reduction
    sig (float):              optional regularization weight (float or None)
    verbose (int):            verbosity flag, 0 = no output, 
                                              1 = end of iteration summary, 
                                              2 = step-by-step
    e (vcl.Vector):           optional return, residual
    r (vcl.Vector):           optional return, normal residual

    NOTE: e here is b-Ax, not Ax-b.
    '''
    
    try:

        #Initialize:
        # print('#1. $x = 0$')
        x.scale(0.0)
        
        # print('#2. $e = b$')
        if e is None:
            e = vcl.Vector(b.space)
        e.copy(b)
        
        # print('#3. $r = A^Tb$')
        #print('vcalg.conjgrad')
        #print('|b|=' + str(b.norm()))
        p=transp(A)*b
        #print('|p|=' + str(p.norm()))
        #p.myNameIs()
        
        # print('#4. $p = r$')
        if r is None:
            r = vcl.Vector(p.space)
        r.copy(p)
        #print('|r|=' + str(r.norm()))
        
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
    
        # prepare square of reg weight if used
        sigsq = None
        if sig is not None:
            sigsq = sig*sig
            
        #Repeat while $k<k_{\rm max}$, $\|e\|>\epsilon \|d\|$:
        while k<kmax and enorm>eps*enorm0 and rnorm>rho*rnorm0:

            [k, enorm, rnorm, gamma] = cgstep(x, A, k, p, r, e, rnorm, enorm, gamma, sigsq, verbose)
            
        if verbose > 0:
            print('----------------------------------------------------')
            print('  k       |e|     |e|/|e0|      |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e  %10.4e  %10.4e' % (k, enorm, enorm/enorm0, rnorm, rnorm/rnorm0))

    except Exception as ex:
        print(ex)
        raise Exception("called from cg.conjgrad")

class cgne(vcl.LSSolver):
    '''
    object interface for CGNE algorithm, for use in VPM and other
    applications requiring solution of reg. least squares problem:
    min_x |Ax-b|^2 + sig^2|x|^2.
    in terms of this description, the args are
    op = A
    rhs = b
    Constructor stores CG parameters. Solve method sanity-checks and
    calls CG algorithm, then returns solution x and residual b-Ax as
    a vector list. Regularization weight sig should be set to None if 
    (Tihonov) regularization not desired, to avoid useless flops.
    '''
    def __init__(self, kmax, eps, rho, sig=None, verbose=0):
        self.kmax = kmax
        self.eps = eps
        self.rho = rho
        self.sig = sig
        self.verbose = verbose

    def solve(self, op, rhs):
        try:
            if not isinstance(op, vcl.LinearOperator):
                raise Exception('first arg not LinearOperator')
            if not isinstance(rhs, vcl.Vector):
                raise Exception('second arg not Vector')
            if rhs.space != op.getRange():
                raise Exception('second arg not in range of first arg')
            x = vcl.Vector(op.getDomain())
            e = vcl.Vector(op.getRange())
            conjgrad(x, rhs, op, self.kmax, self.eps, self.rho, sig=self.sig,
                        verbose=self.verbose, e=e, r=None)
            return [x, e]
        except Exception as ex:
            print(ex)
            raise Exception('called from vcalg.cgnefcn.solve')

    def myNameIs(self):
        print('Conjugate Gradients for the Normal Equations')
        print('  max iterations = ' + str(self.kmax))
        print('  residual tol   = ' + str(self.eps))
        print('  normal res sol = ' + str(self.rho))
        print('  verbosity flag = ' + str(self.verbose))
        

########### TRUST-RADIUS MODIFIED CONJUGATE GRADIENT ITERATION FOR NORMAL EQUATION ##########

def trcgstep(x, A, k, p, r, gamma,
                 Delta=None, verbose=0):
    try:

        # print('k='+str(k))

        # print('#1. $q = Ap$')
        # A.applyFwd(p,q)
        q = A*p
    
        # print('#2. $s = A^Tq$')
        # A.applyAdj(q,s)
        s = transp(A)*q
        
        # print('#3. $\alpha = gamma / \langle p, s\rangle$')
        alpha = gamma/p.dot(s)
        #print('alpha = ' + str(alpha))
        
        # print('#4. $x \leftarrow x+\alpha p$')
        x.linComb(alpha,p)

        cgend = False
        xnorm = x.norm()
        #print('xnorm = ' + str(xnorm))
        
        if Delta is not None:
            if xnorm > Delta:
                x.scale(Delta/xnorm)
                cgend = True
                if verbose > 0:
                    print('trust radius exceeded - return scaled ' + \
                            'soln of length ' + str(Delta))
                        
        if not cgend:
                
            # print('#6. $r \leftarrow r-\alpha s$')
            r.linComb(-alpha,s)
                
            # print('#7. $\delta = \langle r, r \rangle$')
            gammaold = gamma
            gamma = r.dot(r)
            #print('gamma = ' + str(gamma))
            
            # print('#8. $\beta = \delta / \gamma$')
            beta = gamma/gammaold
            #print('beta = ' + str(beta))
                
            # print('#9. $p \leftarrow r + \beta p$')
            p.linComb(1.0,r,beta)
                
            # print('#11. $k \leftarrow k+1$')
            k=k+1
            
            # print('#12. print $k$ (iteration) $\|e\|$ (residual norm),
            # $\|r\|$ (normal residual norm)')
            #if verbose > 1:
            #    print('%3d  %10.4e' % (k, rnorm))
                
        del q
        del s
        # force garbage collection
        n=gc.collect()
        
        return [x, k, p, r, gamma, cgend]
    
    except Exception as ex:
        print(ex)
        raise Exception('called from trcgstep')
    
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
        p = vcl.Vector(r.space)
        p.copy(r)

        # print('#5. $\gamma_0 = \langle r, r \rangle$')
        gamma0 = r.dot(r)
        
        # print('#6. $\gamma = \gamma_0$')
        gamma=gamma0
        
        # print('#7. $k=0$')
        k=0
        rnorm0=r.norm()
        rnorm=rnorm0

        # flag to indicate trust region truncation, indicating
        # return to gn loop
        cgend = False

        if verbose > 0:
            print('-------------------->>> CG >>>------------------------')
            print('  k       |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e' % (k, rnorm, rnorm/rnorm0))            

        #Repeat while $k<k_{\rm max}$, rnorm > rho*rnorm0, |x|<Delta
        while k<kmax and rnorm>rho*rnorm0 and not cgend:
            #print('gamma before: ' + str(gamma))
            [x, k, p, r, gamma, cgend] = trcgstep(x, A, k, p, r, gamma, Delta, verbose)
            #print('gamma after:  ' + str(gamma))

            rnorm = math.sqrt(gamma)

            if verbose > 0 and not cgend:
                print('%3d  %10.4e  %10.4e' % (k, rnorm, rnorm/rnorm0))

        if verbose > 0:            
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
            print('\n GN Iteration ' + str(i))
            print('  i      J        |grad J|      Delta')
            print('%3d  %10.4e  %10.4e  %10.4e' % (i, Jc, gnorm, Delta))

        while i<imax and gnorm>eps*gnorm0:

            # compute step
            k = trconjgrad(s, res, DFx, kmax, rho, Delta, cgverbose)
            ktot += k
#            print('executed ' + str(k) + ' CG steps, total so far = ' + str(ktot))
            
            actred = 0.0
            predred=0.5*s.dot(grad)

            # trial update, actred, predred
            # xp.copy(x)
            xp = vcl.Vector(x.space)
            xp.copy(x)
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

            # if still not there, bag it
            if actred < gammared*predred:
                print('******** trgn: backtrack loop failed *******')
                return 

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

# single steps are separated from loops as separate functions in order
# to force garbage-collection.

###################### BASIC CONJUGATE GRADIENT ITERATION ######################

# Basic algorithm, as opposed to CGNE implemented above

# args:
#  x = estimated solution on return
#  b = data
#  A = SPD operator
#  kmax = max iterations for termination
#  rho = relative residual decrease for termination
#  r = residual on return (optional)
#  verbose = verbosity level
#    0 - no print output
#    1 - print number of its performed, residual norm
#    2 - (or more) print it, res at every it

def bcgstep(x, A, k, p, r, gamma, Delta=None, verbose=0):
    '''
    Iteration step, basic CG algorithm for SPD systemts.

    Parameters:
    x (vcl.Vector):           solution estimate to be updated
    A (vcl.LinearOperator):   SPD operator
    k (int):                  iteration number
    p (vcl.Vector):           update direction
    r (vcl.Vector):           residual on return (optional)
    gamma (float):            residual norm squared
    Delta (float):            optional trust radius constraint
    verbose (int):            verbosity level
                              0 - no print output
                              1 - print number of its performed, residual norm
                              2 - (or more) print it number, res at every it

    Return: [x, k, p, r, gamma, cgend]
    all updated parameters except
    cgend (Boolean):          True if trust region transgressed, else False

    '''

    try:
            
        # print('#2. $s = A^Tq$')
        s = A*p
            
        # print('#3. $\alpha = gamma / \langle p, s\rangle$')
        alpha = gamma/p.dot(s)
        #print('alpha = ' + str(alpha))
        
        # print('#4. $x \leftarrow x+\alpha p$')
        x.linComb(alpha,p)

        cgend = False
        xnorm = x.norm()

        if Delta is not None:
            if xnorm > Delta:
                x.scale(Delta/xnorm)
                cgend = True
                if verbose > 0:
                    print('trust radius exceeded - return scaled ' + \
                              'soln of length ' + str(Delta))

        if not cgend:
            
            # print('#6. $r \leftarrow r-\alpha s$')
            r.linComb(-alpha,s)
        
            # print('#7. $\delta = \langle r, r \rangle$')
            gammaold = gamma
            gamma = r.dot(r)
            #print('gamma = ' + str(gamma))

            # print('#8. $\beta = \delta / \gamma$')
            beta = gamma/gammaold
            #print('beta = ' + str(beta))            
        
            # print('#9. $p \leftarrow r + \beta p$')
            p.linComb(1.0,r,beta)
        
            # print('#11. $k \leftarrow k+1$')
            k=k+1
        
            # print('#12. print $k$ (iteration) $\|e\|$ (residual norm),
            # $\|r\|$ (normal residual norm)')
            #if verbose > 1:
            #    print('%3d  %10.4e' % (k, math.sqrt(gamma)))

        # without the following code, the vector destructor is not called
        # until the program exits
        # reduce ref count for q, s to 1, which is unreachable - the
        # reference created by invoking the constructor 
        del s
        # force garbage collection
        n=gc.collect()

        return [x, k, p, r, gamma, cgend]

    except Exception as ex:
        print(ex)
        raise Exception('called from bcgstep')
                
def bcg(x, b, A, kmax, rho, verbose=0, r=None, Delta=None):
    '''
    Standard CG algorithm for solution of SPD linear systems Ax=b. 
    Initial estimate is zero vector. Residual vector is optional argument.
    Optionally implements trust radius truncation. Note that if truncation
    is active, then residual has not been calculated for truncated solution
    on return - that task is left to the calling unit.
    
    Parameters:
    x (vcl.Vector): solution estimate
    b (vcl.Vector): right-hand side vector
    A (vcl.LinearOperator): SPD linear map
    kmax (int): iteration limit
    rho (float): residual tolerance
    verbose (int): verbosity level - 0 = no output, 1 = summary at 
        end if iteration, > 1 = it count and residual at every iteration.
    r (vcl.Vector): residual vector, optional
    Delta (float): trust radius, optional

    Return value: k (int) = number of iterations

    '''
    try:

        #Initialize:
        # print('#1. $x = 0$')
        x.scale(0.0)

        # initial search direction = initial negative residual = b
        p = vcl.Vector(x.space)
        p.copy(b)
        
        # print('#4. $p = r$')
        if r is None:
            r = vcl.Vector(x.space)
        r.copy(p)
        #print('|r|=' + str(r.norm()))
        
        # print('#7. $\gamma_0 = \langle r, r \rangle$')
        gamma0 = r.dot(r)
        
        # print('#8. $\gamma = \gamma_0$')
        gamma=gamma0
        
        # print('#9. $k=0$')
        k=0
        rnorm = math.sqrt(gamma)
        rnorm0=rnorm
    
        cgend = False
        
        if verbose > 0:
            print('-------------------->>> CG >>>------------------------')
            print('  k       |r|        |r|/r0|')
            print('%3d  %10.4e  %10.4e' % (k, rnorm, rnorm/rnorm0))

        #Repeat while $k<k_{\rm max}$, $\|e\|>\epsilon \|d\|$:
        while k<kmax and rnorm>rho*rnorm0 and not cgend:
            #print('gamma before: ' + str(gamma))
            [x, k, p, r, gamma, cgend] = \
                bcgstep(x, A, k, p, r, gamma, Delta, verbose)
            #print('gamma after:  ' + str(gamma))

            rnorm = math.sqrt(gamma)            

            if verbose > 0 and not cgend:
                print('%3d  %10.4e  %10.4e' % (k, rnorm, rnorm/rnorm0))

        if verbose > 0:            
            print('--------------------<<< CG <<<------------------------')

    except Exception as ex:
        print(ex)
        raise Exception("called from bcg")

    else:
        return k

def trnewt(x, J, newtmax, newteps, cgmax, cgeps, Delta, mured=0.5, muinc=1.8, \
             gammared=0.1, gammainc=0.9, nverbose=0, cgverbose=0, \
             maxreds=0, gnorm0=None, jetargs=None):
    '''
    Generic Newton or Newton-esque algorithm for minimization of a scalar
    function. The objective function is represented by its ScalarJet class
    J. Since the class definition (which has no persistent state controlled 
    at runtime) is passed, provision is made for additional keyword-indexed 
    arguments to the ScalarJet subclass constructor.

    Parameters:
    x (vcl.Vector): solution estimate, Vector in domain of J
    J (class): ScalarJet subclass
    newtmax (int): max Newton iterates
    newteps (float): gradient reduction - GN stopping criterion
    cgmax (int): max CG iterates per GN iteration
    cgeps (float): normal residual reduction - CG stopping criterion
    Delta (float): trust radius
    mured (float): Delta reduction factor
    muinc (float): Delta increase factor
    gammared (float): G-A reduction threshhold
    gammainc (float): G-A increase threshhold
    maxreds (int): max trust radius reductions
    nverbose (int): Newton verbosity flag
    cgverbose (int): CG verbosity flag
    gnorm0 (float): base gradient norm, defaults to initial, 
        for relative reduction test. Assign saved value to 
        re-start iteration
    jetargs (dict): optional keywrod args to jet constructor
    '''

    try:

        # total CG steps
        ktot = 0
        # total function evals
        jtot = 0
        # total gradient evals
        gtot = 0
            
        # Initialize scalar function jet
        Jx = J(x,**jetargs)
        
        # negative gradient
        grad = Jx.gradient()
        grad.scale(-1.0)
        gtot+= 1

        gnorm = grad.norm()
        if gnorm0 is None:
            gnorm0 = gnorm

        # storage for step
        s=Vector(x.space)
        # storage for trial update
        xp = vcl.Vector(x.space)
        
        # current value
        Jc=Jx.value()
        jtot+= 1
        
        # Newton iteration
        i=0
        if nverbose>0:
            print('\nNewton Iteration ' + str(i))
            print('  i      J        |grad J|      Delta')
            print('%3d  %10.4e  %10.4e  %10.4e' % (i, Jc, gnorm, Delta))

        while i<newtmax and gnorm>newteps*gnorm0:

            # compute step
            k = bcg(s, grad, Jx.Hessian(), cgmax, cgeps,
                        verbose=cgverbose, r=None, Delta=Delta)
            ktot += k
            #print('executed ' + str(k) + ' CG steps, total so far = ' +
            #          str(ktot))
            
            actred = 0.0
            # remember grad is NEGATIVE grad
            predred=0.5*s.dot(grad)

            # trial update, actred, predred

            xp.copy(x)
            xp.linComb(1.0,s)
            Jxp = J(xp,**jetargs)
            Jp = Jxp.value()
            jtot+= 1
            actred=Jc-Jp
            if nverbose>1:
                print('actred=' + str(actred) + ' predred=' + str(predred))

            # backtracking loop
            j = 0
            while actred < gammared*predred and j < maxreds:
                # trust radius reduction
                Delta *= mured
                # scale step
                s.scale(mured)
                # scale predred
                predred *= mured
                # recompute
                xp.copy(x)
                xp.linComb(1.0,s)
                Jxp = J(xp,**jetargs)
                Jp = Jxp.value()
                jtot+= 1
                actred=Jc-Jp
                if nverbose>1:
                    print('backtrack step ' + str(j) + ':')
                    print('actred=' + str(actred) + ' predred=' + str(predred))
                j += 1

            # if still not there, reduce Delta and re-run CG
            if actred < gammared*predred:
                print('******** trnewt: backtrack loop failed *******')
                return

            # update
            else:
                x.copy(xp)
                # res = resp
#                Jx = J(x,**jetargs)
#               Jx = Jxp
#                Jc = Jx.value()
#                jtot += 1
                Jx = Jxp
                Jc = Jp
                grad = Jx.gradient()
                grad.scale(-1.0)
                gtot+= 1
                gnorm = grad.norm()
                # trust radius increase
                if actred > gammainc*predred:
                    Delta *= muinc
                i=i+1
                if nverbose > 0:
#                print('\nNewton Iteration ' + str(i))
                    print('%3d  %10.4e  %10.4e  %10.4e' % (i, Jc, gnorm, Delta))
                
        print('total function evals     = ' + str(jtot))
        print('total gradient evals     = ' + str(gtot))
        print('total CG steps           = ' + str(ktot))

        return [Delta, gnorm]
        
    except Exception as ex:
        print(ex)
        raise Exception('called from trcgnewt')

class SearchDir(ABC):

     # accepts Jet instance,
     # returns search direction
     @abstractmethod
     def Update(self, Jx):
         pass

class SDwgrad(SearchDir):
        
    def __init__(self, Winv=None):
        '''
        Weghted gradient ascent direction.

        Parameters:
        Winv (vcl.LinearOperator): SPD operator (inverse of weight defining inner product)
        
        Update method Returns:
        Winv*g:                    weighted gradient = gradient in weighted norm
        '''
        try:
            if not Winv is None:
                if not isinstance(Winv, vcl.LinearOperator):
                    raise Exception('input Winv not LinOp')
                if not Winv.getDomain() == Winv.getRange():
                    raise Exception('Winv not square')
            self.Winv = Winv
        except Exception as ex:
            print(ex)
            raise Exception('called from vcalg.SDwgrad constructor')

    def Update(self, Jx):
        try:
            if not isinstance(Jx, vcl.ScalarJet):
                raise Exception('input arg not vcl.ScalarJet')
            if self.Winv is None:
                return Jx.gradient()
            else:
                if Jx.gradient().space != self.Winv.getDomain():
                    raise Exception('gradient not in domain of Winv')
                return self.Winv*Jx.gradient()
        except Exception as ex:
            print(ex)
            raise Exception('called from vcalg.SDwgrad.Update')

class LineSearch(ABC):
    '''
    Abstract line search class. Stores generically required
    items

    Constructor Parameters:
    x (vcl.Vector):       current solution estimate
    val (float):          current objective value
    grad (vcl.Vector):    current gradient
    dir (vcl.Vector):     search direction
    step (float):         initial step
    J (vcl.ScalarJet):    objective jet (class)
    jetargs (dict):       kwargs for jet
    verbose (int)         verbosity flag (silent = 0)
    fout (iofile):        output unit, or none for stdout
#    lsargs:               additional arguments (keyword)
    '''

    def __init__(self, x, val, grad, dir, step, J, jetargs, verbose, fout):
        self.x = x
        self.val = val
        self.grad = grad
        self.dir = dir
        self.step = step
        self.J = J
        self.jetargs = jetargs
        self.verbose = verbose
        self.fout = fout
            
    # returns updated jet instance and step
#    @abstractmethod
    def search(self):
        raise Exception('Error: somehow called base class method LineSearch.search')
        
class btls(LineSearch):

    '''
    Simple acktracking line search. 

    Base class parameters:
    x (vcl.Vector):       current solution estimate
    val (float):          current objective value
    grad (vcl.Vector):    current gradient
    dir (vcl.Vector):     search direction
    step (float):         initial step
    J (vcl.ScalarJet):    objective jet (class)
    jetargs (dict):       keyword args for jet constructor
    lsverbose (int):      verbosity flag
    fout (file):          output file - always appended!!!

    Parameters:
    lsmax (int):          max steps
    mured (float):        try shorter step if actred < mured*predred
    muinc (float):        try longer step if actred > muinc*predred
    gammared (float):     step reduction factor
    gammainc (float):     step increase factor

    return from search():
    Jx (vcl.ScalarJet):   evaluated jet at final step, or None for failure
    step (float):         updated step
    '''
    
    def __init__(self, x, val, grad, dir, step, J, jetargs=None, lsverbose=0, fout=None,
             lsmax=1, mured=0.5, muinc=1.8, gammared=0.1, gammainc=0.9):

#        self.x = x
#        self.val = val
#        self.grad = grad
#        self.dir = dir
#        self.step = step
#        self.J = J
#        self.jetargs = jetargs

# note explicit invokation of superclass constructor with self pointer
# quite different from C++!
        if fout is None:
            raise Exception('output file object not provided (fout)')
        LineSearch.__init__(self,x, val, grad, dir, step, J, jetargs, lsverbose, fout)

        self.lsmax = lsmax
        self.mured = mured
        self.muinc = muinc
        self.gammared = gammared
        self.gammainc = gammainc

    def search(self):
        try:
            k = 0
            takestep = 0
        
            while k < self.lsmax and takestep==0:
                if self.verbose > 0:
                    print('\n    Line Search Step ' + str(k), file=self.fout)
                    print('        update step, jet', file=self.fout)
                xtest = vcl.Vector(self.x.space)
                xtest.copy(self.x)
                xtest.linComb(-self.step,self.dir)
                Jx = self.J(xtest, **self.jetargs)
                vtest = Jx.value()
                if self.verbose > 0:
                    print('        step = %10.4e val = %10.4e' % (self.step, vtest), file=self.fout)
                # G-A test
                actred = self.val-vtest
                predred = self.step*self.dir.dot(self.grad)
                if self.verbose > 0:
                    print('        actred = %10.4e predred = %10.4e' % (actred,predred), file=self.fout)
                if self.gammainc*predred < actred and k < self.lsmax-1:
                    if self.verbose > 0:
                        print('        try longer step', file=self.fout)
                    takestep = 0
                    self.step *= self.muinc
                    k += 1
                elif self.gammared*predred < actred:
                    if self.verbose > 0:
                        print('        in G-A range', file=self.fout)
                    takestep = 1
                elif k == self.lsmax-1:
                    if self.verbose > 0:
                        print('        last allowed step - save estimate', file=self.fout)
                    takestep = 1
                else:
                    if self.verbose > 0:
                        print('        try shorter step', file=self.fout)
                    takestep = 0
                    self.step *= self.mured
                    k += 1

            if takestep == 1:
                return [Jx, self.step]
            else:
                return [None, self.step]

        except Exception as ex:
            print(ex)
            raise Exception('called from vcalg.btls')

import sys

# lsargs - see above
def lsopt(x, J, SD=None, LS=None, descmax=0, desceps=0.01, descverbose=0, descout=None,
              lsverbose=0, lsargs=None, jetargs=None, ddargs=None, archargs=None ):
    '''
    Line seach optimization algorithm. Two principal components:
    - SD: SearchDir class: Update methodcomputes search direction from 
    current gradient and other parameters (ddargs)
    - LS: line search function, uses direction from SD and other parameters (lsargs)
    
    Parameters:
    x (vcl.Vector):     initial estimate on call
    J (vcl.ScalarJet):  class: jet of objective function 
    SD (SearchDir):     class: returns search direction based on gradient and other params (ddargs)
    LS (LineSearch):    class: implements line search 
    descmax (int):      max number of descent steps
    desceps (float):    terminates if gradient falls below this proportion of initial
    descverbose (int):  verbosity flag (default = 0)
    descout (string):   output file (default = None)

    keyword dictionaries:
    lsargs:             line search params
    jetargs:            jet params
    ddargs:             search direction params
    archargs:           archive params

    Returns:
    Jx (vcl.ScalarJet): final instance of jet class - final state of search
    '''
    
    try:
        # total function evals
        jtot = 0
        # total gradient evals
        gtot = 0
        #
        initprop = 0.125

        if descout is not None:
            # check for legit
            if not isinstance(descout, str):
                raise Exception('descout arg not string')
#            print('vcalg.lsopt: decout = ' + descout)
            fout   = open(descout, 'w')
        else:
            fout   = sys.stdout
            
        if descverbose != 0:
            print('\nLine Search Optimization', file=fout)
            print('\ninitialize jet', file=fout)

        Jx = J(x, **jetargs)

        if descverbose != 0:
            print('compute initial descent direction', file=fout)
        if SD is None:
            raise Exception('Descent Direction class not provided')
        SDinstance = SD(**ddargs)
        ddir = SDinstance.Update(Jx)
        gtot += 1

        if descverbose != 0:
            print('compute initial ascent rate', file=fout)
        dr = Jx.gradient().dot(ddir)

        if descverbose != 0:
            print('sanity check for sufficient ascent', file=fout)
        if 1.0 + dr <= 1.0:
            raise Exception('ascent rate = ' + str(dr) + ' insufficient or negative')

        # initial step calculation
        # assumes target value = 0
        print('compute initial step', file=fout)
        step = initprop*Jx.value()/dr
        jtot += 1
        
        # descent iteration
        i = 0
        # archive
        for k in archargs.keys():
            Jx.archive(k,str(i),archargs[k])
        
        if descverbose !=0:
            print('initial value = %10.4e'  % (Jx.value()), file=fout)
            print('initial step  = %10.4e' % (step), file=fout)
            print('initial ascent rate = %10.4e' % (dr), file=fout)
            fout.flush()
            
        more = True
        drinit = dr
        while i < descmax and dr > desceps*drinit and more:
            if descverbose > 0:
                print('\nIteration ' + str(i), file=fout)
            LSinstance = LS(Jx.point(), Jx.value(), Jx.gradient(), ddir, step,
                                J, jetargs=jetargs, lsverbose=lsverbose,
                                fout=fout, **lsargs)
            [Jxtest, steptest] = LSinstance.search()
            # test for successful step
            if Jxtest is not None:
                if descverbose > 0:
                    print('\naccept new jet and step', file=fout)
                Jx = Jxtest
                step = steptest
                if descverbose > 0:
                    print('value = %10.4e step = %10.4e' % (Jx.value(),step), file=fout)
                if descverbose > 0:
                    print('update search direction', file=fout)
                ddir = SDinstance.Update(Jx)
                # update ascent rate
                dr = Jx.gradient().dot(ddir)
                # update counter
                i += 1
                # archive
                for k in archargs.keys():
                    Jx.archive(k,str(i),archargs[k])
                fout.flush()
            else:
                # bail out
                more = False
                if descverbose > 0:
                    print('line search failed, exit lsopt', file=fout)
                    
        if more == True:
            if descverbose > 0:
                if dr <= drinit:
                    print('\nachieved prescribed ascent rate reduction:', file=fout)
                if i==descmax:
                    print('\nreached iteration limit', file=fout)
                print('final value = %10.4e'  % (Jx.value()), file=fout)
                print('final step  = %10.4e' % (step), file=fout)
                print('finial ascent rate = %10.4e' % (dr), file=fout)
                print('exit lsopt', file=fout)

        fout.close()
        return Jx
    
    except Exception as ex:
        print(ex)
        raise Exception('called from vcalg.lsopt')

def alphaupdate(alpha=0.0, e=0.0, ap=1.0, eplus=1.1, eminus=0.9):
    '''
    alpha update rule for implementing discrepancy principle, following
    Fu & Symes Geophysics 2017 and Symes et al. IP 2023. Assumes that 
    (in notation of Symes et al.) e = 0.5*\|F(x)-d\|^2/\|d\|^2 is 
    relative data fit error, and ap = 0.5*alpha^2*\|Ax\|^2/\|d\|^2 is 
    relative regularization error (both squared, of course). eplus is 
    upper permitted half relative fit error squared, eminus is lower. 
    alpha is current value of alpha.

    return value is [alpha, update], with alpha = updated value, update
    = True if update performed, False is alpha unchanged

    alpha <- sqrt(alpha^2 + (eplus - e)/(2*p))

    (Symes et al. equation 17)

    where p = ap/alpha. Rearranging,

    alpha <- alpha*sqrt(1 + (eplus - e)/(2*ap))
    '''

    if e < eminus or e > eplus:
        try:
            recip = 1.0/(2.0*ap)
        except ZeroDivisionError as e:
            raise Exception('zerodivide in alpha update formula')
    
        return alpha*sqrt(1 + (eplus - e)*recip)
    else:
        return [alpha, False]
        

    





        



