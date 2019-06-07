// cgnealg_mod.hh
// Author: Mario J. Bencomo
// last modified 06/18/18

// Basically a modification of cgnealg.hh


#ifndef __CGNEALG_MOD_HH
#define __CGNEALG_MOD_HH

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"
#ifdef IWAVE_USE_MPI
#include "mpi.h"
#endif


using namespace RVLAlg;

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;    

  // forward declaration
  template<typename Scalar>
  class CGNEAlg_mod;

  /** Single step of conjugate gradient iteration for the normal
      equations.

      Modifications were made to monitor xnorm=norm(x) and xerrs=norm(x-x0) 
      for a given x0.
  */
  template<typename Scalar>
  class CGNEStep_mod : public Algorithm {

    friend class CGNEAlg_mod<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    Vector<Scalar> const & x0;
    atype & rnorm;
    atype & nrnorm;
    atype & xnorm;
    atype & xerrs;
    
  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> g;    // gradient
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction
    atype gamma;         // ||g||^2
    atype x0norm;

  public:

    /**
     * Constructor
     */
    CGNEStep_mod(LinearOp<Scalar> const & _A,
		 Vector<Scalar> & _x,
		 Vector<Scalar> const & _b,
		 Vector<Scalar> const & _x0,
		 atype & _rnorm, 
		 atype & _nrnorm,
		 atype & _xnorm,
		 atype & _xerrs)
      : A(_A), x(_x), b(_b), x0(_x0),
	rnorm(_rnorm), nrnorm(_nrnorm), 
	xnorm(_xnorm), xerrs(_xerrs), 
	r(A.getRange()), g(A.getDomain()),
	q(A.getRange()), p(A.getDomain()) { 

      // NOTE: initial x assumed to be zero vector
      r.copy(b);
      rnorm=r.norm();
      A.applyAdjOp(r,g);
      p.copy(g);
      nrnorm=g.norm();
      gamma = nrnorm*nrnorm;
      
      x0norm = x0.norm();      
      xnorm = x.norm();
      Vector<Scalar> tmp(x);
      tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
      atype tmp2=tmp.norm();
      if (ProtectedDivision<atype>(tmp2,x0norm,xerrs))
	x0norm=ScalarFieldTraits<Scalar>::One();
    }
      
    /**
     * Run a single step of CGNE
     */
    void run() {
      try {

	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e << "Error: CGNEStep_mod::run() from ProtectedDivision: alpha = gamma/qtq\n"
	    << "  gamma = "<< gamma <<"\n"
	    << "  qtq = "<< qtq <<"\n";
	  throw e;
	}
	Scalar alpha=absalpha;

	x.linComb(alpha,p);
	r.linComb(-alpha,q);
	A.applyAdjOp(r,g);

	atype newgamma = g.normsq();
	atype absbeta;
	if (ProtectedDivision<atype>(newgamma,gamma,absbeta)) {
	  RVLException e;
	  e << "Error: CGNEStep_mod::run() from ProtectedDivision: beta=newgamma/gamma\n"
	    << "  newgamma = "<< newgamma <<"\n"
	    << "  gamma = "<< gamma <<"\n";
	  throw e;
	}
	Scalar beta = absbeta;

	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);

	xnorm=x.norm();
	Vector<Scalar> tmp(x);
	tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
	xerrs=tmp.norm()/x0norm;	

      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEStep_mod::run()\n";
	throw e;
      }
     
    }

    ~CGNEStep_mod() {}
  };
  
  /** Preconditioned conjugate gradient iteration for the normal
      equations.

      Preconditioner is of the form M=(L^T L)^-1 for a given invertible 
      linear operator L.
      Modifications were made to monitor xnorm=norm(x), xerrs=norm(x-x0),
      Lxnorm=norm(Lx), and Lxerrs=norm(L(x-x0)) for given x0 and L.
  */
  template<typename Scalar>
  class PCGNEStep_mod : public Algorithm {

    friend class CGNEAlg_mod<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    LinearOpWithInverse<Scalar> const & L;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    Vector<Scalar> const & x0;
    atype & rnorm;
    atype & nrnorm;
    atype & xnorm;
    atype & Lxnorm;
    atype & xerrs;
    atype & Lxerrs;
    ostream & str;

  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> ATr;  // normal residual
    Vector<Scalar> s;    // s = LinvT * AT * r
    Vector<Scalar> g;    // gradient = preconditioned normal residual
                         // g = Linv * LinvT * AT * r
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction    
    Vector<Scalar> Lx;  
    atype gamma;         
    atype x0norm;
    atype Lx0norm;

  public:
  
    /**
     * Constructor
     */
    PCGNEStep_mod(LinearOp<Scalar> const & _A,
		  LinearOpWithInverse<Scalar> const & _L,
		  Vector<Scalar> & _x,
		  Vector<Scalar> const & _b,
		  Vector<Scalar> const & _x0,
		  atype & _rnorm, 
		  atype & _nrnorm,
		  atype & _xnorm,
		  atype & _Lxnorm,
		  atype & _xerrs,
		  atype & _Lxerrs,
		  ostream & _str)
      : A(_A), L(_L), x(_x), b(_b), x0(_x0), 
	rnorm(_rnorm), nrnorm(_nrnorm), 
	xnorm(_xnorm), Lxnorm(_Lxnorm),
	xerrs(_xerrs), Lxerrs(_Lxerrs),
	str(_str),
	r(A.getRange()),   ATr(A.getDomain()), 
	g(A.getDomain()),  s(A.getDomain()), 
	p(A.getDomain()),  q(A.getRange()),
	Lx(A.getDomain()) { 

      // sanity test, L should be a "square" operator compatible with A
      if ((L.getDomain() != A.getDomain()) ||
	  (L.getRange()  != A.getDomain())) {
	RVLException e;
	e << "Error PCGNEStep_mod constructor:\n"
	  << "Tikhonov operator L is not compatible with system operator A.\n";
	e << "Printing Tikhonov operator L:\n";
	L.write(e);
	e << "Printing system operator A:\n";
	A.write(e);
	throw e;
      }

      // NOTE: initial x assumed to be zero vector
      r.copy(b);
      rnorm=r.norm();
      A.applyAdjOp(r,ATr);
      L.applyInvAdjOp(ATr,s);
      L.applyInvOp(s,g);
      p.copy(g);
      gamma = s.normsq();
      nrnorm=sqrt(gamma);

      x0norm = x0.norm();
      L.applyOp(x0,Lx);
      Lx0norm = Lx.norm();

      xnorm = x.norm();
      L.applyOp(x,Lx);
      Lxnorm = Lx.norm();

      Vector<Scalar> tmp(x);
      tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
      atype tmp2 = tmp.norm();
      if (ProtectedDivision<atype>(tmp2,x0norm,xerrs))
	x0norm=ScalarFieldTraits<Scalar>::One();

      L.applyOp(tmp,Lx);
      tmp2 = Lx.norm();
      if (ProtectedDivision<atype>(tmp2,Lx0norm,Lxerrs))
	Lx0norm=ScalarFieldTraits<Scalar>::One();
    }
      
    /**
     *  Run a single step of preconditioned CGNE
     */
    void run() {
      try {

	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e << "Error: PCGNEStep_mod::run() from ProtectedDivision: alpha = gamma/qtq\n"
	    << "  gamma = "<< gamma <<"\n"
	    << "  qtq = "<< qtq <<"\n";
	  throw e;
	}
	Scalar alpha=absalpha;

	x.linComb(alpha,p);
	r.linComb(-alpha,q);
	A.applyAdjOp(r,ATr);
	L.applyInvAdjOp(ATr,s);
	L.applyInvOp(s,g);
	
	atype newgamma = s.normsq();
	atype absbeta;
	if (ProtectedDivision<atype>(newgamma,gamma,absbeta)) {
	  RVLException e;
	  e << "Error: PCGNEStep_mod::run() from ProtectedDivision: beta=newgamma/gamma\n"
	    << "  newgamma = "<< newgamma <<"\n"
	    << "  gamma = "<< gamma <<"\n";
	  throw e;
	}
	Scalar beta = absbeta;

	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);

	xnorm = x.norm();
	L.applyOp(x,Lx);
	Lxnorm = Lx.norm();
	Vector<Scalar> tmp(x);
	tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
	xerrs = tmp.norm()/x0norm;
	L.applyOp(tmp,Lx);
	Lxerrs = Lx.norm()/Lxnorm;

      }
      catch (RVLException & e) {
	e<<"\ncalled from PCGNEStep_mod::run()\n";
	throw e;
      }
     
    }

    ~PCGNEStep_mod() {}
  };

  /** Single step of conjugate gradient iteration for the normal
      equations with Tikhonov regularization, i.e.,
        (A^T A + mu*I)x = A^T b
      with regularization parameter mu>=0.
  */
  template<typename Scalar>
  class RegCGNEStep_mod : public Algorithm {

    friend class CGNEAlg_mod<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    Vector<Scalar> const & x0;
    atype & mu;
    atype & rnorm;
    atype & nrnorm;
    atype & xnorm;
    atype & xerrs;
    
  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> g;    // gradient
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction
    atype gamma;         // ||g||^2
    atype x0norm;

  public:

    /**
     * Constructor
     */
    RegCGNEStep_mod(LinearOp<Scalar> const & _A,
		    Vector<Scalar> & _x,
		    Vector<Scalar> const & _b,
		    Vector<Scalar> const & _x0,
		    atype & _mu,
		    atype & _rnorm, 
		    atype & _nrnorm,
		    atype & _xnorm,
		    atype & _xerrs)
      : A(_A), x(_x), b(_b), x0(_x0), mu(_mu),
	rnorm(_rnorm), nrnorm(_nrnorm), 
	xnorm(_xnorm), xerrs(_xerrs), 
	r(A.getRange()), g(A.getDomain()),
	q(A.getRange()), p(A.getDomain()) { 

      // NOTE: initial x assumed to be zero vector
      r.copy(b);
      rnorm=r.norm();
      A.applyAdjOp(r,g);
      p.copy(g);
      nrnorm=g.norm();
      gamma = nrnorm*nrnorm;
      
      x0norm = x0.norm();      
      xnorm = x.norm();
      Vector<Scalar> tmp(x);
      tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
      atype tmp2=tmp.norm();
      if (ProtectedDivision<atype>(tmp2,x0norm,xerrs))
	x0norm=ScalarFieldTraits<Scalar>::One();
    }
      
    /**
     * Run a single step of CGNE
     */
    void run() {
      try {

	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype ptp = p.normsq();
	atype a = qtq + mu*ptp;
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,a,absalpha)) {
	  RVLException e;
	  e << "Error: RegCGNEStep_mod::run() from ProtectedDivision: alpha = gamma/qtq\n"
	    << "  gamma = "<< gamma <<"\n"
	    << "  qtq = "<< qtq <<"\n";
	  throw e;
	}
	Scalar alpha=absalpha;

	x.linComb(alpha,p);
	r.linComb(-alpha,q);
	A.applyAdjOp(r,g);
	g.linComb(-mu,x);

	atype newgamma = g.normsq();
	atype absbeta;
	if (ProtectedDivision<atype>(newgamma,gamma,absbeta)) {
	  RVLException e;
	  e << "Error: RegCGNEStep_mod::run() from ProtectedDivision: beta=newgamma/gamma\n"
	    << "  newgamma = "<< newgamma <<"\n"
	    << "  gamma = "<< gamma <<"\n";
	  throw e;
	}
	Scalar beta = absbeta;

	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);

	xnorm=x.norm();
	Vector<Scalar> tmp(x);
	tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
	xerrs=tmp.norm()/x0norm;	

      }
      catch (RVLException & e) {
	e<<"\ncalled from RegCGNEStep_mod::run()\n";
	throw e;
      }
     
    }

    ~RegCGNEStep_mod() {}
  };

  /** Preconditioned conjugate gradient iteration for the normal
      equations with Tikhonov regularization, i.e.,
        (A^T A + mu*I)x = A^T b
      with regularization parameter mu>=0.
  */
  template<typename Scalar>
  class RegPCGNEStep_mod : public Algorithm {

    friend class CGNEAlg_mod<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    LinearOpWithInverse<Scalar> const & L;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    Vector<Scalar> const & x0;
    atype & mu;
    atype & rnorm;
    atype & nrnorm;
    atype & xnorm;
    atype & Lxnorm;
    atype & xerrs;
    atype & Lxerrs;
    ostream & str;

  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> ATr;  // normal residual
    Vector<Scalar> s;    // s = LinvT * AT * r
    Vector<Scalar> g;    // gradient = preconditioned normal residual
                         // g = Linv * LinvT * AT * r
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction    
    Vector<Scalar> Lx;  
    atype gamma;         
    atype x0norm;
    atype Lx0norm;

  public:
  
    /**
     * Constructor
     */
    RegPCGNEStep_mod(LinearOp<Scalar> const & _A,
		     LinearOpWithInverse<Scalar> const & _L,
		     Vector<Scalar> & _x,
		     Vector<Scalar> const & _b,
		     Vector<Scalar> const & _x0,
		     atype & _mu,
		     atype & _rnorm, 
		     atype & _nrnorm,
		     atype & _xnorm,
		     atype & _Lxnorm,
		     atype & _xerrs,
		     atype & _Lxerrs,
		     ostream & _str)
      : A(_A), L(_L), x(_x), b(_b), x0(_x0), mu(_mu),
	rnorm(_rnorm), nrnorm(_nrnorm), 
	xnorm(_xnorm), Lxnorm(_Lxnorm),
	xerrs(_xerrs), Lxerrs(_Lxerrs),
	str(_str),
	r(A.getRange()),   ATr(A.getDomain()), 
	g(A.getDomain()),  s(A.getDomain()), 
	p(A.getDomain()),  q(A.getRange()),
	Lx(A.getDomain()) { 

      // sanity test, L should be a "square" operator compatible with A
      if ((L.getDomain() != A.getDomain()) ||
	  (L.getRange()  != A.getDomain())) {
	RVLException e;
	e << "Error RegPCGNEStep_mod constructor:\n"
	  << "Tikhonov operator L is not compatible with system operator A.\n";
	e << "Printing Tikhonov operator L:\n";
	L.write(e);
	e << "Printing system operator A:\n";
	A.write(e);
	throw e;
      }

      // NOTE: initial x assumed to be zero vector
      r.copy(b);
      rnorm=r.norm();
      A.applyAdjOp(r,ATr);
      L.applyInvAdjOp(ATr,s);
      L.applyInvOp(s,g);
      p.copy(g);
      gamma = s.normsq();
      nrnorm=sqrt(gamma);

      x0norm = x0.norm();
      L.applyOp(x0,Lx);
      Lx0norm = Lx.norm();

      xnorm = x.norm();
      L.applyOp(x,Lx);
      Lxnorm = Lx.norm();

      Vector<Scalar> tmp(x);
      tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
      atype tmp2 = tmp.norm();
      if (ProtectedDivision<atype>(tmp2,x0norm,xerrs))
	x0norm=ScalarFieldTraits<Scalar>::One();

      L.applyOp(tmp,Lx);
      tmp2 = Lx.norm();
      if (ProtectedDivision<atype>(tmp2,Lx0norm,Lxerrs))
	Lx0norm=ScalarFieldTraits<Scalar>::One();
    }
      
    /**
     *  Run a single step of preconditioned CGNE
     */
    void run() {
      try {

	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype ptp = p.normsq();
	atype a = qtq + mu*ptp;
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,a,absalpha)) {
	  RVLException e;
	  e << "Error: RegPCGNEStep_mod::run() from ProtectedDivision: alpha = gamma/qtq\n"
	    << "  gamma = "<< gamma <<"\n"
	    << "  qtq = "<< qtq <<"\n";
	  throw e;
	}
	Scalar alpha=absalpha;

	x.linComb(alpha,p);
	r.linComb(-alpha,q);
	A.applyAdjOp(r,ATr);
	ATr.linComb(-mu,x);
	L.applyInvAdjOp(ATr,s);
	L.applyInvOp(s,g);
	
	atype newgamma = s.normsq();
	atype absbeta;
	if (ProtectedDivision<atype>(newgamma,gamma,absbeta)) {
	  RVLException e;
	  e << "Error: RegPCGNEStep_mod::run() from ProtectedDivision: beta=newgamma/gamma\n"
	    << "  newgamma = "<< newgamma <<"\n"
	    << "  gamma = "<< gamma <<"\n";
	  throw e;
	}
	Scalar beta = absbeta;

	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);

	xnorm = x.norm();
	L.applyOp(x,Lx);
	Lxnorm = Lx.norm();
	Vector<Scalar> tmp(x);
	tmp.linComb(-ScalarFieldTraits<Scalar>::One(),x0);
	xerrs = tmp.norm()/x0norm;
	L.applyOp(tmp,Lx);
	Lxerrs = Lx.norm()/Lxnorm;

      }
      catch (RVLException & e) {
	e<<"\ncalled from RegPCGNEStep_mod::run()\n";
	throw e;
      }
     
    }

    ~RegPCGNEStep_mod() {}
  };
  
  /** Conjugate gradient algorithm - efficient implementation for
      normal equations
      \f[ A^{\prime} A x = A^{\prime} b\f]
      for solving the linear least squares problem
      \f[ \min_{x} \vert A x - b \vert^2 \f].

      See cgnealg.hh for more info.
  */
  template<typename Scalar>
  class CGNEAlg_mod: public Algorithm, public Terminator {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:

    atype & rnorm;   // residual norm
    atype & nrnorm;  // gradient norm
    atype & xnorm;
    atype & Lxnorm;
    atype & xerrs;
    atype & Lxerrs;
    atype rtol;       // tolerance residual norm
    atype nrtol;      // tolerance gradient norm 
    atype xtol;
    atype Lxtol;                   
    atype maxstep;     // upper bound for net step x-x0
    int maxcount;      // upper bound for iteration count
    int count;         // actual iteration count
    mutable bool proj; // whether step is projected onto TR boundary
    ostream & str;     // stream for report output
    Algorithm * step;  // CGNE or PCGNE step

  public:

    /** Constructor - basic algorithm
	See cgnealg.hh for more info.
    */
    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		RVL::LinearOp<Scalar> const & A, 
		RVL::Vector<Scalar> const & b, 
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _xerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), nrnorm(_nrnorm),
	xnorm(_xnorm), Lxnorm(_xnorm),
	xerrs(_xerrs), Lxerrs(_xerrs),
	rtol(_rtol),   nrtol(_nrtol), 
	xtol(_xtol),   Lxtol(_xtol), 
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    {
      try {
	x.zero();
	Vector<Scalar> x0(A.getDomain(),true);
	step = new CGNEStep_mod<Scalar>(A,x,b,x0,
					rnorm,nrnorm,
					xnorm,xerrs);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (not preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		RVL::LinearOp<Scalar> const & A, 
		RVL::Vector<Scalar> const & b, 
		RVL::Vector<Scalar> const & x0, 
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _xerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), nrnorm(_nrnorm),
	xnorm(_xnorm), Lxnorm(_xnorm),
	xerrs(_xerrs), Lxerrs(_xerrs),
	rtol(_rtol),   nrtol(_nrtol), 
	xtol(_xtol),   Lxtol(_xtol), 
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    {
      try {
	x.zero();
	step = new CGNEStep_mod<Scalar>(A,x,b,x0,
					rnorm,nrnorm,
					xnorm,xerrs);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (not preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		RVL::LinearOp<Scalar> const & A, 
		RVL::Vector<Scalar> const & b, 
		atype &_mu,
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _xerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), nrnorm(_nrnorm),
	xnorm(_xnorm), Lxnorm(_xnorm),
	xerrs(_xerrs), Lxerrs(_xerrs),
	rtol(_rtol),   nrtol(_nrtol), 
	xtol(_xtol),   Lxtol(_xtol), 
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    {
      try {
	x.zero();
	Vector<Scalar> x0(A.getDomain(),true);
	step = new RegCGNEStep_mod<Scalar>(A,x,b,x0,_mu,
					   rnorm,nrnorm,
					   xnorm,xerrs);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (not preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		RVL::LinearOp<Scalar> const & A, 
		RVL::Vector<Scalar> const & b, 
		RVL::Vector<Scalar> const & x0,
		atype & _mu,
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _xerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), nrnorm(_nrnorm),
	xnorm(_xnorm), Lxnorm(_xnorm),
	xerrs(_xerrs), Lxerrs(_xerrs),
	rtol(_rtol),   nrtol(_nrtol), 
	xtol(_xtol),   Lxtol(_xtol), 
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    {
      try {
	x.zero();
	step = new RegCGNEStep_mod<Scalar>(A,x,b,x0,_mu,
					   rnorm,nrnorm,
					   xnorm,xerrs);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (not preconditioned)\n";
	throw e;
      }
    }

    /** Constructor - preconditioned algorithm
	See cgneal.hh for more info.
    */
    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		LinearOp<Scalar> const & A, 
		LinearOpWithInverse<Scalar> const & L, 
		Vector<Scalar> const & b, 
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _Lxnorm,
		atype & _xerrs,
		atype & _Lxerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		atype _Lxtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), 	nrnorm(_nrnorm), 
	xnorm(_xnorm),	Lxnorm(_Lxnorm),
	xerrs(_xerrs),  Lxerrs(_Lxerrs),
	rtol(_rtol), 	nrtol(_nrtol), 
	xtol(_xtol),	Lxtol(_Lxtol),
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    { 
      try {
	x.zero();
	Vector<Scalar> x0(A.getDomain(),true);
	step = new PCGNEStep_mod<Scalar>(A,L,x,b,x0,
					 rnorm,nrnorm,
					 xnorm,Lxnorm,
					 xerrs,Lxerrs,
					 str);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		LinearOp<Scalar> const & A, 
		LinearOpWithInverse<Scalar> const & L, 
		Vector<Scalar> const & b, 
		Vector<Scalar> const & x0, 
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _Lxnorm,
		atype & _xerrs,
		atype & _Lxerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		atype _Lxtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), 	nrnorm(_nrnorm), 
	xnorm(_xnorm),	Lxnorm(_Lxnorm),
	xerrs(_xerrs),  Lxerrs(_Lxerrs),
	rtol(_rtol), 	nrtol(_nrtol), 
	xtol(_xtol),	Lxtol(_Lxtol),
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    { 
      try {
	x.zero();
	step = new PCGNEStep_mod<Scalar>(A,L,x,b,x0,
					 rnorm,nrnorm,
					 xnorm,Lxnorm,
					 xerrs,Lxerrs,
					 str);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		LinearOp<Scalar> const & A, 
		LinearOpWithInverse<Scalar> const & L, 
		Vector<Scalar> const & b, 
		atype & _mu,
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _Lxnorm,
		atype & _xerrs,
		atype & _Lxerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		atype _Lxtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), 	nrnorm(_nrnorm), 
	xnorm(_xnorm),	Lxnorm(_Lxnorm),
	xerrs(_xerrs),  Lxerrs(_Lxerrs),
	rtol(_rtol), 	nrtol(_nrtol), 
	xtol(_xtol),	Lxtol(_Lxtol),
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    { 
      try {
	x.zero();
	Vector<Scalar> x0(A.getDomain(),true);
	step = new RegPCGNEStep_mod<Scalar>(A,L,x,b,x0,_mu,
					    rnorm,nrnorm,
					    xnorm,Lxnorm,
					    xerrs,Lxerrs,
					    str);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (preconditioned)\n";
	throw e;
      }
    }

    CGNEAlg_mod(RVL::Vector<Scalar> & x, 
		LinearOp<Scalar> const & A, 
		LinearOpWithInverse<Scalar> const & L, 
		Vector<Scalar> const & b, 
		Vector<Scalar> const & x0,
		atype & _mu,
		atype & _rnorm,
		atype & _nrnorm,
		atype & _xnorm,
		atype & _Lxnorm,
		atype & _xerrs,
		atype & _Lxerrs,
		atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
		atype _xtol = numeric_limits<atype>::max(),
		atype _Lxtol = numeric_limits<atype>::max(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : rnorm(_rnorm), 	nrnorm(_nrnorm), 
	xnorm(_xnorm),	Lxnorm(_Lxnorm),
	xerrs(_xerrs),  Lxerrs(_Lxerrs),
	rtol(_rtol), 	nrtol(_nrtol), 
	xtol(_xtol),	Lxtol(_Lxtol),
	maxstep(_maxstep), 
	maxcount(_maxcount), 
	count(0), 
	proj(false), 
	str(_str)
    { 
      try {
	x.zero();
	step = new RegPCGNEStep_mod<Scalar>(A,L,x,b,x0,_mu,
					    rnorm,nrnorm,
					    xnorm,Lxnorm,
					    xerrs,Lxerrs,
					    str);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod constructor (preconditioned)\n";
	throw e;
      }
    }

    ~CGNEAlg_mod() {
      delete step;
    }

    bool query() { return proj; }

    void run() { 
      try {
	// access to internals
	CGNEStep_mod<Scalar>  * cg  = dynamic_cast<CGNEStep_mod<Scalar> *>(step);
	PCGNEStep_mod<Scalar> * pcg = dynamic_cast<PCGNEStep_mod<Scalar> *>(step);
	RegCGNEStep_mod<Scalar> * rcg = dynamic_cast<RegCGNEStep_mod<Scalar> *>(step);
	RegPCGNEStep_mod<Scalar> * rpcg = dynamic_cast<RegPCGNEStep_mod<Scalar> *>(step);

	if ((!cg && !pcg && !rcg && !rpcg) || (cg && pcg && rcg && rpcg)) {
	  RVLException e;
	  e<<"Error: CGNEAlg_mod::run\n";
	  e<<"  unable to determine whether step data member\n";
	  e<<"  is preconditioned or not. PANIC!\n";
	  throw e;
	}
	
	// terminator for CGNE iteration
	vector<string> names(6);
	vector<atype *> nums(6);
	vector<atype> tols(6);
	atype rnorm0=rnorm;
	atype nrnorm0=nrnorm;
          
	names[0]="Residual Norm";          nums[0]=&rnorm;  tols[0]=rtol*rnorm0;
	names[1]="Gradient Norm";          nums[1]=&nrnorm; tols[1]=nrtol*nrnorm0;
	names[2]="Sol Norm |x|";           nums[2]=&xnorm;  tols[2]=xtol;
	names[3]="Sol Norm |Lx|";          nums[3]=&Lxnorm; tols[3]=Lxtol;
	names[4]="Sol Err Norm |x-x0|";    nums[4]=&xerrs;  tols[4]=xtol;
	names[5]="Sol Err Norm |L(x-x0)|"; nums[5]=&Lxerrs; tols[5]=Lxtol;

	if (pcg||rpcg) str<<"========================== BEGIN PCGNE =========================\n";
	else 	str<<"========================== BEGIN CGNE =========================\n";
	VectorCountingThresholdIterationTable<atype> stop1(maxcount+1,names,nums,tols,str);
	stop1.init();

	Terminator * stop2 = NULL;
	if (cg)  stop2 = new BallProjTerminator<Scalar>( cg->x,maxstep,str);  
	if (pcg) stop2 = new BallProjTerminator<Scalar>(pcg->x,maxstep,str); 
	if (rcg)  stop2 = new BallProjTerminator<Scalar>( rcg->x,maxstep,str);  
	if (rpcg) stop2 = new BallProjTerminator<Scalar>(rpcg->x,maxstep,str); 
	
	// terminate if either
	OrTerminator stop(stop1,*stop2);
	// loop
	LoopAlg doit(*step,stop);
	doit.run();

	// must recompute residual if scaling occured 
	proj = stop2->query();
	if (proj) {
	  if (cg||rcg) {
	    Vector<Scalar> temp((cg->A).getRange());
	    (cg->A).applyOp((cg->x),temp);
	    temp.linComb(-1.0,(cg->b));
	    rnorm=temp.norm();
	    Vector<Scalar> temp1((cg->A).getDomain());
	    (cg->A).applyAdjOp(temp,temp1);
	    nrnorm=temp1.norm();
	  }
	  if (pcg||rpcg) {
	    Vector<Scalar> temp((pcg->A).getRange());
	    (pcg->A).applyOp((pcg->x),temp);
	    temp.linComb(-1.0,(pcg->b));
	    rnorm=temp.norm();
	    Vector<Scalar> temp1((pcg->A).getDomain());
	    Vector<Scalar> temp2((pcg->A).getDomain());
	    Vector<Scalar> temp3((pcg->A).getDomain());
	    (pcg->A).applyAdjOp(temp,temp1);
	    (pcg->L).applyInvAdjOp(temp1,temp2);
	    (pcg->L).applyInvOp(temp2,temp3);
	    atype tmpgamma = abs(temp1.inner(temp3));
	    nrnorm=sqrt(tmpgamma);
	  }
	}

	count = stop1.getCount();
	if (pcg||rpcg) str<<"=========================== END PCGNE ==========================\n";
	else str<<"=========================== END CGNE ==========================\n";
	// display results
	str<<"\n ****************** CGNE summary *****************  \n"
	   <<"Initial residual norm      = "<< rnorm0 <<"\n"
	   <<"Residual norm              = "<< rnorm  <<"\n"
	   <<"Residual redn              = "<< rnorm/rnorm0 <<"\n"
	   <<"Initial gradient norm      = "<< nrnorm0 <<"\n"
	   <<"Gradient norm              = "<< nrnorm  <<"\n"
	   <<"Gradient redn              = "<< nrnorm/nrnorm0 <<"\n"
	   <<"Sol norm |x|               = "<< xnorm << "\n"
	   <<"Sol norm |Lx|              = "<< Lxnorm << "\n"
	   <<"Sol err norm |x-x0|        = "<< xnorm << "\n"
	   <<"Sol err norm |L(x-x0)|     = "<< Lxnorm << "\n";	
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg_mod::run\n";
	throw e;
      }
    }

    int getCount() const { return count; }


    // disable default, copy constructors
    CGNEAlg_mod();
    CGNEAlg_mod(CGNEAlg_mod<Scalar> const &);

  };

  /** 
   * Data class for CGNE policy
   */
  template<typename Scalar>
  class CGNEPolicyData {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
  
  public:

    atype rtol;
    atype nrtol;
    atype xtol;
    atype Lxtol;
    atype Delta;
    int maxcount;
    bool verbose;

    CGNEPolicyData(atype _rtol = numeric_limits<atype>::max(),
		   atype _nrtol = numeric_limits<atype>::max(),
		   atype _xtol = numeric_limits<atype>::max(),
		   atype _Lxtol = numeric_limits<atype>::max(),
		   atype _Delta = numeric_limits<atype>::max(),
		   int _maxcount = 0,
		   bool _verbose = false)
      : rtol(_rtol), nrtol(_nrtol), xtol(_xtol), Lxtol(_Lxtol),
	Delta(_Delta), maxcount(_maxcount), verbose(_verbose) {}
      
    CGNEPolicyData(CGNEPolicyData<Scalar> const & a) 
      : rtol(a.rtol), nrtol(a.nrtol), xtol(a.xtol), Lxtol(a.Lxtol),
	Delta(a.Delta), maxcount(a.maxcount), verbose(a.verbose) {}

    ostream & write(ostream & str) const {
      str<<"\n==============================================\n"
	 <<"CGNEPolicyData: \n"
	 <<"rtol      = "<< rtol     <<"\n"
	 <<"nrtol     = "<< nrtol    <<"\n"
	 <<"xtol      = "<< xtol     <<"\n"
	 <<"Lxtol     = "<< Lxtol    <<"\n"
	 <<"Delta     = "<< Delta    <<"\n"
	 <<"maxcount  = "<< maxcount <<"\n"
	 <<"verbose   = "<< verbose  <<"\n"
	 <<"\n==============================================\n";
      return str;
    }
  };

  /** policy class for creation of CGNEAlg in trust region solver and 
      any other algorithm needing a least squares solver component - build
      method creates CGNEAlg with these attributes:

      rtol     = residual threshhold for convergence
      nrtol    = normal residual (LS gradient) tolerance for convergence
      Delta    = trust radius - truncate iteration when reached
      maxcount = max number of iterations permitted

      Default values set to cause immediate return from CGNEAlg::run.

      Other attributes are arguments of build method.

      Conforms to specifications described in TRGNAlg docs.

      Usage: use as policy type, i.e. class name as template parameter
      to TRGNAlg constructor. After construction of TRGNAlg, which IS
      a CGNEPolicy by inheritance, call the assign method on it to set
      rtol, nrtol, and maxcount, BEFORE calling TRGNAlg::run() - else
      you will get immediate termination, as intended.
  */

  template<typename Scalar> 
  class CGNEPolicy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:
    mutable atype rtol;
    mutable atype nrtol;
    mutable atype xtol;    
    mutable atype Lxtol;
    mutable int maxcount;
    mutable bool verbose;
    mutable std::ostream nullstr;

  public:
    /** only Delta need be changed repeatedly, as opposed
	to set post-construction. Simplest way to do this - make
	it public
    */
    mutable atype Delta;

    /**
       build method - see TRGNAlg specs

       @param x - solution vector, initialize to zero on input,
       estimated solution on output 

       @param A - Linear Operator of least squares problem

       @param d - data vector of least squares problem

       @param rnorm - reference to residual norm scalar, norm of RHS
       on input, of residual on output

       @param nrnorm - reference to normal residual norm scalar, norm
       of normal residual (least squares gradient) on input, updated
       to estimated solution on output
       
       @param str - verbose output stream
    */
    CGNEAlg_mod<Scalar> * build(Vector<Scalar> & x, 
				LinearOp<Scalar> const & A,
				Vector<Scalar> const & b,
				atype & rnorm,
				atype & nrnorm,
				atype & xnorm,
				atype & xerrs,
				ostream & str) const {
      try {
	if (verbose) 
	  return new CGNEAlg_mod<Scalar>(x,A,b,
					 rnorm,nrnorm,
					 xnorm,xerrs,
					 rtol,nrtol,xtol,
					 maxcount,Delta,str);
	else
	  return new CGNEAlg_mod<Scalar>(x,A,b,
					 rnorm,nrnorm,
					 xnorm,xerrs,
					 rtol,nrtol,xtol,
					 maxcount,Delta,nullstr);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEPolicy::build\n";
	e<<"inputs: \n";
	e<<"**** x:\n";
	x.write(e);
	e<<"**** A:\n";
	A.write(e);
	e<<"**** b:\n";
	b.write(e);
	throw e;
      }
    }
      
    CGNEAlg_mod<Scalar> * build(Vector<Scalar> & x,
				LinearOp<Scalar> const & A,
				LinearOpWithInverse<Scalar> const & L,
				Vector<Scalar> const & b,
				atype & rnorm,
				atype & nrnorm,
				atype & xnorm,
				atype & Lxnorm,
				atype & xerrs,
				atype & Lxerrs,
				ostream & str) const {
      try {
        if (verbose)
	  return new CGNEAlg_mod<Scalar>(x,A,L,b,
					 rnorm,nrnorm,
					 xnorm,Lxnorm,
					 xerrs,Lxerrs,
					 rtol,nrtol,xtol,Lxtol,
					 maxcount,Delta,str);
        else
	  return new CGNEAlg_mod<Scalar>(x,A,L,b,
					 rnorm,nrnorm,
					 xnorm,Lxnorm,
					 xerrs,Lxerrs,
					 rtol,nrtol,xtol,Lxtol,
					 maxcount,Delta,nullstr);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEPolicy::build\n";
	e<<"inputs: \n";
	e<<"**** x:\n";
	x.write(e);
	e<<"**** A:\n";
	A.write(e);
	e<<"**** b:\n";
	b.write(e);
	throw e;
      }
    }
      
    /** post-construction initialization
	@param _rtol - residual norm stopping threshhold
	@param _nrtol - normal residual (LS gradient) norm stopping threshhold
	@param _maxcount - max number of permitted iterations
    */
    void assign(atype _rtol, atype _nrtol, atype _xtol, atype _Lxtol,
		atype _Delta, int _maxcount, bool _verbose) {
      rtol=_rtol; nrtol=_nrtol; xtol=_xtol; Lxtol=_Lxtol;
      Delta=_Delta; maxcount=_maxcount; verbose=_verbose;
    }

    /** parameter table overload */
    void assign(Table const & t) {
      rtol=getValueFromTable<atype>(t,"CGNE_ResTol");
      nrtol=getValueFromTable<atype>(t,"CGNE_GradTol");
      xtol=getValueFromTable<atype>(t,"CGNE_xTol");  
      Lxtol=getValueFromTable<atype>(t,"CGNE_LxTol"); 
      Delta=getValueFromTable<atype>(t,"TR_Delta");
      maxcount=getValueFromTable<int>(t,"CGNE_MaxItn"); 
      verbose=getValueFromTable<bool>(t,"CGNE_Verbose");
    }

    /** data struct overload */
    void assign(CGNEPolicyData<Scalar> const & s) {
      rtol=s.rtol;
      nrtol=s.nrtol;
      xtol=s.xtol;
      Lxtol=s.Lxtol;
      Delta=s.Delta;
      maxcount=s.maxcount;
      verbose=s.verbose;
    }

    /** main constructor - acts as default. Default values of
	parameters set to result in immediate return, no
	iteration. Note that policy design requires that default
	construction must be valid, and all run-time instance data be
	initiated post-construction, in this case by the assign
	function, to be called by drivers of user classes (subclassed
	from this and with this as template param).
    */
    CGNEPolicy(atype _rtol = numeric_limits<atype>::max(),
	       atype _nrtol = numeric_limits<atype>::max(),
	       atype _xtol = numeric_limits<atype>::max(),
	       atype _Lxtol = numeric_limits<atype>::max(),
	       atype _Delta = numeric_limits<atype>::max(),
	       int _maxcount = 0,
	       bool _verbose = true)
      : rtol(_rtol), nrtol(_nrtol), xtol(_xtol), Lxtol(_Lxtol),
	Delta(_Delta), maxcount(_maxcount), verbose(_verbose), nullstr(0) {}

    CGNEPolicy(CGNEPolicy<Scalar> const & p)
      : rtol(p.rtol), nrtol(p.nrtol), xtol(p.xtol), Lxtol(p.Lxtol),
	Delta(p.Delta), maxcount(p.maxcount), verbose(p.verbose), nullstr(0) {}
  };
}

#endif









